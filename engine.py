"""
PlasmidGuard core engine.

Fixes applied over the original blueprint:
  1. LOCAL alignment (Smith-Waterman via PairwiseAligner mode='local'),
     not global — a short Sanger read should not be forced to span
     the whole plasmid.
  2. Strand check — tries the read and its reverse complement, keeps
     whichever aligns better, before calling anything a mutation.
  3. Quality trimming for .ab1 reads — strips low-confidence bases at
     both ends using Phred scores before alignment.
  4. Circularity handling for the design — the design is extended with
     a wraparound copy of its own prefix so origin-spanning reads still
     seed and align correctly; mutation positions from the wrapped
     region are folded back into real plasmid coordinates afterward.
  5. Seed-and-extend alignment — a k-mer index locates the approximate
     match location in O(design_length), then exact DP alignment runs
     only on a small window around it. This is what makes runtime
     independent of design size (a 500kb design behaves like a 5kb one);
     the remaining cost driver is read length, not design length.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import io

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


# ---------- Step 1: Parsing ----------

def _clean_parser_error(e: ValueError, file_kind: str) -> str:
    """Biopython's parser errors are often multi-paragraph explanations
    aimed at library users, not end users. Collapse anything long/technical
    into a short, actionable message instead of dumping it verbatim."""
    msg = str(e).strip()
    if len(msg) > 150:
        return f"Couldn't parse this as a valid {file_kind} file — check it isn't corrupted or in an unexpected format."
    return f"Couldn't read the {file_kind} file: {msg.splitlines()[0]}"


@dataclass
class DesignRecord:
    seq: str
    features: list  # list of dicts: {type, name, start, end, strand}
    is_circular: bool = True


@dataclass
class RealityRead:
    seq: str
    qualities: Optional[list] = None  # per-base Phred scores, or None for fasta


# Length sanity caps. With seed-and-extend alignment (see align_best_orientation
# below), cost no longer scales with design length. The remaining cost driver is
# the windowed exact-DP step, roughly O(read_length^2): measured ~1.6s at a
# 6,000bp read, but ~5.5s at 8,000bp, so MAX_READ_LENGTH is 6,000 (real single
# Sanger/NGS reads are almost always under 1.5kb).
#
# MAX_DESIGN_LENGTH is 100,000 because real plasmids are 2-20kb and this tool may
# be run as a public hosted demo: measured peak memory per analysis is ~25MB at
# a 100kb design but ~214MB at 500kb, which adds up quickly under public traffic.
MAX_DESIGN_LENGTH = 100_000
MAX_READ_LENGTH = 6_000


class PlasmidGuardError(ValueError):
    """A garbage/invalid-input condition we detected and want to surface
    to the user with a clear message, as opposed to a raw parser traceback."""


def parse_design(file_bytes: bytes, filename: str) -> DesignRecord:
    """Parse a .gb/.gbk or .fasta design file."""
    if not file_bytes.strip():
        raise PlasmidGuardError("The design file is empty.")

    fmt = "genbank" if filename.lower().endswith((".gb", ".gbk")) else "fasta"
    try:
        text = file_bytes.decode("utf-8")
    except UnicodeDecodeError:
        raise PlasmidGuardError(
            "This file isn't valid UTF-8 text — it may have been saved with a different encoding, "
            "or it's a binary file that doesn't belong in this slot. Re-save it as plain UTF-8 text."
        )

    if fmt == "fasta" and text.lstrip().startswith("LOCUS"):
        raise PlasmidGuardError(
            "This looks like a GenBank file but doesn't have a .gb/.gbk extension — "
            "rename it (or re-export it) so we parse it correctly."
        )

    handle = io.StringIO(text)
    try:
        record = SeqIO.read(handle, fmt)
    except ValueError as e:
        if "more than one record" in str(e).lower():
            raise PlasmidGuardError(
                "The design file has more than one sequence in it. Upload a file with a single record."
            ) from e
        raise PlasmidGuardError(_clean_parser_error(e, "design")) from e

    seq = str(record.seq).upper()
    if not seq or set(seq) - set("ACGTN"):
        bad_chars = sorted(set(seq) - set("ACGTN"))
        raise PlasmidGuardError(
            f"The design sequence contains non-DNA characters ({', '.join(bad_chars[:5])}) — "
            f"is this actually a DNA file?"
        )
    if len(seq) > MAX_DESIGN_LENGTH:
        raise PlasmidGuardError(
            f"This design is {len(seq):,} bp — larger than the {MAX_DESIGN_LENGTH:,} bp this tool "
            f"is built for (single plasmids, not genomes/chromosomes). Trim it to the relevant region."
        )

    features = []
    is_circular = True  # sane default: most plasmid designs are circular
    if fmt == "genbank":
        topology = record.annotations.get("topology", "").lower()
        if topology == "linear":
            is_circular = False
        for f in record.features:
            if f.type in ("CDS", "gene", "promoter", "misc_feature"):
                if f.location is None:
                    # Biopython sets this to None for features it couldn't
                    # resolve — e.g. an annotation like "8569..276" that
                    # implies wrapping the origin, on a record whose
                    # topology is linear/undefined so it can't safely
                    # compute what that wrap means. Real Addgene/SnapGene
                    # exports hit this. Skip rather than crash on .start.
                    continue
                start, end = int(f.location.start), int(f.location.end)
                if start >= end or start < 0 or end > len(seq):
                    # Malformed coordinates (start>=end, negative, or past the
                    # end of the sequence) — skip rather than let it silently
                    # corrupt frame/codon math downstream.
                    continue
                features.append({
                    "type": f.type,
                    "name": (
                        f.qualifiers.get("gene", f.qualifiers.get(
                            "label", f.qualifiers.get(
                                "locus_tag", f.qualifiers.get("product", [f.type]))))[0]
                    ),
                    "start": start,
                    "end": end,
                    "strand": f.location.strand or 1,
                    "compound_location": len(f.location.parts) > 1,
                })
    return DesignRecord(seq=seq, features=features, is_circular=is_circular)


def parse_reality(file_bytes: bytes, filename: str) -> RealityRead:
    """Parse a .ab1 (Sanger trace), .fasta, or .fastq reality file."""
    if not file_bytes.strip():
        raise PlasmidGuardError("The sequencing file is empty.")

    name = filename.lower()
    try:
        if name.endswith(".ab1"):
            handle = io.BytesIO(file_bytes)
            record = SeqIO.read(handle, "abi")
            quals = record.letter_annotations.get("phred_quality")
        elif name.endswith((".fastq", ".fq")):
            text = file_bytes.decode("utf-8")
            record = SeqIO.read(io.StringIO(text), "fastq")
            quals = record.letter_annotations.get("phred_quality")
        else:
            text = file_bytes.decode("utf-8")
            if text.lstrip().startswith("LOCUS"):
                raise PlasmidGuardError(
                    "This looks like a GenBank file uploaded as the sequencing result — "
                    "that slot expects the raw read (.ab1/.fasta/.fastq)."
                )
            record = SeqIO.read(io.StringIO(text), "fasta")
            quals = None
    except ValueError as e:
        if isinstance(e, PlasmidGuardError):
            raise
        if isinstance(e, UnicodeDecodeError):
            raise PlasmidGuardError(
                "This file isn't valid UTF-8 text — it may have been saved with a different encoding, "
                "or it's a binary file that doesn't belong in this slot. Re-save it as plain UTF-8 text."
            ) from e
        if "more than one record" in str(e).lower():
            raise PlasmidGuardError(
                "The sequencing file has more than one read in it. Upload a single read at a time."
            ) from e
        raise PlasmidGuardError(_clean_parser_error(e, "sequencing")) from e

    seq = str(record.seq).upper()
    if not seq.replace("N", ""):
        raise PlasmidGuardError(
            "This sequencing read is entirely 'N' (no confident base calls) — "
            "the sequencing run likely failed. There's nothing to compare against the design."
        )
    if len(seq) > MAX_READ_LENGTH:
        raise PlasmidGuardError(
            f"This read is {len(seq):,} bp — longer than the {MAX_READ_LENGTH:,} bp expected for a "
            f"single sequencing read. Check you uploaded a read and not a full reference/assembly."
        )
    return RealityRead(seq=seq, qualities=quals)



# ---------- Step 2: Quality trimming ----------

def quality_trim(read: RealityRead, threshold: int = 20, window: int = 10) -> RealityRead:
    """Trim low-confidence bases off both ends. No-op if no quality scores."""
    if not read.qualities:
        return read
    q = read.qualities
    n = len(q)

    def first_good_index():
        for i in range(n - window):
            if sum(q[i:i + window]) / window >= threshold:
                return i
        return 0

    def last_good_index():
        for i in range(n, window, -1):
            if sum(q[i - window:i]) / window >= threshold:
                return i
        return n

    start = first_good_index()
    end = last_good_index()
    if end <= start:
        start, end = 0, n  # give up trimming rather than return empty
    return RealityRead(seq=read.seq[start:end], qualities=q[start:end])


# ---------- Step 3: Strand-aware local alignment ----------

def _make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -3
    aligner.extend_gap_score = -1
    return aligner


def _reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


@dataclass
class AlignmentResult:
    score: float
    design_seq: str          # design sequence actually used (may be doubled)
    read_seq: str             # read sequence actually used (may be revcomp)
    read_qualities: Optional[list]  # quality scores aligned to read_seq order
    strand: int                # 1 or -1 relative to design
    blocks: list               # list of (d_start, d_end, r_start, r_end) gap-free blocks
    used_doubled_design: bool
    coverage: float = 0.0       # fraction of the read that participated in the alignment
    identity: float = 0.0       # fraction of aligned (non-gap) bases that match exactly


def _score_blocks(blocks: list, design_seq: str, read_seq: str, read_len: int) -> tuple[float, float]:
    if not blocks:
        return 0.0, 0.0
    aligned_read_bases = sum(r1 - r0 for (_, _, r0, r1) in blocks)
    matches = sum(
        1
        for (d0, d1, r0, r1) in blocks
        for offset in range(d1 - d0)
        if design_seq[d0 + offset] == read_seq[r0 + offset]
    )
    coverage = aligned_read_bases / read_len if read_len else 0.0
    identity = matches / aligned_read_bases if aligned_read_bases else 0.0
    return coverage, identity


# Below these thresholds we don't trust the alignment enough to report
# mutations from it — a real Sanger read against its correct design should
# clear both easily; garbage/unrelated/wrong-file inputs typically won't.
MIN_COVERAGE = 0.5
MIN_IDENTITY = 0.75

# k-mer size for seeding. Long enough that a random 4-letter-alphabet match
# is vanishingly unlikely by chance (4^16 ≈ 4.3 billion), short enough that
# real reads with a few percent sequencing error still share plenty of
# exact k-mers with the correct design region.
SEED_K = 16


def _build_kmer_index(seq: str, k: int) -> dict:
    index: dict = {}
    for i in range(len(seq) - k + 1):
        index.setdefault(seq[i:i + k], []).append(i)
    return index


def _best_seed_diagonal(design_index: dict, query: str, k: int) -> Optional[tuple[int, int]]:
    """Find the (design_pos - query_pos) diagonal with the most shared
    k-mers between query and design. Returns (diagonal, seed_count) or
    None if the query shares no k-mer with the design at all — which
    means this orientation almost certainly doesn't belong to this design."""
    if len(query) < k:
        return None
    diag_counts: dict = {}
    for i in range(len(query) - k + 1):
        for d_pos in design_index.get(query[i:i + k], ()):
            diag = d_pos - i
            diag_counts[diag] = diag_counts.get(diag, 0) + 1
    if not diag_counts:
        return None
    best_diag = max(diag_counts, key=diag_counts.get)
    return best_diag, diag_counts[best_diag]


def align_best_orientation(design: DesignRecord, read: RealityRead) -> AlignmentResult:
    """Seed-and-extend alignment: find the approximate match location via
    shared k-mers (cheap, O(design_length)), then run exact DP alignment
    only on a small window around that location (cheap, independent of
    design length). This is what lets the tool handle large designs
    without the runtime scaling with design size — the earlier full-design
    DP approach was O(design_length x read_length) per orientation, which
    became a multi-second-to-multi-minute hang on anything beyond a small
    plasmid.

    Circularity is handled by extending the design with a wraparound copy
    of its own prefix before indexing/seeding, so a read spanning the
    origin still finds seeds and gets a contiguous window to align against.
    Mutation positions from the wrapped region are folded back into real
    plasmid coordinates by the caller (run_pipeline), not here — this
    function's `design_seq` field may report positions >= len(design.seq).
    """
    read_len = len(read.seq)
    wrap_pad = min(len(design.seq), max(read_len, 200)) if design.is_circular else 0
    design_ext = design.seq + design.seq[:wrap_pad]
    design_index = _build_kmer_index(design_ext, SEED_K)

    def orientations():
        yield 1, read.seq, read.qualities
        rc_quals = list(reversed(read.qualities)) if read.qualities else None
        yield -1, _reverse_complement(read.seq), rc_quals

    aligner = _make_aligner()
    candidates = []

    for strand, seq, quals in orientations():
        seed = _best_seed_diagonal(design_index, seq, SEED_K)
        if seed is None:
            candidates.append((0, float("-inf"), strand, seq, quals, None, 0.0, 0.0))
            continue
        diagonal, seed_count = seed

        # Window around the seeded location, padded generously to absorb
        # indels without needing to re-seed.
        padding = max(50, int(0.25 * len(seq)))
        window_start = max(0, diagonal - padding)
        window_end = min(len(design_ext), diagonal + len(seq) + padding)
        design_window = design_ext[window_start:window_end]

        alns = aligner.align(design_window, seq)
        try:
            top = next(iter(alns))
        except StopIteration:
            candidates.append((seed_count, float("-inf"), strand, seq, quals, None, 0.0, 0.0))
            continue

        blocks = [(int(d0) + window_start, int(d1) + window_start, int(r0), int(r1))
                  for (d0, d1), (r0, r1) in zip(top.aligned[0], top.aligned[1])]
        coverage, identity = _score_blocks(blocks, design_ext, seq, read_len)
        candidates.append((seed_count, top.score, strand, seq, quals, blocks, coverage, identity))

    # Prefer more seed support first (a stronger signal that this is the
    # right orientation/location at all), then alignment score as tiebreak.
    candidates.sort(key=lambda c: (c[0], c[1]), reverse=True)
    seed_count, score, strand, seq, quals, blocks, coverage, identity = candidates[0]

    if blocks is None:
        return AlignmentResult(
            score=0.0, design_seq=design_ext, read_seq=seq, read_qualities=quals,
            strand=strand, blocks=[], used_doubled_design=False,
            coverage=0.0, identity=0.0,
        )

    return AlignmentResult(
        score=score, design_seq=design_ext, read_seq=seq, read_qualities=quals,
        strand=strand, blocks=blocks,
        used_doubled_design=any(d1 > len(design.seq) for (_, d1, _, _) in blocks),
        coverage=coverage, identity=identity,
    )


# ---------- Step 4: Mutation extraction ----------

IUPAC_AMBIGUITY = {
    "N": "ACGT", "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
    "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
}


def _compatible(read_char: str, design_char: str) -> bool:
    """True if the (possibly ambiguous) read call is consistent with the
    design base — an ambiguous call that's consistent with the design is
    NOT a mutation, just an uncertain read."""
    if read_char == design_char:
        return True
    possible = IUPAC_AMBIGUITY.get(read_char)
    return possible is not None and design_char in possible


@dataclass
class Mutation:
    position: int          # 0-indexed position in the design coordinate space
    kind: str              # "substitution" | "deletion" | "insertion"
    ref: str
    alt: str
    confidence: Optional[str] = None  # "high" | "low" | None (no quality data)
    min_quality: Optional[int] = None


def extract_mutations(result: AlignmentResult, quality_threshold: int = 20) -> list[Mutation]:
    """Walk the gap-free aligned blocks directly against the original
    sequences (not the printed alignment strings) so multi-base indels
    come out as single events, and ambiguous-but-consistent basecalls
    (e.g. read 'R' vs design 'A') are never reported as mutations."""
    mutations: list[Mutation] = []
    blocks = result.blocks
    quals = result.read_qualities

    for bi, (d0, d1, r0, r1) in enumerate(blocks):
        # substitutions within this gap-free block
        for offset in range(d1 - d0):
            d_pos = d0 + offset
            r_pos = r0 + offset
            d_char = result.design_seq[d_pos]
            r_char = result.read_seq[r_pos]
            if not _compatible(r_char, d_char):
                q = quals[r_pos] if quals else None
                mutations.append(Mutation(
                    position=d_pos, kind="substitution", ref=d_char, alt=r_char,
                    confidence=(None if q is None else ("high" if q >= quality_threshold else "low")),
                    min_quality=q,
                ))

        # gap between this block and the next -> indel
        if bi + 1 < len(blocks):
            nd0, _, nr0, _ = blocks[bi + 1]
            if nd0 > d1:
                # design has bases the read doesn't -> deletion
                ref = result.design_seq[d1:nd0]
                q_flank = None
                if quals:
                    flank = [q for q in (quals[max(0, r1 - 1)], quals[min(len(quals) - 1, r1)]) if q is not None]
                    q_flank = min(flank) if flank else None
                mutations.append(Mutation(
                    position=d1, kind="deletion", ref=ref, alt="-",
                    confidence=(None if q_flank is None else ("high" if q_flank >= quality_threshold else "low")),
                    min_quality=q_flank,
                ))
            elif nr0 > r1:
                # read has bases the design doesn't -> insertion
                alt = result.read_seq[r1:nr0]
                q_ins = min(quals[r1:nr0]) if quals else None
                mutations.append(Mutation(
                    position=d1, kind="insertion", ref="-", alt=alt,
                    confidence=(None if q_ins is None else ("high" if q_ins >= quality_threshold else "low")),
                    min_quality=q_ins,
                ))

    return mutations


# ---------- Step 5: Biological impact ----------

@dataclass
class ImpactReport:
    mutation: Mutation
    feature_name: Optional[str]
    frameshift: bool
    premature_stop: bool
    detail: str


def _find_feature(pos: int, features: list) -> Optional[dict]:
    """Prefer CDS over other annotation types when features overlap —
    CDS is what drives frame/stop-codon analysis, so it should win over
    a promoter or misc_feature covering the same span."""
    matches = [f for f in features if f["start"] <= pos < f["end"]]
    if not matches:
        return None
    cds_matches = [f for f in matches if f["type"] == "CDS"]
    return cds_matches[0] if cds_matches else matches[0]


def _codon_at_position(position: int, feature: dict, design_seq: str):
    """Return (plus_codon, sense_codon, plus_offset, codon_start_design) for
    the codon containing `position`, correctly handling reverse-strand (-1)
    features where the gene reads 5'->3' from feature['end']-1 back down to
    feature['start']. Returns None if the codon would run past the
    annotated feature boundary (e.g. a mutation in the last partial codon),
    since we can't safely translate a partial codon.
    """
    strand = feature.get("strand", 1)
    if strand == 1:
        codon_start = feature["start"] + ((position - feature["start"]) // 3) * 3
        codon_end = codon_start + 3
        if codon_end > feature["end"]:
            return None
        plus_codon = design_seq[codon_start:codon_end]
        plus_offset = position - codon_start
        return plus_codon, plus_codon, plus_offset, codon_start
    else:
        sense_pos = (feature["end"] - 1) - position  # distance from the gene's 5' end
        codon_index = sense_pos // 3
        codon_end_design = feature["end"] - 3 * codon_index
        codon_start_design = codon_end_design - 3
        if codon_start_design < feature["start"]:
            return None
        plus_codon = design_seq[codon_start_design:codon_end_design]
        sense_codon = str(Seq(plus_codon).reverse_complement())
        plus_offset = position - codon_start_design
        return plus_codon, sense_codon, plus_offset, codon_start_design


def assess_impact(mutations: list[Mutation], design: DesignRecord,
                   observed_bases: Optional[dict] = None) -> list[ImpactReport]:
    """observed_bases maps design position -> the base the read actually
    has there (for positions covered by a gap-free aligned block). This is
    what lets two substitutions inside the SAME codon be evaluated
    TOGETHER — e.g. AAC -> TAC and AAC -> AAA each individually translate
    fine, but combined (AAC -> TAA) create a real premature stop that
    neither substitution shows on its own. Looking up each codon position
    independently from the read (rather than applying one Mutation's alt
    base to an otherwise-unchanged design codon) catches this correctly."""
    observed_bases = observed_bases or {}
    reports = []
    # net indel length so far, per feature, to detect frame-breaking totals
    indel_running_total = {}

    for m in mutations:
        feature = _find_feature(m.position, design.features)
        fname = feature["name"] if feature else None

        frameshift = False
        premature_stop = False
        detail = ""

        if feature and feature["type"] == "CDS" and isinstance(feature.get("start"), int):
            if feature.get("compound_location"):
                detail = (
                    f"{m.kind.capitalize()} at position {m.position} inside {fname}, which has a "
                    f"multi-segment location (e.g. spans the plasmid origin, or has an intron/exon "
                    f"structure) — frame and stop-codon analysis is skipped here; verify manually."
                )
            else:
                key = fname
                delta = 0
                if m.kind == "deletion":
                    delta = -len(m.ref)
                elif m.kind == "insertion":
                    delta = len(m.alt)
                indel_running_total[key] = indel_running_total.get(key, 0) + delta

                if delta != 0 and indel_running_total[key] % 3 != 0:
                    frameshift = True
                    detail = (
                        f"{m.kind.capitalize()} at position {m.position} inside {fname} "
                        f"shifts the reading frame — downstream amino acids will be wrong."
                    )
                elif m.kind == "substitution":
                    codon_info = _codon_at_position(m.position, feature, design.seq)
                    if codon_info:
                        plus_codon, sense_codon, plus_offset, codon_start = codon_info
                        try:
                            # Reconstruct the codon from what the read ACTUALLY
                            # shows at all three positions, not just this one
                            # mutation's alt base — catches compound effects.
                            observed_plus = "".join(
                                observed_bases.get(codon_start + k, plus_codon[k])
                                for k in range(3)
                            )
                            mutated_sense = (
                                observed_plus if feature.get("strand", 1) == 1
                                else str(Seq(observed_plus).reverse_complement())
                            )
                            original_aa = str(Seq(sense_codon).translate())
                            mutated_aa = str(Seq(mutated_sense).translate())
                            if mutated_aa == "*" and original_aa != "*":
                                premature_stop = True
                                detail = (
                                    f"Substitution at position {m.position} in {fname} "
                                    f"turns codon {sense_codon}->{mutated_sense}, creating a premature stop codon."
                                )
                        except Exception:
                            pass

        elif feature and feature["type"] == "gene":
            # A "gene" annotation with no co-located CDS feature at this
            # position. We deliberately do NOT treat the gene span as if
            # it were the CDS — genes can include UTRs/introns that would
            # throw off frame math — but silently falling through to a
            # bare "substitution inside X" message would look identical to
            # a real "checked, nothing wrong" result. Say plainly that
            # frame/stop-codon analysis wasn't attempted here.
            detail = (
                f"{m.kind.capitalize()} at position {m.position} inside gene {fname}, but no "
                f"separate CDS feature was found here — frame and stop-codon analysis needs a CDS "
                f"annotation and was NOT run. If this gene is coding, verify manually."
            )

        if not detail:
            location = f"inside {fname}" if fname else "outside any annotated feature"
            detail = f"{m.kind.capitalize()} at position {m.position} {location}."

        if m.confidence == "low":
            detail += (
                f" ⚠ Low base-call confidence (Phred {m.min_quality}) at this position — "
                f"could be a sequencing artifact rather than a real mutation. Verify against the trace."
            )

        reports.append(ImpactReport(
            mutation=m, feature_name=fname, frameshift=frameshift,
            premature_stop=premature_stop, detail=detail,
        ))
    return reports


# ---------- Orchestration ----------

def run_pipeline(design_bytes: bytes, design_filename: str,
                  reality_bytes: bytes, reality_filename: str) -> dict:
    design = parse_design(design_bytes, design_filename)
    read = parse_reality(reality_bytes, reality_filename)
    read = quality_trim(read)

    alignment = align_best_orientation(design, read)

    base_result = {
        "alignment_score": alignment.score,
        "strand": alignment.strand,
        "used_circular_wrap": alignment.used_doubled_design,
        "coverage": alignment.coverage,
        "identity": alignment.identity,
    }

    if alignment.coverage < MIN_COVERAGE or alignment.identity < MIN_IDENTITY:
        # Don't report mutations off an alignment we don't trust — this is
        # what an all-N read, an unrelated random sequence, or a wrong-file
        # upload all look like. A wrong "no mutations" or a made-up mutation
        # list is worse than telling the user we can't confidently compare.
        return {
            **base_result,
            "status": "no_confident_match",
            "mutation_count": None,
            "impacts": [],
        }

    mutations = extract_mutations(alignment)

    # Mutation positions may fall in the wrapped-copy region of design_seq
    # (design.seq + design.seq[:wrap_pad]) when the read spans the origin.
    # Fold them back into real plasmid coordinates before feature lookup —
    # design_ext[pos] == design.seq[pos % design_len] by construction, so
    # this only changes the reported coordinate, not the ref/alt content.
    design_len = len(design.seq)
    for m in mutations:
        m.position %= design_len

    # Build a design-position -> actual-observed-base map from the winning
    # alignment's gap-free blocks, so assess_impact can reconstruct a full
    # codon from what the read really shows (catching compound mutations
    # within one codon) instead of changing one base at a time.
    observed_bases = {}
    for (d0, d1, r0, r1) in alignment.blocks:
        for offset in range(d1 - d0):
            observed_bases[(d0 + offset) % design_len] = alignment.read_seq[r0 + offset]

    impacts = assess_impact(mutations, design, observed_bases)

    return {
        **base_result,
        "status": "ok",
        "mutation_count": len(mutations),
        "impacts": [
            {
                "position": r.mutation.position,
                "kind": r.mutation.kind,
                "ref": r.mutation.ref,
                "alt": r.mutation.alt,
                "feature": r.feature_name,
                "frameshift": r.frameshift,
                "premature_stop": r.premature_stop,
                "confidence": r.mutation.confidence,
                "min_quality": r.mutation.min_quality,
                "detail": r.detail,
            }
            for r in impacts
        ],
    }
