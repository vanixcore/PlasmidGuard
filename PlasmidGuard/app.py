"""
PlasmidGuard — Streamlit frontend.

Run locally with: python -m streamlit run app.py
"""

from pathlib import Path

import streamlit as st
import pandas as pd

from engine import run_pipeline, PlasmidGuardError

st.set_page_config(page_title="PlasmidGuard", page_icon="🧬", layout="wide")

EXAMPLE_DIR = Path(__file__).parent / "examples"

# label -> (design file, sequencing file, what the visitor should expect)
EXAMPLES = {
    "Perfect match": (
        "test_design.gb", "test_reality_perfect.fasta",
        "The read is an exact copy of part of the design. Expected: no mutations detected.",
    ),
    "Frameshift deletion in a gene": (
        "test_design.gb", "test_reality_mutated.fasta",
        "One base is missing inside a gene. Expected: a deletion at position 220 flagged as a frameshift.",
    ),
    "Wrong file (unrelated sequence)": (
        "test_design.gb", "test_reality_garbage.fasta",
        "A random sequence unrelated to the design. Expected: rejected as 'doesn't match', "
        "not a made-up mutation list.",
    ),
    "Real plasmid: frameshift in bla (reverse-strand gene)": (
        "real_plasmid1_design.gb", "real_plasmid1_mutated.fasta",
        "A real Addgene lentiviral vector with one base removed inside its ampicillin-resistance gene. "
        "Expected: deletion at position 5400 in bla, flagged as a frameshift.",
    ),
    "Real plasmid: substitution in pla": (
        "real_plasmid2_design.gb", "real_plasmid2_mutated.fasta",
        "A real NCBI Yersinia pestis plasmid with one base changed inside the pla gene. "
        "Expected: substitution at position 6913 in pla, not flagged as a frameshift or stop.",
    ),
}


def available_examples() -> dict:
    return {
        label: spec for label, spec in EXAMPLES.items()
        if (EXAMPLE_DIR / spec[0]).exists() and (EXAMPLE_DIR / spec[1]).exists()
    }


def analyze_and_render(design_bytes, design_name, reality_bytes, reality_name) -> None:
    with st.spinner("Aligning and scanning for mutations..."):
        try:
            result = run_pipeline(design_bytes, design_name, reality_bytes, reality_name)
        except PlasmidGuardError as e:
            st.error(str(e))
            return
        except Exception as e:
            st.error(f"Something went wrong processing these files: {e}")
            return

    if result["status"] == "no_confident_match":
        st.error(
            f"This sequencing read doesn't confidently match the design "
            f"({result['identity']*100:.0f}% identity over {result['coverage']*100:.0f}% of the read). "
            f"Double-check you uploaded the right files — this isn't a mutation report, "
            f"it's a 'these two sequences don't line up' signal."
        )
        return

    strand_label = "forward" if result["strand"] == 1 else "reverse complement"
    wrap_note = " (read spans the plasmid origin)" if result["used_circular_wrap"] else ""
    st.caption(f"Best alignment: {strand_label} strand{wrap_note} · alignment score {result['alignment_score']:.0f}")

    n = result["mutation_count"]
    if n == 0:
        st.success("No mutations detected. The sequencing result matches the design.")
        return

    high_conf_danger = [
        i for i in result["impacts"]
        if (i["frameshift"] or i["premature_stop"]) and i["confidence"] != "low"
    ]
    low_conf = [i for i in result["impacts"] if i["confidence"] == "low"]

    if high_conf_danger:
        st.error(f"🚨 {len(high_conf_danger)} mutation(s) with likely functional impact.")
    else:
        st.warning(f"{n} mutation(s) found, none flagged as high-confidence functional breaks.")

    if low_conf:
        st.caption(
            f"⚠ {len(low_conf)} of these are at low-confidence base calls — "
            f"could be sequencing noise rather than real mutations. Check the trace before acting on them."
        )

    df = pd.DataFrame(result["impacts"])
    df = df.rename(columns={
        "position": "Position", "kind": "Type", "ref": "Ref", "alt": "Alt",
        "feature": "Feature", "frameshift": "Frameshift", "premature_stop": "Premature Stop",
        "confidence": "Confidence", "min_quality": "Min Phred Q", "detail": "Detail",
    })
    df = df.sort_values("Position")

    st.dataframe(
        df,
        width="stretch",
        hide_index=True,
        column_config={
            "Frameshift": st.column_config.CheckboxColumn(),
            "Premature Stop": st.column_config.CheckboxColumn(),
        },
    )
    st.caption(
        "Positions are 0-based (the first base is position 0). Sequence viewers such as SnapGene "
        "show 1-based positions, so add 1 when comparing."
    )

    with st.expander("Full mutation detail"):
        for impact in sorted(result["impacts"], key=lambda x: x["position"]):
            icon = "🚨" if (impact["frameshift"] or impact["premature_stop"]) else "•"
            st.markdown(f"{icon} **Position {impact['position']}** — {impact['detail']}")


st.title("🧬 PlasmidGuard")
st.caption("Compare your plasmid design with your sequencing result and see exactly what differs.")

st.info(
    "Early prototype. It has been tested on real plasmid files and a real Sanger trace with mutations "
    "introduced in software, but **not yet validated against independently confirmed lab results**. "
    "Treat every result as something to double-check."
)

with st.expander("🔒 Privacy — where do my files go?"):
    st.markdown(
        "- The code holds your files **in server memory only, for your session**. It never writes them "
        "to disk, a database, or any outside service, and does no logging or caching of their contents.\n"
        "- **If this is a hosted copy**, your upload is still sent over the internet to the machine "
        "running it, which belongs to the hosting provider. This tool cannot control what the provider's "
        "infrastructure does.\n"
        "- **For unpublished or sensitive sequences, run PlasmidGuard on your own computer** "
        "(see the README). Then nothing leaves your machine. On a public copy, use the built-in examples."
    )

examples = available_examples()
mode = "Upload my files"
if examples:
    mode = st.radio("What would you like to do?", ["Upload my files", "Try an example"], horizontal=True)

if mode == "Upload my files":
    st.caption("Design = the plasmid as intended. Sequencing result = what came back from the sequencer.")
    col1, col2 = st.columns(2)
    with col1:
        design_file = st.file_uploader(
            "Design file (.gb, .gbk, .fasta)",
            type=["gb", "gbk", "fasta", "fa"],
            help="Your intended plasmid sequence — ideally a GenBank file so we can map mutations to features.",
        )
    with col2:
        reality_file = st.file_uploader(
            "Sequencing result (.ab1, .fasta, .fastq)",
            type=["ab1", "fasta", "fa", "fastq", "fq"],
            help="What actually came back from the sequencer.",
        )

    run_clicked = st.button("Analyze", type="primary", disabled=not (design_file and reality_file))

    if not (design_file and reality_file):
        st.info("Upload both files to run the comparison.")

    if run_clicked and design_file and reality_file:
        analyze_and_render(
            design_file.getvalue(), design_file.name,
            reality_file.getvalue(), reality_file.name,
        )
else:
    choice = st.selectbox("Example", list(examples))
    design_name, reality_name, blurb = examples[choice]
    st.caption(blurb)
    st.caption("All example mutations were introduced by software for demonstration — they are not results from a real experiment.")
    if st.button("Run example", type="primary"):
        analyze_and_render(
            (EXAMPLE_DIR / design_name).read_bytes(), design_name,
            (EXAMPLE_DIR / reality_name).read_bytes(), reality_name,
        )

st.divider()
st.caption(
    "PlasmidGuard does not save your files to disk or send them to any third party. "
    "If you are using a hosted copy, your upload still passes through the host's servers — "
    "for unpublished sequences, run it locally (see README)."
)
