# PlasmidGuard

**Live demo: https://plasmidguard.streamlit.app/** (use "Try an example", no files needed)

Checks whether a sequenced plasmid matches the one you designed. It finds every difference and,
for mutations inside annotated genes, flags **frameshifts** and **premature stop codons**,
including for genes on the reverse strand.

> **Status: early prototype.** Tested on real plasmid files (Addgene, NCBI) and a real Sanger
> trace, with mutations introduced in software. **Not yet validated against independently
> confirmed lab results.** If you have a sequenced plasmid where the answer is already known and
> are willing to test it, please open an issue or get in touch.

## Run it on your own computer

1. **Unzip the folder first** (do not run it from inside the zip).
2. Install Python 3.10 or newer if you don't have it: https://www.python.org/downloads/
   (on Windows, tick "Add python.exe to PATH").
3. Start it:
   - **Windows:** double-click `run_windows.bat`
   - **Mac / Linux:** open a terminal in this folder and run `bash run_mac_linux.sh`
4. The first run downloads what it needs (1-3 minutes, needs internet). Later runs start in seconds.
5. A browser tab opens. Keep the black window open while you use it; close it to stop.

If the tab doesn't open, go to http://localhost:8501

## Try it with the examples (folder `examples/`)

Upload the **design** file in the left box and the **sequencing result** in the right box, then click Analyze.

| Design | Sequencing result | You should see |
|---|---|---|
| `test_design.gb` | `test_reality_perfect.fasta` | Green: no mutations detected |
| `test_design.gb` | `test_reality_mutated.fasta` | 1 mutation: deletion at position 220, frameshift, in `testGFP` |
| `test_design.gb` | `test_reality_garbage.fasta` | "Doesn't confidently match the design" (not a mutation list) |
| `real_plasmid1_design.gb` | `real_plasmid1_mutated.fasta` | Deletion at position 5400 in `bla`, frameshift (reverse-strand gene) |
| `real_plasmid2_design.gb` | `real_plasmid2_mutated.fasta` | Substitution at position 6913 in `pla` |

The two `real_plasmid` designs are genuine files (an Addgene lentiviral vector and an NCBI
Yersinia pestis plasmid). The mutations in the matching result files were introduced in software.

## Your own files

- Design: GenBank (`.gb`, `.gbk`) or FASTA. Use GenBank to get gene-level impact calls.
- Sequencing result: `.ab1` (Sanger trace), `.fasta`, or `.fastq`. One read per file.
- Positions are 0-based (the first base is position 0). Sequence viewers such as SnapGene show 1-based positions, so add 1 when comparing.

## Privacy

- Runs on your own computer and is set to accept connections from this computer only.
- Usage statistics to Streamlit are switched off.
- Your files are held in memory for the session. The code never writes them to disk, a database
  or any outside service, and does no logging or caching of their contents.
- If you instead put this on a website, uploads pass through that host's servers. Keep sensitive
  or unpublished sequences on a local copy.

## Known limitations (please read)

- **Not yet validated against independently confirmed lab results.** It has been tested on real
  plasmid files and a real Sanger trace, with mutations introduced in software.
  Treat every result as something to confirm, not as a final answer.
- One design vs. one read at a time. It does not combine several reads covering one construct.
- Where two coding features overlap on opposite strands, it does not choose between them.
- Coding features that span the plasmid origin, or that are annotated only as `gene` (no `CDS`),
  get a "verify manually" note instead of frameshift/stop-codon analysis.
- Size limits: design up to 100,000 bp, read up to 6,000 bp.
- "No confident match" means the two files don't line up. It is not a mutation report.

## Troubleshooting

- **Windows says "Python was not found":** install Python (step 2) and double-click again.
- **Setup failed:** check your internet connection, delete the `.venv` folder here, and try again.
- **Port already in use:** close other PlasmidGuard windows and start again.
- **Mac blocks the script:** use the terminal command in step 3 instead of double-clicking.
