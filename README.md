# CrisprGuard

**Zero-dependency Cas9/Cas12a CRISPR Guide RNA Designer, CFD Off-Target Kinetics Engine & Bae MMEJ Microhomology Frameshift Predictor**

CrisprGuard is a pure Python 3.10+ bioinformatics toolkit, CLI, and dark-lab interactive web workbench for designing high-efficiency, high-specificity CRISPR-Cas9 and CRISPR-Cas12a guide RNAs while predicting post-cleavage DNA repair frameshift outcomes.

---

## Key Scientific Capabilities

### 1. Dual-Endonuclease Bidirectional PAM Scanning (`design.py`)
- **SpCas9 (*Streptococcus pyogenes*)**: Scans both `+` and `-` strands for canonical `5'-[20 nt spacer]-NGG-3'` loci and computes exact blunt double-strand break (DSB) coordinates at **-3 bp upstream of the PAM**.
- **LbCas12a / Cpf1 (*Lachnospiraceae bacterium*)**: Scans both strands for `5'-TTTV-[23 nt spacer]-3'` (`V ∈ {A, C, G}`) loci and tracks the staggered **5-nt 5' overhang** cleavage sites (**+18 nt** non-target strand / **+23 nt** target strand downstream of the PAM).

### 2. Doench Rule-Set / Kim DeepCpf1-Inspired On-Target Scoring (`0–100`)
Each candidate spacer is scored via an interpretable biophysical logit model mapped smoothly to `0–100`:
- **GC Balance**: Rewards optimal 40–65% GC (SpCas9) or 35–60% GC (LbCas12a) while penalizing extreme GC content.
- **Position-Specific Nucleotide Preferences**:
  - *SpCas9*: Favors `G` at position 20 (immediately upstream of `NGG`), `C` at position 18, `5' G` for U6 transcription initiation, and `CGG`/`AGG` PAM context while penalizing `C` at position 20 and PAM-proximal poly-U runs (positions 17–20).
  - *LbCas12a*: Models `TTTC`/`TTTA` > `TTTG` preference and 5' seed position biases (strongly disfavoring `T` at position 1).
- **Pol III Terminator Detection**: Heavily penalizes and flags internal `TTTT` motifs that prematurely terminate U6/H1 RNA Polymerase III transcription.
- **Secondary Structure & Homopolymers**: Computes internal self-folding hairpin stem lengths & estimated $\Delta G$, spacer-scaffold complementarity, melting temperature ($T_m$), and `GGGG` G-quadruplex / homopolymer runs.

### 3. PAM-Proximal Seed-Weighted CFD Off-Target Kinetics (`offtarget.py`)
Implements a Cutting Frequency Determination (CFD) model grounded in Doench et al. (2016) and Kleinstiver/Kim Cas12a kinetics:
- Scans both strands of any multi-FASTA reference background for PAM-adjacent loci with up to `max_mismatches` (default `<= 4`) mismatches.
- **Opposite Seed Polarity**:
  - *SpCas9*: Steep mismatch penalties across the **3' PAM-proximal seed** (positions 11–20; core 13–20).
  - *LbCas12a*: Steep mismatch penalties across the **5' PAM-proximal seed** (positions 1–8).
- **Thermodynamic RNA:DNA Mismatch Pair Kinetics**: Models exact guide-RNA to target-strand DNA base-pair mismatches (`rG:dT` and `rU:dG` wobble mismatches retain substantially higher cleavage tolerance than bulky purine-purine or `rC:dC` transversions), plus consecutive seed-mismatch cooperativity.
- **Non-Canonical PAM Weights**: Scores canonical `NGG` (`1.00`), `NAG` (`0.259`), `NGA` (`0.069`), `NGC`/`NGT`/`NCG` for SpCas9, and `TTTV` (`1.00`), `TTTT` (`0.20`), `CTTV`/`TCTV`/`TTCV` for LbCas12a.
- **Aggregate MIT/CFD Specificity Score (`0–100`)**:

$$\text{Specificity} = \frac{100}{1 + \sum_{i} \text{CFD}_i}$$

### 4. Bae et al. Microhomology-Mediated End Joining (MMEJ) Frameshift Predictor (`design.py`)
Predicts gene knockout efficiency directly from sequence flanking the exact DSB cut site (±25 bp):
- Enumerates maximal 2–8 bp flanking microhomology repeats bridging the DSB.
- Scores each microhomology deletion pattern using hydrogen-bond GC weighting, cooperative annealing length scaling, and exponential deletion-length decay:

$$S_{\text{pattern}} = 100 \cdot (k + n_{\text{GC}}) \cdot 1.85^{k-2} \cdot e^{-\Delta L / 15}$$

- Distinguishes **frameshift (`ΔL mod 3 ≠ 0`)** vs. **in-frame (`ΔL mod 3 = 0`)** deletions to report **Out-of-Frame Knockout Probability (%)** alongside character-aligned post-deletion sequence diagrams.

### 5. Interactive Dark-Lab Web Workbench (`web/`)
- **Bilingual EN / TR interface** with instant switching.
- **Curated Biological Presets**: Human **TP53 Exon 5**, **BRCA1 Exon 11**, **PCSK9 Exon 1**, and **eGFP Reporter Knockout**, each paired with homologous paralog/pseudogene off-target reference FASTAs.
- **Interactive HTML5 Canvas 2D Genomic Target Track**: Zoomable double-stranded DNA map with strand-packed sgRNA/crRNA arrows color-coded by composite score, PAM caps, and exact DSB cut-line overlays.
- **Guide Inspector Drawer**: Circular efficiency/specificity/MMEJ gauges, biophysical logit breakdown, clickable Microhomology Deletion Alignment diagram, Off-Target Mismatch Matrix with seed/wobble highlighting, and one-click CSV/FASTA export.

---

## Installation

CrisprGuard requires **Python 3.10+** and has **zero third-party runtime dependencies** (standard library only).

```bash
pip install -e .
```

Or run directly from source without installation:

```bash
PYTHONPATH=src python -m crisprguard --help
```

---

## CLI Quick Start

### Design & rank guides on a curated human exon preset
```bash
# SpCas9 guides for Human TP53 Exon 5 (with paired TP63/TP73 off-target background)
crisprguard design --sample tp53_exon5 --enzyme spcas9 --top 10

# LbCas12a (TTTV) guides for Human BRCA1 Exon 11 exported as CSV or FASTA
crisprguard design --sample brca1_exon11 --enzyme cas12a --format csv -o brca1_cas12a.csv
crisprguard design --sample egfp_reporter --enzyme spcas9 --format fasta -o egfp_guides.fasta
```

### Run standalone CFD off-target kinetics
```bash
crisprguard offtarget \
  --guide GACCTGCCCTGTGCAGCTGT \
  --sample tp53_exon5 \
  --enzyme spcas9 \
  --max-mismatches 4
```

### Predict Bae MMEJ frameshift probability around a DSB cut site
```bash
crisprguard mmej --sample pcsk9_exon1 --cut-pos 142 --window 25
```

### Launch the interactive Dark-Lab Web Workbench
```bash
crisprguard serve --host 127.0.0.1 --port 8765
```

---

## Python API Example

```python
from crisprguard import design_guides, get_sample, predict_mmej, score_offtarget_cfd

preset = get_sample("tp53_exon5")
report = design_guides(
    sequence=preset.target_fasta,
    enzyme="spcas9",
    background=preset.background_fasta,
    max_mismatches=4,
    top_n=5,
)

top = report.guides[0]
print(f"{top.id} ({top.strand}) {top.spacer} [{top.pam}]")
print(
    f"Composite={top.composite_score} | On-Target={top.on_target_score} | "
    f"Specificity={top.specificity_score} | MMEJ OOF={top.mmej.out_of_frame_probability}%"
)

# Direct CFD comparison: PAM-proximal seed mismatch vs PAM-distal mismatch
cfd_seed, *_ = score_offtarget_cfd("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGA", "AGG")
cfd_distal, *_ = score_offtarget_cfd("ACGTACGTACGTACGTACGT", "TCGTACGTACGTACGTACGT", "AGG")
print(f"Seed mismatch CFD={cfd_seed:.4f} < Distal mismatch CFD={cfd_distal:.4f}")
```

---

## Project Structure

```text
crisprguard/
├── pyproject.toml
├── LICENSE
├── README.md
├── .github/workflows/ci.yml
├── src/crisprguard/
│   ├── __init__.py       # Public Python API exports
│   ├── __main__.py       # `python -m crisprguard` runner
│   ├── offtarget.py      # CFD mismatch kinetics & bidirectional off-target scanner
│   ├── design.py         # SpCas9/LbCas12a PAM scanner, on-target scorer & Bae MMEJ predictor
│   ├── samples.py        # Curated TP53, BRCA1, PCSK9, eGFP presets + off-target decoys
│   ├── cli.py            # CLI subcommands (design, offtarget, mmej, serve)
│   └── web/
│       ├── server.py     # Pure-stdlib HTTP server & testable JSON API router
│       └── static/       # Dark-lab interactive UI (index.html, styles.css, app.js)
└── tests/
    └── test_crisprguard.py
```

---

## Running Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

---

## Scientific References

1. **Doench, J. G. et al.** (2016). *Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9.* Nature Biotechnology, 34(2), 184–191. (Rule Set 2 & CFD off-target scoring)
2. **Doench, J. G. et al.** (2014). *Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation.* Nature Biotechnology, 32(12), 1262–1267.
3. **Bae, S., Kweon, J., Kim, H. S., & Kim, J.-S.** (2014). *Microhomology-based choice of Cas9 nuclease target sites.* Nature Methods, 11(7), 705–706. (MMEJ microhomology out-of-frame predictor)
4. **Kim, H. K. et al.** (2018). *Deep learning improves prediction of CRISPR-Cpf1 guide RNA activity.* Nature Biotechnology, 36(3), 239–241. (LbCas12a/Cpf1 `TTTV` activity rules)
5. **Kleinstiver, B. P. et al.** (2016). *Genome-wide specificities of CRISPR-Cas Cpf1 nucleases in human cells.* Nature Biotechnology, 34(8), 869–874.
6. **Hsu, P. D. et al.** (2013). *DNA targeting specificity of RNA-guided Cas9 nucleases.* Nature Biotechnology, 31(9), 827–832.

---

## License

MIT License — Copyright (c) 2026 Muhammed Emir Esti.
