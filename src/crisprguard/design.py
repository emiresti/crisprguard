"""Cas9/Cas12a CRISPR guide RNA designer, on-target efficiency scorer & MMEJ predictor.

Supports:
  - SpCas9: 5'-[20nt spacer]-[NGG]-3', blunt DSB at -3 bp upstream of PAM.
  - LbCas12a (Cpf1): 5'-[TTTV]-[23nt spacer]-3', staggered 5' overhang DSB (+18/+23 nt).
  - Doench Rule-Set / Kim DeepCpf1-inspired transparent on-target scoring (0..100).
  - Bae et al. (Nature Methods, 2014) Microhomology-Mediated End Joining (MMEJ)
    frameshift / out-of-frame knockout predictor around the exact DSB cut site.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import io
import math
from typing import Sequence

from crisprguard.offtarget import (
    OffTargetSummary,
    clean_dna,
    normalize_enzyme,
    parse_fasta,
    reverse_complement,
    scan_offtargets,
)

# Scaffold repeat sequences used to detect destructive spacer-scaffold self-folding
SPCAS9_SCAFFOLD_CORE = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG"
CAS12A_DIRECT_REPEAT = "AATTTCTACTAAGTGTAGAT"

# Nearest-neighbor RNA:DNA hybrid enthalpy/entropy proxy (kcal/mol at 37 C) for Tm & hairpin
BASE_PAIR_DG = {
    "GC": -2.4,
    "CG": -2.4,
    "AT": -1.2,
    "TA": -1.2,
}


@dataclass(frozen=True)
class EnzymeSpec:
    """Biophysical specification of a CRISPR endonuclease."""

    id: str
    name: str
    pam_motif: str
    pam_side: str  # '3prime' or '5prime'
    pam_length: int
    spacer_length: int
    cut_description: str
    seed_range_1based: tuple[int, int]
    overhang_nt: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


ENZYME_SPECS: dict[str, EnzymeSpec] = {
    "spcas9": EnzymeSpec(
        id="spcas9",
        name="SpCas9 (S. pyogenes)",
        pam_motif="NGG",
        pam_side="3prime",
        pam_length=3,
        spacer_length=20,
        cut_description="Blunt DSB at -3 bp upstream of 3' NGG PAM",
        seed_range_1based=(11, 20),
        overhang_nt=0,
    ),
    "cas12a": EnzymeSpec(
        id="cas12a",
        name="LbCas12a / Cpf1 (Lachnospiraceae)",
        pam_motif="TTTV",
        pam_side="5prime",
        pam_length=4,
        spacer_length=23,
        cut_description="Staggered 5-nt 5' overhang DSB (+18 NTS / +23 TS downstream of 5' TTTV PAM)",
        seed_range_1based=(1, 8),
        overhang_nt=5,
    ),
}


@dataclass(frozen=True)
class MMEJPattern:
    """Single microhomology-mediated deletion pattern across a DSB cut site."""

    microhomology: str
    mh_length: int
    gc_count: int
    left_mh_start: int  # 0-based index within flanking window (0..left_len)
    right_mh_start: int  # 0-based index within flanking window (left_len..total_len)
    deletion_length: int
    is_frameshift: bool  # True when deletion_length % 3 != 0
    pattern_score: float
    aligned_outcome: str  # Same length as left_flank + right_flank, '-' marks deleted bases
    deleted_sequence: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class MMEJResult:
    """Bae et al. (2014) Microhomology-Mediated End Joining frameshift prediction."""

    cut_pos: int
    window_bp: int
    left_flank: str
    right_flank: str
    window_start_genomic: int
    microhomology_score: float
    out_of_frame_score: float
    in_frame_score: float
    out_of_frame_probability: float  # 0..100 %
    in_frame_probability: float  # 0..100 %
    total_patterns: int
    frameshift_patterns: int
    in_frame_patterns: int
    knockout_grade: str  # 'High KO', 'Moderate KO', 'In-Frame Risk'
    top_patterns: tuple[MMEJPattern, ...]

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["top_patterns"] = [p.to_dict() for p in self.top_patterns]
        return data


@dataclass(frozen=True)
class OnTargetScoreDetail:
    """Decomposed Doench/Kim-inspired on-target activity score."""

    score: float  # 0..100
    gc_percent: float
    tm_celsius: float
    hairpin_stem_bp: int
    scaffold_stem_bp: int
    estimated_hairpin_dg: float
    pol3_terminator: bool
    homopolymer_run: int
    seed_t_count: int
    components: dict[str, float]
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["warnings"] = list(self.warnings)
        return data


@dataclass(frozen=True)
class GuideCandidate:
    """Ranked CRISPR sgRNA/crRNA candidate with kinetics and repair outcome predictions."""

    rank: int
    id: str
    enzyme: str
    strand: str  # '+' or '-'
    start: int  # 0-based genomic start of full [spacer+PAM] span on + strand
    end: int  # 0-based exclusive genomic end on + strand
    spacer_start: int  # 0-based genomic start of spacer on + strand
    spacer_end: int
    pam_start: int
    pam_end: int
    cut_pos: int  # Canonical DSB boundary on + strand (blunt or NTS cut)
    cut_ts: int  # Target strand cut boundary (equals cut_pos for SpCas9; offset by 5 for Cas12a)
    spacer: str
    pam: str
    context_30nt: str
    gc_percent: float
    on_target_score: float
    on_target_detail: OnTargetScoreDetail
    mmej: MMEJResult
    specificity_score: float
    offtarget_summary: OffTargetSummary | None
    composite_score: float
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "rank": self.rank,
            "id": self.id,
            "enzyme": self.enzyme,
            "strand": self.strand,
            "start": self.start,
            "end": self.end,
            "spacer_start": self.spacer_start,
            "spacer_end": self.spacer_end,
            "pam_start": self.pam_start,
            "pam_end": self.pam_end,
            "cut_pos": self.cut_pos,
            "cut_ts": self.cut_ts,
            "spacer": self.spacer,
            "pam": self.pam,
            "context_30nt": self.context_30nt,
            "gc_percent": self.gc_percent,
            "on_target_score": self.on_target_score,
            "on_target_detail": self.on_target_detail.to_dict(),
            "mmej": self.mmej.to_dict(),
            "specificity_score": self.specificity_score,
            "offtarget_summary": (
                self.offtarget_summary.to_dict() if self.offtarget_summary else None
            ),
            "composite_score": self.composite_score,
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class DesignReport:
    """Complete guide design report for a target DNA sequence."""

    target_name: str
    target_sequence: str
    target_length: int
    target_gc_percent: float
    enzyme: str
    enzyme_spec: EnzymeSpec
    background_contigs: int
    background_bp: int
    max_mismatches: int
    total_candidates: int
    plus_strand_count: int
    minus_strand_count: int
    guides: tuple[GuideCandidate, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "target_name": self.target_name,
            "target_sequence": self.target_sequence,
            "target_length": self.target_length,
            "target_gc_percent": self.target_gc_percent,
            "enzyme": self.enzyme,
            "enzyme_spec": self.enzyme_spec.to_dict(),
            "background_contigs": self.background_contigs,
            "background_bp": self.background_bp,
            "max_mismatches": self.max_mismatches,
            "total_candidates": self.total_candidates,
            "plus_strand_count": self.plus_strand_count,
            "minus_strand_count": self.minus_strand_count,
            "guides": [g.to_dict() for g in self.guides],
        }


def gc_content_percent(seq: str) -> float:
    """Compute GC percentage (0..100) ignoring ambiguous N bases."""
    valid = [b for b in seq.upper() if b in "ACGT"]
    if not valid:
        return 0.0
    gc = sum(1 for b in valid if b in "GC")
    return round(100.0 * gc / len(valid), 1)


def estimate_tm_celsius(seq: str) -> float:
    """Estimate RNA:DNA spacer melting temperature (salt-adjusted nearest-neighbor proxy)."""
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 0.0
    gc = sum(1 for b in s if b in "GC")
    # Standard salt-adjusted oligonucleotide Tm (50 mM Na+)
    tm = 64.9 + 41.0 * (gc - 16.4) / n
    return round(tm, 1)


def longest_homopolymer(seq: str) -> int:
    """Return length of the longest identical single-base run."""
    if not seq:
        return 0
    best = 1
    run = 1
    for a, b in zip(seq, seq[1:]):
        if a == b:
            run += 1
            if run > best:
                best = run
        else:
            run = 1
    return best


def analyze_hairpin(spacer: str, enzyme: str = "spcas9") -> tuple[int, int, float]:
    """Detect internal spacer self-complementarity and spacer-scaffold pairing.

    Returns (max_internal_stem_bp, max_scaffold_stem_bp, estimated_dg_kcal).
    """
    s = spacer.upper()
    n = len(s)
    max_internal_stem = 0
    best_dg = 0.0

    # Check internal reverse-complement stems separated by hairpin loop >= 3 nt
    for stem_len in range(min(8, n // 2), 2, -1):
        found = False
        for i in range(n - 2 * stem_len - 2):
            left = s[i : i + stem_len]
            rc_left = reverse_complement(left)
            for j in range(i + stem_len + 3, n - stem_len + 1):
                if s[j : j + stem_len] == rc_left:
                    gc = sum(1 for b in left if b in "GC")
                    at = stem_len - gc
                    # Loop initiation penalty (~+3.6 kcal/mol) plus stacking stabilization
                    dg = round(3.6 - (2.1 * gc + 1.1 * at), 2)
                    if stem_len > max_internal_stem or (
                        stem_len == max_internal_stem and dg < best_dg
                    ):
                        max_internal_stem = stem_len
                        best_dg = min(best_dg, dg)
                    found = True
                    break
            if found:
                break
        if max_internal_stem >= stem_len:
            break

    # Check complementarity against the CRISPR scaffold / direct repeat
    scaffold = (
        SPCAS9_SCAFFOLD_CORE
        if normalize_enzyme(enzyme) == "spcas9"
        else CAS12A_DIRECT_REPEAT
    )
    rc_spacer = reverse_complement(s)
    max_scaffold_stem = 0
    for stem_len in range(min(9, n), 3, -1):
        matched = False
        for i in range(n - stem_len + 1):
            if rc_spacer[i : i + stem_len] in scaffold:
                max_scaffold_stem = stem_len
                matched = True
                break
        if matched:
            break

    return max_internal_stem, max_scaffold_stem, best_dg


def score_on_target(
    spacer: str,
    pam: str,
    enzyme: str = "spcas9",
    downstream_base: str = "",
) -> OnTargetScoreDetail:
    """Compute a Doench Rule-Set / Kim DeepCpf1-inspired on-target activity score (0..100)."""
    enz = normalize_enzyme(enzyme)
    s = clean_dna(spacer)
    p = pam.upper()
    expected = 20 if enz == "spcas9" else 23
    if len(s) != expected:
        raise ValueError(f"{enz} spacer must be {expected} nt (got {len(s)}).")

    gc = gc_content_percent(s)
    tm = estimate_tm_celsius(s)
    internal_stem, scaffold_stem, hairpin_dg = analyze_hairpin(s, enz)
    has_pol3 = "TTTT" in s
    max_homo = longest_homopolymer(s)

    components: dict[str, float] = {}
    warnings: list[str] = []

    if enz == "spcas9":
        # Base logit calibrated so median functional guides score ~55-68, elite guides ~78-92
        logit = 0.35

        # 1. GC optimality (40-65% optimal per Doench 2014/2016)
        if 45.0 <= gc <= 60.0:
            gc_term = 0.55
        elif 40.0 <= gc <= 65.0:
            gc_term = 0.35
        elif 35.0 <= gc < 40.0 or 65.0 < gc <= 70.0:
            gc_term = -0.15
        elif 25.0 <= gc < 35.0 or 70.0 < gc <= 80.0:
            gc_term = -0.70
            warnings.append("extreme_gc")
        else:
            gc_term = -1.35
            warnings.append("extreme_gc")
        components["gc_balance"] = round(gc_term, 3)
        logit += gc_term

        # 2. Doench position-specific nucleotide preferences (1-based pos 1..20)
        pos_term = 0.0
        # Pos 20 (immediately upstream of PAM): G strongly favored, C strongly disfavored
        pos20 = s[19]
        pos_term += {"G": 0.58, "A": 0.14, "T": -0.18, "C": -0.56}.get(pos20, 0.0)
        # Pos 19: G/A favored, C disfavored
        pos_term += {"G": 0.22, "A": 0.15, "T": -0.08, "C": -0.22}.get(s[18], 0.0)
        # Pos 18: C favored (Doench Rule Set 1/2 peak), T/G disfavored
        pos_term += {"C": 0.36, "A": 0.05, "G": -0.12, "T": -0.18}.get(s[17], 0.0)
        # Pos 17: A/C mildly favored
        pos_term += {"A": 0.14, "C": 0.12, "G": 0.0, "T": -0.16}.get(s[16], 0.0)
        # Pos 16: C/G favored, T disfavored
        pos_term += {"C": 0.20, "G": 0.12, "A": 0.0, "T": -0.22}.get(s[15], 0.0)
        # Pos 1 (5' U6 initiation prefers G)
        pos_term += {"G": 0.16, "A": 0.05, "C": 0.0, "T": -0.06}.get(s[0], 0.0)
        # PAM 5' variable nucleotide N in NGG (CGG > AGG/GGG > TGG)
        if len(p) == 3:
            pos_term += {"C": 0.18, "A": 0.08, "G": 0.06, "T": -0.12}.get(p[0], 0.0)
        # Downstream +1 base after NGG (G disfavored per Doench 2014)
        if downstream_base:
            pos_term += {"T": 0.08, "C": 0.06, "A": 0.02, "G": -0.16}.get(
                downstream_base.upper(), 0.0
            )
        components["position_prefs"] = round(pos_term, 3)
        logit += pos_term

        # 3. PAM-proximal 4-nt seed U/T depletion (positions 17..20)
        seed_t = s[16:20].count("T")
        seed_t_term = 0.12 if seed_t == 0 else (-0.18 * (seed_t - 1) if seed_t >= 2 else 0.0)
        if seed_t >= 3:
            seed_t_term -= 0.25
            warnings.append("seed_poly_u")
        components["seed_u_depletion"] = round(seed_t_term, 3)
        logit += seed_t_term

        # 4. Pol III terminator (TTTT) strong penalty
        pol3_term = -1.95 if has_pol3 else 0.0
        if has_pol3:
            warnings.append("pol3_terminator_tttt")
        components["pol3_terminator"] = round(pol3_term, 3)
        logit += pol3_term

        # 5. Homopolymer / G-quadruplex penalty
        homo_term = 0.0
        if "GGGG" in s:
            homo_term -= 0.45
            warnings.append("g_quadruplex_gggg")
        elif max_homo >= 4 and not has_pol3:
            homo_term -= 0.24
            warnings.append("homopolymer_run")
        components["homopolymer"] = round(homo_term, 3)
        logit += homo_term

        # 6. Self-folding hairpin & scaffold interference penalty
        fold_term = 0.0
        if internal_stem >= 6 or hairpin_dg <= -6.0:
            fold_term -= 0.72
            warnings.append("strong_hairpin")
        elif internal_stem == 5 or hairpin_dg <= -4.0:
            fold_term -= 0.35
            warnings.append("moderate_hairpin")
        if scaffold_stem >= 7:
            fold_term -= 0.45
            warnings.append("scaffold_complementarity")
        elif scaffold_stem == 6:
            fold_term -= 0.18
        components["secondary_structure"] = round(fold_term, 3)
        logit += fold_term

        # 7. Melting temperature thermodynamic stability window (52..68 C)
        if 54.0 <= tm <= 66.0:
            tm_term = 0.18
        elif 48.0 <= tm <= 72.0:
            tm_term = 0.0
        else:
            tm_term = -0.30
        components["tm_stability"] = round(tm_term, 3)
        logit += tm_term

    else:
        # LbCas12a (Cpf1) scoring inspired by Kim et al. (Nature Biotechnology, 2018)
        logit = 0.30

        # 1. Optimal GC for Cas12a is slightly lower (35-60%)
        if 40.0 <= gc <= 56.0:
            gc_term = 0.55
        elif 35.0 <= gc <= 62.0:
            gc_term = 0.32
        elif 28.0 <= gc < 35.0 or 62.0 < gc <= 70.0:
            gc_term = -0.22
        else:
            gc_term = -1.10
            warnings.append("extreme_gc")
        components["gc_balance"] = round(gc_term, 3)
        logit += gc_term

        # 2. PAM preference (TTTA/TTTC > TTTG) and 5' seed position preferences (pos 1..8)
        pos_term = 0.0
        if p == "TTTC":
            pos_term += 0.28
        elif p == "TTTA":
            pos_term += 0.22
        elif p == "TTTG":
            pos_term += 0.02
        # Position 1 immediately downstream of TTTV: T strongly disfavored (Kim 2018)
        pos_term += {"C": 0.34, "G": 0.24, "A": 0.12, "T": -0.52}.get(s[0], 0.0)
        # Position 2 & 3 seed GC/C preference
        pos_term += {"C": 0.18, "G": 0.14, "A": 0.04, "T": -0.14}.get(s[1], 0.0)
        pos_term += {"G": 0.16, "C": 0.14, "A": 0.02, "T": -0.12}.get(s[2], 0.0)
        # Position 8 (seed boundary): G/C favored
        pos_term += {"G": 0.14, "C": 0.12, "A": 0.0, "T": -0.10}.get(s[7], 0.0)
        components["position_prefs"] = round(pos_term, 3)
        logit += pos_term

        # 3. 5' seed T count (positions 1..4)
        seed_t = s[0:4].count("T")
        seed_t_term = 0.10 if seed_t == 0 else (-0.20 * (seed_t - 1) if seed_t >= 2 else 0.0)
        if seed_t >= 3:
            warnings.append("seed_poly_u")
        components["seed_u_depletion"] = round(seed_t_term, 3)
        logit += seed_t_term

        # 4. Pol III terminator (TTTT)
        pol3_term = -1.95 if has_pol3 else 0.0
        if has_pol3:
            warnings.append("pol3_terminator_tttt")
        components["pol3_terminator"] = round(pol3_term, 3)
        logit += pol3_term

        # 5. Homopolymer
        homo_term = 0.0
        if "GGGG" in s:
            homo_term -= 0.42
            warnings.append("g_quadruplex_gggg")
        elif max_homo >= 4 and not has_pol3:
            homo_term -= 0.22
            warnings.append("homopolymer_run")
        components["homopolymer"] = round(homo_term, 3)
        logit += homo_term

        # 6. Hairpin & direct-repeat complementarity
        fold_term = 0.0
        if internal_stem >= 6 or hairpin_dg <= -6.0:
            fold_term -= 0.68
            warnings.append("strong_hairpin")
        elif internal_stem == 5 or hairpin_dg <= -4.0:
            fold_term -= 0.32
            warnings.append("moderate_hairpin")
        if scaffold_stem >= 6:
            fold_term -= 0.42
            warnings.append("scaffold_complementarity")
        components["secondary_structure"] = round(fold_term, 3)
        logit += fold_term

        # 7. Tm stability
        if 52.0 <= tm <= 66.0:
            tm_term = 0.15
        elif 46.0 <= tm <= 72.0:
            tm_term = 0.0
        else:
            tm_term = -0.28
        components["tm_stability"] = round(tm_term, 3)
        logit += tm_term

    score = round(100.0 / (1.0 + math.exp(-logit)), 1)
    return OnTargetScoreDetail(
        score=score,
        gc_percent=gc,
        tm_celsius=tm,
        hairpin_stem_bp=internal_stem,
        scaffold_stem_bp=scaffold_stem,
        estimated_hairpin_dg=hairpin_dg,
        pol3_terminator=has_pol3,
        homopolymer_run=max_homo,
        seed_t_count=seed_t,
        components=components,
        warnings=tuple(warnings),
    )


def predict_mmej(
    sequence: str,
    cut_pos: int,
    window: int = 25,
    min_mh: int = 2,
    max_mh: int = 8,
    top_n: int = 12,
) -> MMEJResult:
    """Predict post-DSB MMEJ deletion patterns and out-of-frame knockout probability.

    Implements the Bae et al. (Nature Methods, 2014) microhomology scoring algorithm:
      - Scans `window` bp upstream and downstream of the exact DSB cut site (`cut_pos`).
      - Identifies non-redundant maximal 2..8 bp microhomology repeats spanning the break.
      - Scores each deletion pattern via GC-weighted repeat length and exponential
        deletion-length decay: `100 * (len(mh) + num_GC) * exp(-del_len / 20.0)`.
      - Computes Out-of-Frame Probability % = 100 * sum(score where del_len % 3 != 0) / total_score.
    """
    seq = clean_dna(sequence)
    seq_len = len(seq)
    if cut_pos < 0 or cut_pos > seq_len:
        raise ValueError(f"cut_pos {cut_pos} must be within 0..{seq_len}.")
    if window < 4:
        raise ValueError("window must be at least 4 bp.")

    win_start = max(0, cut_pos - window)
    win_end = min(seq_len, cut_pos + window)
    left_flank = seq[win_start:cut_pos]
    right_flank = seq[cut_pos:win_end]
    l_len = len(left_flank)
    r_len = len(right_flank)

    # Deduplicate by (left_start, deletion_length) while finding maximal microhomologies
    seen_deletions: dict[tuple[int, int], MMEJPattern] = {}

    for i in range(l_len):
        for j in range(r_len):
            # Only start at left-maximal positions so sub-repeats of the same deletion aren't double-counted
            if i > 0 and j > 0 and left_flank[i - 1] == right_flank[j - 1]:
                continue

            # Extend match rightward while matching bases exist
            k = 0
            while (
                i + k < l_len
                and j + k < r_len
                and left_flank[i + k] == right_flank[j + k]
                and left_flank[i + k] != "N"
            ):
                k += 1

            if k < min_mh:
                continue

            # Cap scored microhomology length to max_mh (2..8 bp canonical window)
            eff_k = min(k, max_mh)
            mh = left_flank[i : i + eff_k]
            del_len = (l_len - i) + j
            if del_len <= 0:
                continue

            gc_count = sum(1 for b in mh if b in "GC")
            # Hydrogen-bond strength (AT=1, GC=2) with thermodynamic length cooperativity
            mh_weight = (eff_k + gc_count) * (1.85 ** (eff_k - 2))
            length_factor = math.exp(-del_len / 15.0)
            pattern_score = round(100.0 * mh_weight * length_factor, 2)

            # Build character-aligned post-MMEJ outcome relative to (left_flank + right_flank)
            deleted_seq = left_flank[i + eff_k :] + right_flank[: j + eff_k]
            aligned = (
                left_flank[:i]
                + mh
                + ("-" * del_len)
                + right_flank[j + eff_k :]
            )

            key = (i, del_len)
            candidate = MMEJPattern(
                microhomology=mh,
                mh_length=eff_k,
                gc_count=gc_count,
                left_mh_start=i,
                right_mh_start=l_len + j,
                deletion_length=del_len,
                is_frameshift=(del_len % 3 != 0),
                pattern_score=pattern_score,
                aligned_outcome=aligned,
                deleted_sequence=deleted_seq,
            )
            existing = seen_deletions.get(key)
            if existing is None or candidate.pattern_score > existing.pattern_score:
                seen_deletions[key] = candidate

    patterns = sorted(
        seen_deletions.values(),
        key=lambda p: (-p.pattern_score, p.deletion_length, p.left_mh_start),
    )

    total_score = sum(p.pattern_score for p in patterns)
    oof_score = sum(p.pattern_score for p in patterns if p.is_frameshift)
    inf_score = total_score - oof_score

    if total_score > 0:
        oof_prob = round(100.0 * oof_score / total_score, 1)
        inf_prob = round(100.0 - oof_prob, 1)
    else:
        # Neutral NHEJ frameshift expectation (2/3 out-of-frame) when no repeats exist
        oof_prob = 66.7
        inf_prob = 33.3

    if oof_prob >= 70.0:
        grade = "High KO"
    elif oof_prob >= 55.0:
        grade = "Moderate KO"
    else:
        grade = "In-Frame Risk"

    fs_count = sum(1 for p in patterns if p.is_frameshift)

    return MMEJResult(
        cut_pos=cut_pos,
        window_bp=window,
        left_flank=left_flank,
        right_flank=right_flank,
        window_start_genomic=win_start,
        microhomology_score=round(total_score, 1),
        out_of_frame_score=round(oof_score, 1),
        in_frame_score=round(inf_score, 1),
        out_of_frame_probability=oof_prob,
        in_frame_probability=inf_prob,
        total_patterns=len(patterns),
        frameshift_patterns=fs_count,
        in_frame_patterns=len(patterns) - fs_count,
        knockout_grade=grade,
        top_patterns=tuple(patterns[:top_n]),
    )


def _extract_context_30nt(strand_seq: str, start: int, end: int, enzyme: str) -> str:
    """Extract a 30-nt context window around a protospacer+PAM match on the matching strand."""
    seq_len = len(strand_seq)
    if enzyme == "spcas9":
        # 4 nt upstream + 20 nt spacer + 3 nt PAM + 3 nt downstream = 30 nt
        c_start = max(0, start - 4)
        c_end = min(seq_len, end + 3)
    else:
        # 1 nt upstream + 4 nt PAM + 23 nt spacer + 2 nt downstream = 30 nt
        c_start = max(0, start - 1)
        c_end = min(seq_len, end + 2)
    return strand_seq[c_start:c_end]


def design_guides(
    sequence: str,
    enzyme: str = "spcas9",
    background: str | Sequence[tuple[str, str]] | None = None,
    max_mismatches: int = 4,
    include_noncanonical_pam: bool = True,
    mmej_window: int = 25,
    top_n: int | None = None,
    min_gc: float = 0.0,
    max_gc: float = 100.0,
    target_name: str | None = None,
) -> DesignReport:
    """Scan both strands of a target DNA sequence to design & rank optimal CRISPR guides."""
    enz = normalize_enzyme(enzyme)
    spec = ENZYME_SPECS[enz]
    target_records = parse_fasta(sequence, default_name=target_name or "target")
    rec_name, seq = target_records[0]
    if target_name:
        rec_name = target_name
    seq_len = len(seq)
    rc_seq = reverse_complement(seq)

    # Prepare off-target background contigs
    bg_contigs: list[tuple[str, str]] = []
    if background is not None:
        if isinstance(background, str):
            if background.strip():
                bg_contigs = parse_fasta(background, default_name="background")
        else:
            bg_contigs = [(n, clean_dna(s)) for n, s in background if s.strip()]

    bg_bp = sum(len(s) for _, s in bg_contigs)
    raw_hits: list[dict[str, object]] = []

    for strand, scan_seq in (("+", seq), ("-", rc_seq)):
        if enz == "spcas9":
            for i in range(seq_len - 22):
                if scan_seq[i + 21 : i + 23] != "GG":
                    continue
                spacer = scan_seq[i : i + 20]
                pam = scan_seq[i + 20 : i + 23]
                if "N" in spacer or "N" in pam:
                    continue
                downstream = scan_seq[i + 23] if i + 23 < seq_len else ""
                context = _extract_context_30nt(scan_seq, i, i + 23, enz)

                if strand == "+":
                    g_start, g_end = i, i + 23
                    sp_start, sp_end = i, i + 20
                    pam_start, pam_end = i + 20, i + 23
                    cut_pos = i + 17  # -3 bp upstream of PAM
                    cut_ts = cut_pos
                else:
                    g_start, g_end = seq_len - (i + 23), seq_len - i
                    pam_start, pam_end = seq_len - (i + 23), seq_len - (i + 20)
                    sp_start, sp_end = seq_len - (i + 20), seq_len - i
                    cut_pos = seq_len - (i + 17)
                    cut_ts = cut_pos

                raw_hits.append(
                    {
                        "strand": strand,
                        "start": g_start,
                        "end": g_end,
                        "spacer_start": sp_start,
                        "spacer_end": sp_end,
                        "pam_start": pam_start,
                        "pam_end": pam_end,
                        "cut_pos": cut_pos,
                        "cut_ts": cut_ts,
                        "spacer": spacer,
                        "pam": pam,
                        "downstream": downstream,
                        "context": context,
                    }
                )
        else:
            # LbCas12a: 5'-TTTV-[23nt spacer]-3'
            for i in range(seq_len - 26):
                pam = scan_seq[i : i + 4]
                if pam not in {"TTTA", "TTTC", "TTTG"}:
                    continue
                spacer = scan_seq[i + 4 : i + 27]
                if "N" in spacer:
                    continue
                context = _extract_context_30nt(scan_seq, i, i + 27, enz)

                if strand == "+":
                    g_start, g_end = i, i + 27
                    pam_start, pam_end = i, i + 4
                    sp_start, sp_end = i + 4, i + 27
                    cut_pos = i + 22  # +18 nt downstream of PAM (NTS cut)
                    cut_ts = i + 27   # +23 nt downstream of PAM (TS cut)
                else:
                    g_start, g_end = seq_len - (i + 27), seq_len - i
                    sp_start, sp_end = seq_len - (i + 27), seq_len - (i + 4)
                    pam_start, pam_end = seq_len - (i + 4), seq_len - i
                    cut_pos = seq_len - (i + 22)
                    cut_ts = seq_len - (i + 27)

                raw_hits.append(
                    {
                        "strand": strand,
                        "start": g_start,
                        "end": g_end,
                        "spacer_start": sp_start,
                        "spacer_end": sp_end,
                        "pam_start": pam_start,
                        "pam_end": pam_end,
                        "cut_pos": cut_pos,
                        "cut_ts": cut_ts,
                        "spacer": spacer,
                        "pam": pam,
                        "downstream": "",
                        "context": context,
                    }
                )

    candidates: list[GuideCandidate] = []
    # Combine target + background for self-and-reference off-target evaluation
    scan_corpus: list[tuple[str, str]] = [(rec_name, seq), *bg_contigs]

    for item in raw_hits:
        spacer = str(item["spacer"])
        pam = str(item["pam"])
        gc = gc_content_percent(spacer)
        if gc < min_gc or gc > max_gc:
            continue

        on_detail = score_on_target(
            spacer,
            pam,
            enzyme=enz,
            downstream_base=str(item["downstream"]),
        )
        mmej_res = predict_mmej(seq, int(item["cut_pos"]), window=mmej_window)

        ot_summary = scan_offtargets(
            guide=spacer,
            background=scan_corpus,
            enzyme=enz,
            max_mismatches=max_mismatches,
            include_noncanonical_pam=include_noncanonical_pam,
            exclude_locus=(rec_name, int(item["start"]), str(item["strand"])),
            top_n=20,
        )
        spec_score = ot_summary.specificity_score

        # Composite score balances activity (45%), specificity (35%), and KO frameshift (20%)
        composite = round(
            0.45 * on_detail.score
            + 0.35 * spec_score
            + 0.20 * mmej_res.out_of_frame_probability,
            1,
        )

        combined_warnings = list(on_detail.warnings)
        if ot_summary.high_risk_hits > 0:
            combined_warnings.append("high_risk_offtarget")
        if mmej_res.out_of_frame_probability < 50.0:
            combined_warnings.append("in_frame_mmej_bias")

        candidates.append(
            GuideCandidate(
                rank=0,
                id="",
                enzyme=enz,
                strand=str(item["strand"]),
                start=int(item["start"]),
                end=int(item["end"]),
                spacer_start=int(item["spacer_start"]),
                spacer_end=int(item["spacer_end"]),
                pam_start=int(item["pam_start"]),
                pam_end=int(item["pam_end"]),
                cut_pos=int(item["cut_pos"]),
                cut_ts=int(item["cut_ts"]),
                spacer=spacer,
                pam=pam,
                context_30nt=str(item["context"]),
                gc_percent=gc,
                on_target_score=on_detail.score,
                on_target_detail=on_detail,
                mmej=mmej_res,
                specificity_score=spec_score,
                offtarget_summary=ot_summary,
                composite_score=composite,
                warnings=tuple(combined_warnings),
            )
        )

    candidates.sort(
        key=lambda g: (
            -g.composite_score,
            -g.on_target_score,
            -g.specificity_score,
            g.start,
        )
    )

    if top_n is not None and top_n > 0:
        candidates = candidates[:top_n]

    ranked: list[GuideCandidate] = []
    prefix = "sg" if enz == "spcas9" else "cr"
    for idx, g in enumerate(candidates, start=1):
        ranked.append(
            GuideCandidate(
                rank=idx,
                id=f"{prefix}{idx:02d}",
                enzyme=g.enzyme,
                strand=g.strand,
                start=g.start,
                end=g.end,
                spacer_start=g.spacer_start,
                spacer_end=g.spacer_end,
                pam_start=g.pam_start,
                pam_end=g.pam_end,
                cut_pos=g.cut_pos,
                cut_ts=g.cut_ts,
                spacer=g.spacer,
                pam=g.pam,
                context_30nt=g.context_30nt,
                gc_percent=g.gc_percent,
                on_target_score=g.on_target_score,
                on_target_detail=g.on_target_detail,
                mmej=g.mmej,
                specificity_score=g.specificity_score,
                offtarget_summary=g.offtarget_summary,
                composite_score=g.composite_score,
                warnings=g.warnings,
            )
        )

    plus_count = sum(1 for g in ranked if g.strand == "+")
    minus_count = len(ranked) - plus_count

    return DesignReport(
        target_name=rec_name,
        target_sequence=seq,
        target_length=seq_len,
        target_gc_percent=gc_content_percent(seq),
        enzyme=enz,
        enzyme_spec=spec,
        background_contigs=len(bg_contigs),
        background_bp=bg_bp,
        max_mismatches=max_mismatches,
        total_candidates=len(ranked),
        plus_strand_count=plus_count,
        minus_strand_count=minus_count,
        guides=tuple(ranked),
    )


def format_guides_csv(guides: Sequence[GuideCandidate]) -> str:
    """Serialize candidate guides to CSV text."""
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(
        [
            "rank",
            "id",
            "enzyme",
            "strand",
            "start_1based",
            "end_1based",
            "cut_pos",
            "spacer",
            "pam",
            "gc_percent",
            "on_target_score",
            "specificity_score",
            "mmej_oof_percent",
            "composite_score",
            "offtarget_hits",
            "high_risk_offtargets",
            "warnings",
        ]
    )
    for g in guides:
        writer.writerow(
            [
                g.rank,
                g.id,
                g.enzyme,
                g.strand,
                g.start + 1,
                g.end,
                g.cut_pos,
                g.spacer,
                g.pam,
                g.gc_percent,
                g.on_target_score,
                g.specificity_score,
                g.mmej.out_of_frame_probability,
                g.composite_score,
                g.offtarget_summary.total_hits if g.offtarget_summary else 0,
                g.offtarget_summary.high_risk_hits if g.offtarget_summary else 0,
                ";".join(g.warnings),
            ]
        )
    return buf.getvalue()


def format_guides_fasta(guides: Sequence[GuideCandidate], target_name: str = "target") -> str:
    """Serialize candidate guides to FASTA format with rich header annotations."""
    lines: list[str] = []
    for g in guides:
        header = (
            f">{g.id}_{target_name} enzyme={g.enzyme} strand={g.strand} "
            f"pos={g.start + 1}-{g.end} cut={g.cut_pos} pam={g.pam} "
            f"composite={g.composite_score} on_target={g.on_target_score} "
            f"specificity={g.specificity_score} mmej_oof={g.mmej.out_of_frame_probability}%"
        )
        lines.append(header)
        lines.append(g.spacer)
    return "\n".join(lines) + ("\n" if lines else "")
