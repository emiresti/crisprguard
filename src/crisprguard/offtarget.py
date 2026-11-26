"""PAM-proximal seed-weighted CFD off-target kinetics engine for SpCas9 and LbCas12a.

Implements a thermodynamic Cutting Frequency Determination (CFD) scoring model inspired
by Doench et al. (Nature Biotechnology, 2016) for SpCas9 (5'-20nt-NGG-3') and
Kim et al. / Kleinstiver et al. for LbCas12a (5'-TTTV-23nt-3').
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import re
from typing import Iterable, Sequence

DNA_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)
VALID_DNA_RE = re.compile(r"^[ACGTN]+$")

# SpCas9: 20 nt spacer, PAM is 3' (positions 11..20 are PAM-proximal seed; 13..20 core seed)
SPCAS9_POS_WEIGHTS: tuple[float, ...] = (
    0.86, 0.84, 0.80, 0.76, 0.72,  # 1..5   PAM-distal
    0.67, 0.62, 0.56, 0.50, 0.44,  # 6..10  Mid-spacer
    0.36, 0.29, 0.23, 0.18, 0.14,  # 11..15 Seed entry
    0.11, 0.085, 0.065, 0.050, 0.040,  # 16..20 Core PAM-proximal seed
)

# LbCas12a: 23 nt spacer, PAM is 5' (positions 1..8 are PAM-proximal seed)
CAS12A_POS_WEIGHTS: tuple[float, ...] = (
    0.045, 0.055, 0.070, 0.090, 0.115, 0.145, 0.180, 0.220,  # 1..8   5' PAM-proximal seed
    0.270, 0.325, 0.385, 0.445, 0.505, 0.560, 0.610, 0.660,  # 9..16  Trunk
    0.705, 0.745, 0.785, 0.820, 0.850, 0.875, 0.895,         # 17..23 3' PAM-distal
)

# RNA:target-DNA mismatch kinetic odds multipliers (Boltzmann / logistic shift)
# Note: guide base is non-target strand (T -> U in RNA); target DNA strand is complement of protospacer.
# E.g. guide='G' (rG) vs off-target protospacer='A' (target DNA='T') -> rG:dT wobble.
MISMATCH_ODDS_MULTIPLIERS: dict[str, float] = {
    "rG:dT": 4.20,  # Canonical wobble — highest mismatch tolerance in Cas9/Cas12a
    "rU:dG": 2.05,  # Secondary wobble
    "rA:dC": 1.18,  # Purine-pyrimidine transition mismatch
    "rC:dA": 1.12,  # Pyrimidine-purine transition mismatch
    "rU:dC": 0.86,  # Pyrimidine-pyrimidine
    "rC:dT": 0.82,  # Pyrimidine-pyrimidine
    "rU:dT": 0.78,  # Pyrimidine-pyrimidine
    "rA:dA": 0.64,  # Purine-purine transversion
    "rA:dG": 0.60,  # Purine-purine transversion
    "rG:dA": 0.58,  # Purine-purine transversion
    "rG:dG": 0.54,  # Bulky purine-purine transversion
    "rC:dC": 0.50,  # Rigid pyrimidine-pyrimidine clash
}

WOBBLE_PAIRS: frozenset[str] = frozenset({"rG:dT", "rU:dG"})

# PAM activity scores (Doench et al. 2016 Table S19 for SpCas9; Kleinstiver/Kim for LbCas12a)
SPCAS9_PAM_SCORES: dict[str, float] = {
    "GG": 1.000,  # NGG canonical
    "AG": 0.259,  # NAG primary non-canonical
    "GA": 0.069,  # NGA secondary non-canonical
    "GC": 0.022,  # NGC weak non-canonical
    "GT": 0.016,  # NGT weak non-canonical
    "CG": 0.011,  # NCG weak non-canonical
}

CAS12A_PAM_SCORES: dict[str, float] = {
    "TTTA": 1.000,
    "TTTC": 1.000,
    "TTTG": 1.000,
    "TTTT": 0.200,
    "CTTA": 0.180,
    "CTTC": 0.180,
    "CTTG": 0.180,
    "TCTA": 0.120,
    "TCTC": 0.120,
    "TCTG": 0.120,
    "TTCA": 0.100,
    "TTCC": 0.100,
    "TTCG": 0.100,
}


def normalize_enzyme(enzyme: str) -> str:
    """Normalize enzyme alias to canonical 'spcas9' or 'cas12a'."""
    key = (enzyme or "spcas9").strip().lower().replace("-", "").replace("_", "")
    if key in {"spcas9", "cas9", "ngg"}:
        return "spcas9"
    if key in {"cas12a", "lbcas12a", "ascas12a", "cpf1", "tttv"}:
        return "cas12a"
    raise ValueError(f"Unsupported enzyme '{enzyme}'. Expected 'spcas9' or 'cas12a'.")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def clean_dna(raw: str) -> str:
    """Clean a raw DNA string or single FASTA entry into uppercase A/C/G/T/N."""
    if not raw or not raw.strip():
        raise ValueError("DNA sequence input is empty.")
    lines: list[str] = []
    for line in raw.strip().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(">") or stripped.startswith(";"):
            continue
        lines.append(stripped)
    cleaned = re.sub(r"[\s0-9]+", "", "".join(lines)).upper().replace("U", "T")
    if not cleaned:
        raise ValueError("No nucleotide bases found after parsing input.")
    if not VALID_DNA_RE.match(cleaned):
        invalid = sorted(set(cleaned) - {"A", "C", "G", "T", "N"})
        raise ValueError(
            f"Invalid nucleotide character(s) found: {', '.join(invalid)}. "
            "Allowed bases: A, C, G, T, U, N."
        )
    return cleaned


def parse_fasta(raw: str, default_name: str = "sequence") -> list[tuple[str, str]]:
    """Parse raw DNA or multi-FASTA text into a list of (record_name, cleaned_seq)."""
    if not raw or not raw.strip():
        raise ValueError("FASTA / DNA input is empty.")

    records: list[tuple[str, list[str]]] = []
    current_name: str | None = None
    current_chunks: list[str] = []

    for line in raw.strip().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith(">"):
            if current_name is not None or current_chunks:
                records.append((current_name or default_name, current_chunks))
            header = stripped[1:].strip()
            current_name = header.split()[0] if header else f"{default_name}_{len(records) + 1}"
            current_chunks = []
        else:
            if current_name is None:
                current_name = default_name
            current_chunks.append(stripped)

    if current_name is not None or current_chunks:
        records.append((current_name or default_name, current_chunks))

    parsed: list[tuple[str, str]] = []
    for idx, (name, chunks) in enumerate(records, start=1):
        joined = re.sub(r"[\s0-9]+", "", "".join(chunks)).upper().replace("U", "T")
        if not joined:
            continue
        if not VALID_DNA_RE.match(joined):
            invalid = sorted(set(joined) - {"A", "C", "G", "T", "N"})
            raise ValueError(
                f"Invalid nucleotide character(s) in record '{name}': {', '.join(invalid)}."
            )
        parsed.append((name or f"{default_name}_{idx}", joined))

    if not parsed:
        raise ValueError("No valid DNA sequence records found in input.")
    return parsed


@dataclass(frozen=True)
class MismatchDetail:
    """Single base-pair mismatch between guide RNA and off-target DNA locus."""

    position: int  # 1-based position along spacer (1..20 SpCas9, 1..23 Cas12a)
    guide_base: str  # Non-target RNA/DNA base (A/C/G/T)
    protospacer_base: str  # Off-target non-target strand DNA base
    target_dna_base: str  # Complement on target DNA strand paired with RNA
    pair: str  # e.g. 'rG:dT'
    is_seed: bool
    is_wobble: bool
    cfd_factor: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class OffTargetHit:
    """Candidate off-target genomic locus with CFD kinetic score."""

    contig: str
    start: int  # 0-based start on + reference strand
    end: int  # 0-based exclusive end on + reference strand
    strand: str  # '+' or '-'
    protospacer: str
    pam: str
    pam_score: float
    alignment: str  # '.' for match, mismatched protospacer base otherwise
    mismatches: int
    seed_mismatches: int
    wobble_mismatches: int
    cfd_score: float
    risk_tier: str  # 'critical', 'high', 'moderate', 'low'
    mismatch_details: tuple[MismatchDetail, ...]

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["mismatch_details"] = [m.to_dict() for m in self.mismatch_details]
        return data


@dataclass(frozen=True)
class OffTargetSummary:
    """Aggregate CFD specificity evaluation for a guide RNA."""

    guide: str
    enzyme: str
    max_mismatches: int
    specificity_score: float  # 0..100 MIT/CFD aggregate score: 100 / (1 + sum(cfd))
    total_cfd_sum: float
    total_hits: int
    high_risk_hits: int  # CFD >= 0.20
    canonical_pam_hits: int
    noncanonical_pam_hits: int
    scanned_bp: int
    scanned_contigs: int
    top_hits: tuple[OffTargetHit, ...]

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["top_hits"] = [h.to_dict() for h in self.top_hits]
        return data


def is_seed_position(position_1based: int, enzyme: str = "spcas9") -> bool:
    """Return True if the 1-based spacer position falls in the PAM-proximal seed region."""
    enz = normalize_enzyme(enzyme)
    if enz == "spcas9":
        return 11 <= position_1based <= 20
    return 1 <= position_1based <= 8


def rna_dna_pair_label(guide_base: str, protospacer_base: str) -> tuple[str, str]:
    """Return (target_dna_base, 'rX:dY') for a guide RNA base vs non-target DNA base."""
    r_base = "U" if guide_base == "T" else guide_base
    target_dna_base = protospacer_base.translate(DNA_COMPLEMENT)
    return target_dna_base, f"r{r_base}:d{target_dna_base}"


def score_cfd_mismatch(
    position_1based: int,
    guide_base: str,
    protospacer_base: str,
    enzyme: str = "spcas9",
) -> MismatchDetail:
    """Compute position-specific CFD tolerance factor for a single RNA:DNA mismatch."""
    enz = normalize_enzyme(enzyme)
    weights = SPCAS9_POS_WEIGHTS if enz == "spcas9" else CAS12A_POS_WEIGHTS
    if not (1 <= position_1based <= len(weights)):
        raise ValueError(
            f"Position {position_1based} out of range for {enz} (1..{len(weights)})."
        )

    g_base = guide_base.upper().replace("U", "T")
    p_base = protospacer_base.upper().replace("U", "T")
    target_dna, pair = rna_dna_pair_label(g_base, p_base)
    in_seed = is_seed_position(position_1based, enz)
    wobble = pair in WOBBLE_PAIRS

    if g_base == p_base:
        factor = 1.0
    elif g_base == "N" or p_base == "N":
        factor = round(weights[position_1based - 1] * 0.5, 4)
    else:
        base_p = weights[position_1based - 1]
        mult = MISMATCH_ODDS_MULTIPLIERS.get(pair, 0.70)
        # Boltzmann odds-ratio transformation keeps strict monotonicity in position & pair type
        factor = (base_p * mult) / (base_p * mult + (1.0 - base_p))
        factor = round(factor, 4)

    return MismatchDetail(
        position=position_1based,
        guide_base=g_base,
        protospacer_base=p_base,
        target_dna_base=target_dna,
        pair=pair,
        is_seed=in_seed,
        is_wobble=wobble,
        cfd_factor=factor,
    )


def pam_activity_score(
    pam: str,
    enzyme: str = "spcas9",
    include_noncanonical: bool = True,
) -> float:
    """Return CFD PAM recognition weight (0.0 if inactive)."""
    enz = normalize_enzyme(enzyme)
    p = pam.upper()
    if enz == "spcas9":
        if len(p) != 3 or "N" in p[1:]:
            return 0.0
        dinuc = p[1:3]
        if dinuc == "GG":
            return SPCAS9_PAM_SCORES["GG"]
        if include_noncanonical:
            return SPCAS9_PAM_SCORES.get(dinuc, 0.0)
        return 0.0
    else:
        if len(p) != 4 or "N" in p:
            return 0.0
        if p in {"TTTA", "TTTC", "TTTG"}:
            return 1.0
        if include_noncanonical:
            return CAS12A_PAM_SCORES.get(p, 0.0)
        return 0.0


def score_offtarget_cfd(
    guide: str,
    protospacer: str,
    pam: str,
    enzyme: str = "spcas9",
    include_noncanonical_pam: bool = True,
) -> tuple[float, float, str, tuple[MismatchDetail, ...]]:
    """Score a single candidate off-target site against a guide RNA.

    Returns (cfd_score, pam_score, alignment_string, mismatch_details).
    """
    enz = normalize_enzyme(enzyme)
    g = clean_dna(guide)
    p = clean_dna(protospacer)
    expected_len = 20 if enz == "spcas9" else 23
    if len(g) != expected_len:
        raise ValueError(f"{enz} guide must be exactly {expected_len} nt (got {len(g)}).")
    if len(p) != expected_len:
        raise ValueError(
            f"{enz} protospacer must be exactly {expected_len} nt (got {len(p)})."
        )

    p_score = pam_activity_score(pam, enz, include_noncanonical=include_noncanonical_pam)
    if p_score <= 0.0:
        return 0.0, 0.0, "." * expected_len, ()

    score = p_score
    align_chars: list[str] = []
    details: list[MismatchDetail] = []
    seed_positions: list[int] = []

    for idx in range(expected_len):
        gb = g[idx]
        pb = p[idx]
        if gb == pb:
            align_chars.append(".")
        else:
            align_chars.append(pb)
            detail = score_cfd_mismatch(idx + 1, gb, pb, enz)
            details.append(detail)
            score *= detail.cfd_factor
            if detail.is_seed:
                seed_positions.append(idx + 1)

    # Consecutive mismatches inside the PAM-proximal seed strongly destabilize R-loop propagation
    adjacent_seed_pairs = sum(
        1 for a, b in zip(seed_positions, seed_positions[1:]) if b == a + 1
    )
    if adjacent_seed_pairs > 0:
        score *= 0.82 ** adjacent_seed_pairs

    return round(score, 5), round(p_score, 4), "".join(align_chars), tuple(details)


def _classify_risk_tier(cfd_score: float) -> str:
    if cfd_score >= 0.50:
        return "critical"
    if cfd_score >= 0.20:
        return "high"
    if cfd_score >= 0.05:
        return "moderate"
    return "low"


def _iter_strand_sites(
    seq: str,
    enzyme: str,
    include_noncanonical_pam: bool,
) -> Iterable[tuple[int, int, str, str, float]]:
    """Yield (start, end, protospacer, pam, pam_score) along a single 5'->3' strand."""
    length = len(seq)
    if enzyme == "spcas9":
        if length < 23:
            return
        allowed_dinucs = (
            SPCAS9_PAM_SCORES if include_noncanonical_pam else {"GG": 1.0}
        )
        for i in range(length - 22):
            dinuc = seq[i + 21 : i + 23]
            p_score = allowed_dinucs.get(dinuc)
            if p_score is not None:
                protospacer = seq[i : i + 20]
                if "N" in protospacer:
                    continue
                pam = seq[i + 20 : i + 23]
                yield i, i + 23, protospacer, pam, p_score
    else:
        if length < 27:
            return
        allowed_pams = (
            CAS12A_PAM_SCORES
            if include_noncanonical_pam
            else {"TTTA": 1.0, "TTTC": 1.0, "TTTG": 1.0}
        )
        for i in range(length - 26):
            pam = seq[i : i + 4]
            p_score = allowed_pams.get(pam)
            if p_score is not None:
                protospacer = seq[i + 4 : i + 27]
                if "N" in protospacer:
                    continue
                yield i, i + 27, protospacer, pam, p_score


def scan_offtargets(
    guide: str,
    background: str | Sequence[tuple[str, str]],
    enzyme: str = "spcas9",
    max_mismatches: int = 4,
    include_noncanonical_pam: bool = True,
    exclude_locus: tuple[str, int, str] | None = None,
    exclude_one_exact: bool = False,
    top_n: int = 25,
) -> OffTargetSummary:
    """Scan both strands of background sequence(s) for PAM-adjacent off-targets <= max_mismatches.

    Args:
        guide: Candidate 20 nt (SpCas9) or 23 nt (LbCas12a) spacer sequence.
        background: Raw DNA/FASTA text or list of (contig_name, cleaned_seq) tuples.
        enzyme: 'spcas9' or 'cas12a'.
        max_mismatches: Maximum Hamming distance to report (default 4).
        include_noncanonical_pam: Whether to include NAG/NGA/etc. (SpCas9) or TTTT/CTTV (Cas12a).
        exclude_locus: Optional (contig, start, strand) tuple identifying the on-target self site.
        exclude_one_exact: If True, skips the first exact canonical match (0 mismatches, PAM=1.0).
        top_n: Maximum number of top-ranked hits to keep in `top_hits` (all hits count in specificity).
    """
    enz = normalize_enzyme(enzyme)
    g = clean_dna(guide)
    expected_len = 20 if enz == "spcas9" else 23
    if len(g) != expected_len:
        raise ValueError(f"{enz} guide must be {expected_len} nt long (got {len(g)}).")
    if max_mismatches < 0 or max_mismatches > 8:
        raise ValueError("max_mismatches must be between 0 and 8.")

    if isinstance(background, str):
        contigs = parse_fasta(background, default_name="background")
    else:
        contigs = [(name, clean_dna(seq)) for name, seq in background]

    hits: list[OffTargetHit] = []
    skipped_exact = False
    scanned_bp = 0

    for contig_name, seq in contigs:
        seq_len = len(seq)
        scanned_bp += seq_len
        rc_seq = reverse_complement(seq)

        for strand, scan_seq in (("+", seq), ("-", rc_seq)):
            for raw_start, raw_end, protospacer, pam, p_score in _iter_strand_sites(
                scan_seq, enz, include_noncanonical_pam
            ):
                # Fast Hamming distance filter before detailed CFD calculation
                mm_count = 0
                for a, b in zip(g, protospacer):
                    if a != b:
                        mm_count += 1
                        if mm_count > max_mismatches:
                            break
                if mm_count > max_mismatches:
                    continue

                if strand == "+":
                    genomic_start, genomic_end = raw_start, raw_end
                else:
                    genomic_start = seq_len - raw_end
                    genomic_end = seq_len - raw_start

                if exclude_locus is not None and (
                    contig_name == exclude_locus[0]
                    and genomic_start == exclude_locus[1]
                    and strand == exclude_locus[2]
                ):
                    continue

                if exclude_one_exact and not skipped_exact and mm_count == 0 and p_score >= 1.0:
                    skipped_exact = True
                    continue

                cfd, pam_w, alignment, details = score_offtarget_cfd(
                    g,
                    protospacer,
                    pam,
                    enzyme=enz,
                    include_noncanonical_pam=include_noncanonical_pam,
                )
                seed_mm = sum(1 for d in details if d.is_seed)
                wobble_mm = sum(1 for d in details if d.is_wobble)
                hits.append(
                    OffTargetHit(
                        contig=contig_name,
                        start=genomic_start,
                        end=genomic_end,
                        strand=strand,
                        protospacer=protospacer,
                        pam=pam,
                        pam_score=pam_w,
                        alignment=alignment,
                        mismatches=mm_count,
                        seed_mismatches=seed_mm,
                        wobble_mismatches=wobble_mm,
                        cfd_score=cfd,
                        risk_tier=_classify_risk_tier(cfd),
                        mismatch_details=details,
                    )
                )

    hits.sort(key=lambda h: (-h.cfd_score, h.mismatches, h.contig, h.start))
    total_cfd = sum(h.cfd_score for h in hits)
    specificity = round(100.0 / (1.0 + total_cfd), 2)
    high_risk = sum(1 for h in hits if h.cfd_score >= 0.20)
    canonical_hits = sum(1 for h in hits if h.pam_score >= 1.0)
    noncanonical_hits = len(hits) - canonical_hits

    return OffTargetSummary(
        guide=g,
        enzyme=enz,
        max_mismatches=max_mismatches,
        specificity_score=specificity,
        total_cfd_sum=round(total_cfd, 5),
        total_hits=len(hits),
        high_risk_hits=high_risk,
        canonical_pam_hits=canonical_hits,
        noncanonical_pam_hits=noncanonical_hits,
        scanned_bp=scanned_bp,
        scanned_contigs=len(contigs),
        top_hits=tuple(hits[:top_n]),
    )
