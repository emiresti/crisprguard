"""Command-line interface for CrisprGuard."""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import sys
from typing import Sequence

from crisprguard import __version__
from crisprguard.design import (
    DesignReport,
    MMEJResult,
    design_guides,
    format_guides_csv,
    format_guides_fasta,
    predict_mmej,
)
from crisprguard.offtarget import OffTargetSummary, clean_dna, scan_offtargets
from crisprguard.samples import get_sample, list_samples
from crisprguard.web.server import serve


def _load_target_and_bg(
    seq: str | None,
    fasta_path: str | None,
    sample_id: str | None,
    bg_fasta_path: str | None = None,
    bg_seq: str | None = None,
    use_sample_bg: bool = True,
) -> tuple[str, str | None, str | None]:
    """Resolve target sequence, optional background FASTA, and optional target label."""
    target_text: str | None = None
    target_label: str | None = None
    background_text: str | None = None

    if sample_id:
        preset = get_sample(sample_id)
        target_text = preset.target_fasta
        target_label = preset.id
        if use_sample_bg:
            background_text = preset.background_fasta
    elif fasta_path:
        path = Path(fasta_path)
        if not path.is_file():
            raise ValueError(f"Target FASTA file not found: {fasta_path}")
        target_text = path.read_text(encoding="utf-8")
        target_label = path.stem
    elif seq:
        target_text = seq
        target_label = "custom_target"
    else:
        raise ValueError("Provide one of --seq, --fasta, or --sample.")

    if bg_fasta_path:
        bg_path = Path(bg_fasta_path)
        if not bg_path.is_file():
            raise ValueError(f"Background FASTA file not found: {bg_fasta_path}")
        background_text = bg_path.read_text(encoding="utf-8")
    elif bg_seq:
        background_text = bg_seq

    return target_text, background_text, target_label


def _write_output(text: str, output_path: str | None) -> None:
    if output_path:
        Path(output_path).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
        if not text.endswith("\n"):
            sys.stdout.write("\n")


def _format_design_table(report: DesignReport) -> str:
    lines: list[str] = [
        f"CrisprGuard v{__version__} | Target: {report.target_name} "
        f"({report.target_length} bp, GC {report.target_gc_percent}%%)".replace("%%", "%"),
        f"Enzyme: {report.enzyme_spec.name} | PAM: {report.enzyme_spec.pam_motif} | "
        f"Candidates: {report.total_candidates} (+{report.plus_strand_count} / -{report.minus_strand_count}) | "
        f"Background: {report.background_contigs} contigs ({report.background_bp} bp)",
        "-" * 116,
        f"{'Rank':<5} {'ID':<5} {'Str':<3} {'Cut':>5} {'Spacer':<24} {'PAM':<5} "
        f"{'GC%':>5} {'OnTgt':>6} {'Spec':>6} {'MMEJ%':>6} {'Comp':>6} {'OTs':>4}  Warnings",
        "-" * 116,
    ]
    if not report.guides:
        lines.append("No candidate guides matched the PAM and filter criteria.")
        return "\n".join(lines) + "\n"

    for g in report.guides:
        ot_count = g.offtarget_summary.total_hits if g.offtarget_summary else 0
        warns = ",".join(g.warnings) if g.warnings else "-"
        lines.append(
            f"{g.rank:<5} {g.id:<5} {g.strand:<3} {g.cut_pos:>5} {g.spacer:<24} {g.pam:<5} "
            f"{g.gc_percent:>5.1f} {g.on_target_score:>6.1f} {g.specificity_score:>6.1f} "
            f"{g.mmej.out_of_frame_probability:>6.1f} {g.composite_score:>6.1f} {ot_count:>4}  {warns}"
        )
    return "\n".join(lines) + "\n"


def _format_offtarget_table(summary: OffTargetSummary) -> str:
    lines: list[str] = [
        f"CrisprGuard CFD Off-Target Kinetics | Guide: {summary.guide} ({summary.enzyme})",
        f"Specificity Score: {summary.specificity_score:.2f} / 100 | "
        f"Total CFD Sum: {summary.total_cfd_sum:.4f} | "
        f"Hits (<= {summary.max_mismatches} mm): {summary.total_hits} "
        f"(High-Risk: {summary.high_risk_hits}) | Scanned: {summary.scanned_bp} bp",
        "-" * 110,
        f"{'#':<4} {'Contig':<24} {'Pos':>10} {'Str':<3} {'Alignment':<24} {'PAM':<5} "
        f"{'MM':>3} {'Seed':>4} {'Wob':>3} {'CFD':>8} {'Risk':<9}",
        "-" * 110,
    ]
    if not summary.top_hits:
        lines.append("No off-target sites found within the specified mismatch threshold.")
        return "\n".join(lines) + "\n"

    for idx, hit in enumerate(summary.top_hits, start=1):
        pos_str = f"{hit.start + 1}-{hit.end}"
        lines.append(
            f"{idx:<4} {hit.contig[:23]:<24} {pos_str:>10} {hit.strand:<3} "
            f"{hit.alignment:<24} {hit.pam:<5} {hit.mismatches:>3} {hit.seed_mismatches:>4} "
            f"{hit.wobble_mismatches:>3} {hit.cfd_score:>8.4f} {hit.risk_tier:<9}"
        )
    return "\n".join(lines) + "\n"


def _format_offtarget_csv(summary: OffTargetSummary) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(
        [
            "contig",
            "start_1based",
            "end_1based",
            "strand",
            "protospacer",
            "alignment",
            "pam",
            "pam_score",
            "mismatches",
            "seed_mismatches",
            "wobble_mismatches",
            "cfd_score",
            "risk_tier",
        ]
    )
    for hit in summary.top_hits:
        writer.writerow(
            [
                hit.contig,
                hit.start + 1,
                hit.end,
                hit.strand,
                hit.protospacer,
                hit.alignment,
                hit.pam,
                hit.pam_score,
                hit.mismatches,
                hit.seed_mismatches,
                hit.wobble_mismatches,
                hit.cfd_score,
                hit.risk_tier,
            ]
        )
    return buf.getvalue()


def _format_mmej_table(result: MMEJResult) -> str:
    ref_line = f"{result.left_flank}|{result.right_flank}"
    lines: list[str] = [
        f"CrisprGuard Bae MMEJ Predictor | Cut Pos: {result.cut_pos} | Window: +/-{result.window_bp} bp",
        f"Out-of-Frame Knockout Probability: {result.out_of_frame_probability:.1f}% ({result.knockout_grade}) | "
        f"In-Frame: {result.in_frame_probability:.1f}% | "
        f"Microhomology Score: {result.microhomology_score:.1f} ({result.total_patterns} patterns)",
        f"DSB Context: {ref_line}",
        "-" * 104,
        f"{'#':<4} {'MH':<9} {'Del':>5} {'Frame':<10} {'Score':>8}  Aligned Deletion Outcome",
        "-" * 104,
    ]
    if not result.top_patterns:
        lines.append("No 2-8 bp flanking microhomology repeats detected around cut site.")
        return "\n".join(lines) + "\n"

    for idx, pat in enumerate(result.top_patterns, start=1):
        frame = "Frameshift" if pat.is_frameshift else "In-Frame"
        lines.append(
            f"{idx:<4} {pat.microhomology:<9} -{pat.deletion_length:<4} {frame:<10} "
            f"{pat.pattern_score:>8.1f}  {pat.aligned_outcome}"
        )
    return "\n".join(lines) + "\n"


def _cmd_design(args: argparse.Namespace) -> int:
    if args.list_samples:
        for preset in list_samples():
            print(f"{preset.id:<16} {preset.name:<26} {preset.locus}")
        return 0

    target_text, bg_text, label = _load_target_and_bg(
        seq=args.seq,
        fasta_path=args.fasta,
        sample_id=args.sample,
        bg_fasta_path=args.background_fasta,
        bg_seq=args.background_seq,
        use_sample_bg=not args.no_sample_bg,
    )
    report = design_guides(
        sequence=target_text,
        enzyme=args.enzyme,
        background=bg_text,
        max_mismatches=args.max_mismatches,
        include_noncanonical_pam=not args.canonical_pam_only,
        mmej_window=args.mmej_window,
        top_n=args.top,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        target_name=label,
    )

    if args.format == "json":
        out = json.dumps(report.to_dict(), indent=2) + "\n"
    elif args.format == "csv":
        out = format_guides_csv(report.guides)
    elif args.format == "fasta":
        out = format_guides_fasta(report.guides, target_name=report.target_name)
    else:
        out = _format_design_table(report)

    _write_output(out, args.output)
    return 0


def _cmd_offtarget(args: argparse.Namespace) -> int:
    if args.background_fasta:
        path = Path(args.background_fasta)
        if not path.is_file():
            raise ValueError(f"Background FASTA file not found: {args.background_fasta}")
        bg_text = path.read_text(encoding="utf-8")
    elif args.background_seq:
        bg_text = args.background_seq
    elif args.sample:
        preset = get_sample(args.sample)
        bg_text = preset.background_fasta
    else:
        raise ValueError("Provide --background-fasta, --background-seq, or --sample.")

    summary = scan_offtargets(
        guide=args.guide,
        background=bg_text,
        enzyme=args.enzyme,
        max_mismatches=args.max_mismatches,
        include_noncanonical_pam=not args.canonical_pam_only,
        exclude_one_exact=args.exclude_one_exact,
        top_n=args.top,
    )

    if args.format == "json":
        out = json.dumps(summary.to_dict(), indent=2) + "\n"
    elif args.format == "csv":
        out = _format_offtarget_csv(summary)
    else:
        out = _format_offtarget_table(summary)

    _write_output(out, args.output)
    return 0


def _cmd_mmej(args: argparse.Namespace) -> int:
    target_text, _, _ = _load_target_and_bg(
        seq=args.seq,
        fasta_path=args.fasta,
        sample_id=args.sample,
        use_sample_bg=False,
    )
    cleaned = clean_dna(target_text)
    cut_pos = args.cut_pos if args.cut_pos is not None else len(cleaned) // 2
    result = predict_mmej(
        sequence=cleaned,
        cut_pos=cut_pos,
        window=args.window,
        top_n=args.top,
    )

    if args.format == "json":
        out = json.dumps(result.to_dict(), indent=2) + "\n"
    else:
        out = _format_mmej_table(result)

    _write_output(out, args.output)
    return 0


def _cmd_serve(args: argparse.Namespace) -> int:
    serve(host=args.host, port=args.port, open_browser=not args.no_open)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the CrisprGuard CLI argument parser (Python 3.14 percent-escape safe)."""
    parser = argparse.ArgumentParser(
        prog="crisprguard",
        description=(
            "CrisprGuard: Cas9/Cas12a CRISPR guide RNA designer, CFD off-target kinetics "
            "engine, and Bae MMEJ microhomology frameshift predictor."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- design ---
    p_design = subparsers.add_parser(
        "design",
        help="Design and rank SpCas9 or LbCas12a guide RNAs with MMEJ and CFD scores.",
    )
    p_design.add_argument("--seq", help="Raw target DNA sequence.")
    p_design.add_argument("--fasta", help="Path to target FASTA file.")
    p_design.add_argument(
        "--sample",
        choices=["tp53_exon5", "brca1_exon11", "pcsk9_exon1", "egfp_reporter"],
        help="Use a curated human/reporter exon preset.",
    )
    p_design.add_argument(
        "--list-samples",
        action="store_true",
        help="List all curated target gene presets and exit.",
    )
    p_design.add_argument(
        "--enzyme",
        choices=["spcas9", "cas12a"],
        default="spcas9",
        help="CRISPR endonuclease (spcas9 for 5'-20nt-NGG-3', cas12a for 5'-TTTV-23nt-3').",
    )
    p_design.add_argument(
        "--background-fasta",
        help="Optional reference FASTA file for CFD off-target scanning.",
    )
    p_design.add_argument(
        "--background-seq",
        help="Optional raw background DNA string for CFD off-target scanning.",
    )
    p_design.add_argument(
        "--no-sample-bg",
        action="store_true",
        help="Disable automatic loading of the preset's paired off-target background FASTA.",
    )
    p_design.add_argument(
        "--max-mismatches",
        type=int,
        default=4,
        help="Maximum mismatches allowed during CFD off-target scan (default: 4).",
    )
    p_design.add_argument(
        "--canonical-pam-only",
        action="store_true",
        help="Restrict off-target scanning to canonical PAMs only (NGG or TTTV).",
    )
    p_design.add_argument(
        "--mmej-window",
        type=int,
        default=25,
        help="Flanking bp window on each side of the DSB cut site for MMEJ analysis (default: 25).",
    )
    p_design.add_argument(
        "--min-gc",
        type=float,
        default=0.0,
        help="Minimum spacer GC percentage filter, 0 to 100 (default: 0.0).",
    )
    p_design.add_argument(
        "--max-gc",
        type=float,
        default=100.0,
        help="Maximum spacer GC percentage filter, 0 to 100 (default: 100.0).",
    )
    p_design.add_argument(
        "--top",
        type=int,
        default=15,
        help="Maximum number of top-ranked guides to output (default: 15).",
    )
    p_design.add_argument(
        "--format",
        choices=["table", "json", "csv", "fasta"],
        default="table",
        help="Output format (default: table).",
    )
    p_design.add_argument("--output", "-o", help="Write output to file path instead of stdout.")
    p_design.set_defaults(func=_cmd_design)

    # --- offtarget ---
    p_ot = subparsers.add_parser(
        "offtarget",
        help="Run PAM-proximal seed-weighted CFD off-target kinetics for a candidate guide.",
    )
    p_ot.add_argument(
        "--guide",
        required=True,
        help="Candidate spacer sequence (20 nt for SpCas9, 23 nt for LbCas12a).",
    )
    p_ot.add_argument("--background-fasta", help="Path to reference background FASTA file.")
    p_ot.add_argument("--background-seq", help="Raw reference background DNA string.")
    p_ot.add_argument(
        "--sample",
        choices=["tp53_exon5", "brca1_exon11", "pcsk9_exon1", "egfp_reporter"],
        help="Scan the paired mini-reference background of a curated sample.",
    )
    p_ot.add_argument(
        "--enzyme",
        choices=["spcas9", "cas12a"],
        default="spcas9",
        help="CRISPR endonuclease model (default: spcas9).",
    )
    p_ot.add_argument(
        "--max-mismatches",
        type=int,
        default=4,
        help="Maximum mismatch count to report (default: 4).",
    )
    p_ot.add_argument(
        "--canonical-pam-only",
        action="store_true",
        help="Only score canonical PAM off-target sites.",
    )
    p_ot.add_argument(
        "--exclude-one-exact",
        action="store_true",
        help="Exclude one exact canonical match (useful when target is inside background).",
    )
    p_ot.add_argument(
        "--top",
        type=int,
        default=20,
        help="Maximum top risk off-target loci to display (default: 20).",
    )
    p_ot.add_argument(
        "--format",
        choices=["table", "json", "csv"],
        default="table",
        help="Output format (default: table).",
    )
    p_ot.add_argument("--output", "-o", help="Write output to file path instead of stdout.")
    p_ot.set_defaults(func=_cmd_offtarget)

    # --- mmej ---
    p_mmej = subparsers.add_parser(
        "mmej",
        help="Predict Bae et al. MMEJ deletion patterns and out-of-frame knockout probability (%%).",
    )
    p_mmej.add_argument("--seq", help="Raw target DNA sequence.")
    p_mmej.add_argument("--fasta", help="Path to target FASTA file.")
    p_mmej.add_argument(
        "--sample",
        choices=["tp53_exon5", "brca1_exon11", "pcsk9_exon1", "egfp_reporter"],
        help="Use a curated sample sequence.",
    )
    p_mmej.add_argument(
        "--cut-pos",
        type=int,
        default=None,
        help="0-based DSB cut boundary index (defaults to sequence midpoint).",
    )
    p_mmej.add_argument(
        "--window",
        type=int,
        default=25,
        help="Flanking bp window on each side of the cut site (default: 25).",
    )
    p_mmej.add_argument(
        "--top",
        type=int,
        default=12,
        help="Number of top microhomology patterns to display (default: 12).",
    )
    p_mmej.add_argument(
        "--format",
        choices=["table", "json"],
        default="table",
        help="Output format (default: table).",
    )
    p_mmej.add_argument("--output", "-o", help="Write output to file path instead of stdout.")
    p_mmej.set_defaults(func=_cmd_mmej)

    # --- serve ---
    p_serve = subparsers.add_parser(
        "serve",
        help="Launch the interactive CrisprGuard dark-lab web workbench.",
    )
    p_serve.add_argument("--host", default="127.0.0.1", help="Bind host (default: 127.0.0.1).")
    p_serve.add_argument("--port", type=int, default=8765, help="Bind port (default: 8765).")
    p_serve.add_argument(
        "--no-open",
        action="store_true",
        help="Do not automatically open a web browser tab.",
    )
    p_serve.set_defaults(func=_cmd_serve)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI main entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except ValueError as exc:
        sys.stderr.write(f"Error: {exc}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
