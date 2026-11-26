"""Comprehensive unit tests for CrisprGuard (design, CFD off-target, MMEJ, web API & CLI)."""

from __future__ import annotations

import contextlib
import io
import json
import unittest

from crisprguard.cli import build_parser, main as cli_main
from crisprguard.design import (
    design_guides,
    format_guides_csv,
    format_guides_fasta,
    predict_mmej,
    score_on_target,
)
from crisprguard.offtarget import (
    clean_dna,
    pam_activity_score,
    parse_fasta,
    reverse_complement,
    scan_offtargets,
    score_cfd_mismatch,
    score_offtarget_cfd,
)
from crisprguard.samples import get_sample, list_samples
from crisprguard.web.server import route_request


class TestSequenceAndPAMDesign(unittest.TestCase):
    """Verify exact SpCas9 (NGG) and LbCas12a (TTTV) PAM detection on + and - strands."""

    def test_spcas9_exact_plus_and_minus_strands(self) -> None:
        # + strand site: 20 nt spacer starting at index 3, PAM 'AGG' at [23, 26), cut at 20 (-3 bp upstream of PAM)
        plus_spacer = "ACGTACGTACGTACGTACGG"
        plus_pam = "AGG"
        seq_plus = "AAA" + plus_spacer + plus_pam + "TTT"
        report_plus = design_guides(seq_plus, enzyme="spcas9")

        plus_matches = [g for g in report_plus.guides if g.strand == "+" and g.spacer == plus_spacer]
        self.assertEqual(len(plus_matches), 1)
        g_plus = plus_matches[0]
        self.assertEqual(g_plus.pam, plus_pam)
        self.assertEqual(g_plus.start, 3)
        self.assertEqual(g_plus.end, 26)
        self.assertEqual(g_plus.pam_start, 23)
        self.assertEqual(g_plus.pam_end, 26)
        self.assertEqual(g_plus.cut_pos, 20)  # 23 - 3
        self.assertEqual(g_plus.cut_ts, 20)

        # - strand site: reverse complement of (minus_spacer + 'CGG') placed at [4, 27)
        minus_spacer = "GCTAGCTAGCTAGCTAGCTA"
        minus_pam = "CGG"
        rc_site = reverse_complement(minus_spacer + minus_pam)  # starts with 'CCG' on + strand
        seq_minus = "ATAT" + rc_site + "TATA"
        report_minus = design_guides(seq_minus, enzyme="spcas9")

        minus_matches = [g for g in report_minus.guides if g.strand == "-" and g.spacer == minus_spacer]
        self.assertEqual(len(minus_matches), 1)
        g_minus = minus_matches[0]
        self.assertEqual(g_minus.pam, minus_pam)
        self.assertEqual(g_minus.start, 4)
        self.assertEqual(g_minus.end, 27)
        self.assertEqual(g_minus.pam_start, 4)
        self.assertEqual(g_minus.pam_end, 7)
        # Cut is 3 bp into spacer from PAM boundary (index 7): 7 + 3 = 10
        self.assertEqual(g_minus.cut_pos, 10)

    def test_cas12a_exact_plus_and_minus_strands_and_tttt_rejection(self) -> None:
        # + strand site: PAM 'TTTC' at [2, 6), 23 nt spacer at [6, 29), NTS cut at 24 (+18), TS cut at 29 (+23)
        plus_pam = "TTTC"
        plus_spacer = "GCATGCTAGCTAGCATGCATGCA"  # 23 nt
        seq_plus = "CC" + plus_pam + plus_spacer + "GG"
        report_plus = design_guides(seq_plus, enzyme="cas12a")

        matches = [g for g in report_plus.guides if g.strand == "+" and g.spacer == plus_spacer]
        self.assertEqual(len(matches), 1)
        g = matches[0]
        self.assertEqual(g.pam, "TTTC")
        self.assertEqual(g.start, 2)
        self.assertEqual(g.end, 29)
        self.assertEqual(g.pam_start, 2)
        self.assertEqual(g.pam_end, 6)
        self.assertEqual(g.cut_pos, 24)  # 6 + 18
        self.assertEqual(g.cut_ts, 29)   # 6 + 23

        # - strand site: reverse complement of ('TTTA' + minus_spacer) placed at [3, 30)
        minus_pam = "TTTA"
        minus_spacer = "CGATCGATCGATCGATCGATCGA"  # 23 nt
        rc_site = reverse_complement(minus_pam + minus_spacer)
        seq_minus = "GGG" + rc_site + "CCC"
        report_minus = design_guides(seq_minus, enzyme="cas12a")

        m_matches = [g for g in report_minus.guides if g.strand == "-" and g.spacer == minus_spacer]
        self.assertEqual(len(m_matches), 1)
        gm = m_matches[0]
        self.assertEqual(gm.pam, "TTTA")
        self.assertEqual(gm.start, 3)
        self.assertEqual(gm.end, 30)
        self.assertEqual(gm.pam_start, 26)
        self.assertEqual(gm.pam_end, 30)
        # NTS cut 18 bp upstream of PAM boundary 26 -> 8; TS cut 23 bp upstream -> 3
        self.assertEqual(gm.cut_pos, 8)
        self.assertEqual(gm.cut_ts, 3)

        # Canonical Cas12a design must NOT treat 'TTTT' as a TTTV PAM
        # Exact 27 nt window: index 0 has PAM 'TTTT' + 23 nt spacer; index 1 ('TTTC') only has 22 nt remaining
        tttt_seq = "TTTT" + "CGCGCGCGCGCGCGCGCGCGCGC"
        report_tttt = design_guides(tttt_seq, enzyme="cas12a")
        self.assertEqual(report_tttt.total_candidates, 0)

    def test_pol3_tttt_terminator_penalized_and_flagged(self) -> None:
        clean_spacer = "GACTGCTACGACTGCTAGCG"  # 20 nt, 55% GC, G at pos 20, C at pos 18, no TTTT
        pol3_spacer = "GACTGTTTTGACTGCTAGCG"   # Identical except internal TTTT motif

        clean_score = score_on_target(clean_spacer, "CGG", enzyme="spcas9")
        pol3_score = score_on_target(pol3_spacer, "CGG", enzyme="spcas9")

        self.assertFalse(clean_score.pol3_terminator)
        self.assertTrue(pol3_score.pol3_terminator)
        self.assertIn("pol3_terminator_tttt", pol3_score.warnings)
        self.assertGreater(clean_score.score - pol3_score.score, 25.0)

        # Same check for LbCas12a
        cas12_clean = score_on_target("CGACGCTAGCTAGCATCGATCGA", "TTTC", enzyme="cas12a")
        cas12_pol3 = score_on_target("CGACGCTTTTTAGCATCGATCGA", "TTTC", enzyme="cas12a")
        self.assertTrue(cas12_pol3.pol3_terminator)
        self.assertGreater(cas12_clean.score - cas12_pol3.score, 25.0)

    def test_position_preferences_and_hairpin_penalty(self) -> None:
        # SpCas9 pos 20: G strongly preferred over C
        g20 = score_on_target("GCTAGCTAGCTAGCTAGCCG", "CGG", enzyme="spcas9")
        c20 = score_on_target("GCTAGCTAGCTAGCTAGCCC", "CGG", enzyme="spcas9")
        self.assertGreater(g20.score, c20.score)

        # Strong self-folding hairpin (6 bp GC-rich stem + 4 nt loop) penalized
        hairpin_spacer = "GCGCGCAAAAGCGCGCTACG"
        normal_spacer = "GACTGCTAGCAGTCATGCCG"
        hp = score_on_target(hairpin_spacer, "AGG", enzyme="spcas9")
        norm = score_on_target(normal_spacer, "AGG", enzyme="spcas9")
        self.assertGreaterEqual(hp.hairpin_stem_bp, 6)
        self.assertIn("strong_hairpin", hp.warnings)
        self.assertLess(hp.score, norm.score)


class TestCFDOffTargetKinetics(unittest.TestCase):
    """Verify PAM-proximal seed mismatch kinetics, rG:dT wobble tolerance & MIT/CFD scoring."""

    def test_seed_mismatch_penalized_more_than_distal_spcas9_and_cas12a(self) -> None:
        guide_cas9 = "GCTAGCTAGCTAGCTAGCTG"
        # Same mismatch type (guide G vs protospacer C -> rG:dG transversion) at pos 1 (distal) vs pos 20 (core seed)
        ot_distal = "CCTAGCTAGCTAGCTAGCTG"
        ot_seed = "GCTAGCTAGCTAGCTAGCTC"

        cfd_distal, _, align_d, det_d = score_offtarget_cfd(guide_cas9, ot_distal, "AGG", enzyme="spcas9")
        cfd_seed, _, align_s, det_s = score_offtarget_cfd(guide_cas9, ot_seed, "AGG", enzyme="spcas9")

        self.assertLess(cfd_seed, cfd_distal)
        self.assertTrue(det_s[0].is_seed)
        self.assertFalse(det_d[0].is_seed)
        self.assertEqual(align_d, "C...................")
        self.assertEqual(align_s, "...................C")

        # LbCas12a has opposite polarity: 5' seed (pos 1..8) vs 3' distal (pos 18..23)
        guide_c12 = "GCTAGCTAGCTAGCTAGCTAGCG"
        c12_seed = "CCTAGCTAGCTAGCTAGCTAGCG"    # mismatch at pos 1 (5' PAM-proximal seed)
        c12_distal = "GCTAGCTAGCTAGCTAGCTAGCC"  # mismatch at pos 23 (3' PAM-distal)
        cfd_c12_seed, _, _, det_c12_s = score_offtarget_cfd(guide_c12, c12_seed, "TTTC", enzyme="cas12a")
        cfd_c12_dist, _, _, det_c12_d = score_offtarget_cfd(guide_c12, c12_distal, "TTTC", enzyme="cas12a")
        self.assertLess(cfd_c12_seed, cfd_c12_dist)
        self.assertTrue(det_c12_s[0].is_seed)
        self.assertFalse(det_c12_d[0].is_seed)

    def test_rg_dt_wobble_scores_higher_than_transversion(self) -> None:
        # At position 16 (seed): guide has 'G' (rG).
        # Protospacer 'A' -> target DNA 'T' -> rG:dT wobble.
        # Protospacer 'C' -> target DNA 'G' -> rG:dG transversion.
        # Protospacer 'T' -> target DNA 'A' -> rG:dA transversion.
        wobble = score_cfd_mismatch(16, "G", "A", enzyme="spcas9")
        transv_gg = score_cfd_mismatch(16, "G", "C", enzyme="spcas9")
        transv_ga = score_cfd_mismatch(16, "G", "T", enzyme="spcas9")

        self.assertEqual(wobble.pair, "rG:dT")
        self.assertTrue(wobble.is_wobble)
        self.assertFalse(transv_gg.is_wobble)
        self.assertGreater(wobble.cfd_factor, transv_gg.cfd_factor)
        self.assertGreater(wobble.cfd_factor, transv_ga.cfd_factor)

    def test_pam_variant_weights_and_aggregate_specificity(self) -> None:
        self.assertEqual(pam_activity_score("AGG", "spcas9"), 1.0)
        self.assertAlmostEqual(pam_activity_score("CAG", "spcas9"), 0.259, places=3)
        self.assertAlmostEqual(pam_activity_score("TGA", "spcas9"), 0.069, places=3)
        self.assertEqual(pam_activity_score("CAG", "spcas9", include_noncanonical=False), 0.0)

        guide = "GCTGATACCGTTACGCAGTG"
        # Embed two known off-targets: one + strand NAG (0 mm), one - strand NGG (1 distal mm)
        ot1_plus = guide + "CAG"  # CFD = 0.259
        ot2_protospacer = "CCTGATACCGTTACGCAGTG"  # pos 1 mismatch
        ot2_minus = reverse_complement(ot2_protospacer + "AGG")
        bg = f">decoy1\nAAAAA{ot1_plus}CCCCC\n>decoy2\nTTTTT{ot2_minus}GGGGG\n"

        summary = scan_offtargets(guide, bg, enzyme="spcas9", max_mismatches=4)
        self.assertEqual(summary.total_hits, 2)
        strands = {h.strand for h in summary.top_hits}
        self.assertEqual(strands, {"+", "-"})

        expected_spec = round(100.0 / (1.0 + sum(h.cfd_score for h in summary.top_hits)), 2)
        self.assertAlmostEqual(summary.specificity_score, expected_spec, places=2)


class TestMMEJFrameshiftPredictor(unittest.TestCase):
    """Verify Bae et al. MMEJ microhomology detection and in-frame vs frameshift scoring."""

    def test_distinguishes_in_frame_vs_frameshift_microhomology(self) -> None:
        # Deliberate 6-bp GC-rich microhomology 'GCCGGC' placed cleanly across the cut site:
        # Frameshift (del_len = 7, 7 % 3 == 1): left ends with 'A' + GCCGGC, right starts with 'T' + GCCGGC
        left_common = "TATACGTACGTATACGTA" + "GCCGGC"  # length 24, cut_pos = 24
        right_fs = "T" + "GCCGGC" + "ATCTAGCATGCTATCAT"
        seq_fs = left_common + right_fs
        res_fs = predict_mmej(seq_fs, cut_pos=24, window=25)

        # In-frame (del_len = 6, 6 % 3 == 0): right starts immediately with GCCGGC -> del_len = 6
        right_if = "GCCGGC" + "TACTAGCATGCTATCAT"
        seq_if = left_common + right_if
        res_if = predict_mmej(seq_if, cut_pos=24, window=25)

        # Verify top pattern detects 'GCCGGC' with exact deletion lengths
        self.assertEqual(res_fs.top_patterns[0].microhomology, "GCCGGC")
        self.assertEqual(res_fs.top_patterns[0].deletion_length, 7)
        self.assertTrue(res_fs.top_patterns[0].is_frameshift)

        self.assertEqual(res_if.top_patterns[0].microhomology, "GCCGGC")
        self.assertEqual(res_if.top_patterns[0].deletion_length, 6)
        self.assertFalse(res_if.top_patterns[0].is_frameshift)

        # Frameshift sequence must yield significantly higher Out-of-Frame Knockout Probability %
        self.assertGreater(res_fs.out_of_frame_probability, 65.0)
        self.assertLess(res_if.out_of_frame_probability, 45.0)
        self.assertGreater(res_fs.out_of_frame_probability, res_if.out_of_frame_probability + 25.0)

        # Verify aligned deletion outcome string preserves total flanking length
        window_len = len(res_fs.left_flank) + len(res_fs.right_flank)
        for pat in res_fs.top_patterns:
            self.assertEqual(len(pat.aligned_outcome), window_len)
            self.assertEqual(pat.aligned_outcome.count("-"), pat.deletion_length)


class TestCuratedSamplesAndExporters(unittest.TestCase):
    """Verify all curated presets work across both enzymes and export cleanly."""

    def test_all_presets_design_both_enzymes(self) -> None:
        samples = list_samples()
        self.assertEqual(len(samples), 4)
        for preset in samples:
            r_cas9 = design_guides(preset.target_fasta, enzyme="spcas9", background=preset.background_fasta)
            r_cas12 = design_guides(preset.target_fasta, enzyme="cas12a", background=preset.background_fasta)
            self.assertGreater(r_cas9.total_candidates, 10, f"{preset.id} SpCas9 guides")
            self.assertGreater(r_cas9.plus_strand_count, 0)
            self.assertGreater(r_cas9.minus_strand_count, 0)
            self.assertGreater(r_cas12.total_candidates, 2, f"{preset.id} Cas12a guides")

            csv_text = format_guides_csv(r_cas9.guides[:3])
            self.assertIn("rank,id,enzyme,strand", csv_text)
            fasta_text = format_guides_fasta(r_cas9.guides[:3], target_name=preset.id)
            self.assertIn(f">sg01_{preset.id}", fasta_text)

    def test_fasta_and_dna_validation(self) -> None:
        records = parse_fasta(">seq1 desc\naugc\n>seq2\nACGT\n")
        self.assertEqual(records, [("seq1", "ATGC"), ("seq2", "ACGT")])
        self.assertEqual(clean_dna("5'-acgu 123-3'".replace("5'-", "").replace("-3'", "")), "ACGT")
        with self.assertRaises(ValueError):
            clean_dna("ACGTXYZ")


class TestWebAPIRouter(unittest.TestCase):
    """Verify web server API routes return 200/400/403/404 cleanly without binding sockets."""

    def test_static_and_health_and_samples(self) -> None:
        status, body, ctype = route_request("GET", "/api/health")
        self.assertEqual(status, 200)
        self.assertIn("application/json", ctype)
        payload = json.loads(body.decode("utf-8"))
        self.assertEqual(payload["status"], "ok")

        status, body, _ = route_request("GET", "/api/samples")
        self.assertEqual(status, 200)
        self.assertEqual(len(json.loads(body.decode("utf-8"))["samples"]), 4)

        for asset, expected_mime in [("/", "text/html"), ("/styles.css", "text/css"), ("/app.js", "javascript")]:
            st, data, ct = route_request("GET", asset)
            self.assertEqual(st, 200)
            self.assertIn(expected_mime, ct)
            self.assertGreater(len(data), 100)

        # Path traversal protection
        st, _, _ = route_request("GET", "/../../pyproject.toml")
        self.assertIn(st, (403, 404))

        # Unknown API route
        st, _, _ = route_request("GET", "/api/does-not-exist")
        self.assertEqual(st, 404)

    def test_api_design_offtarget_mmej_export_200_and_400(self) -> None:
        preset = get_sample("tp53_exon5")

        # Valid POST /api/design
        req = json.dumps({
            "sequence": preset.target_fasta,
            "background": preset.background_fasta,
            "enzyme": "spcas9",
            "max_mismatches": 4,
            "top_n": 5,
        }).encode("utf-8")
        st, body, _ = route_request("POST", "/api/design", req)
        self.assertEqual(st, 200)
        res = json.loads(body.decode("utf-8"))
        self.assertEqual(len(res["guides"]), 5)
        self.assertIn("mmej", res["guides"][0])
        self.assertIn("offtarget_summary", res["guides"][0])

        # Valid POST /api/offtarget
        guide_seq = res["guides"][0]["spacer"]
        ot_req = json.dumps({
            "guide": guide_seq,
            "background": preset.background_fasta,
            "enzyme": "spcas9",
        }).encode("utf-8")
        st, body, _ = route_request("POST", "/api/offtarget", ot_req)
        self.assertEqual(st, 200)
        self.assertIn("specificity_score", json.loads(body.decode("utf-8")))

        # Valid POST /api/mmej
        mmej_req = json.dumps({"sequence": preset.target_fasta, "cut_pos": 100}).encode("utf-8")
        st, body, _ = route_request("POST", "/api/mmej", mmej_req)
        self.assertEqual(st, 200)
        self.assertIn("out_of_frame_probability", json.loads(body.decode("utf-8")))

        # Valid POST /api/export (CSV & FASTA)
        exp_req = json.dumps({"sequence": preset.target_fasta, "format": "csv", "top_n": 3}).encode("utf-8")
        st, body, ct = route_request("POST", "/api/export", exp_req)
        self.assertEqual(st, 200)
        self.assertIn("text/csv", ct)
        self.assertIn("on_target_score", body.decode("utf-8"))

        # Clean 400 Bad Request cases
        bad_payloads = [
            ("/api/design", b"{not-json"),
            ("/api/design", json.dumps({"sequence": ""}).encode("utf-8")),
            ("/api/design", json.dumps({"sequence": "ACGT123XYZ"}).encode("utf-8")),
            ("/api/design", json.dumps({"sequence": "ACGT" * 10, "enzyme": "cas99"}).encode("utf-8")),
            ("/api/offtarget", json.dumps({"guide": "ACGT", "background": "ACGT" * 10}).encode("utf-8")),
            ("/api/mmej", json.dumps({"sequence": "ACGT" * 10}).encode("utf-8")),
            ("/api/export", json.dumps({"sequence": "ACGT" * 10, "format": "xml"}).encode("utf-8")),
        ]
        for endpoint, payload in bad_payloads:
            st, body, _ = route_request("POST", endpoint, payload)
            self.assertEqual(st, 400, f"Expected 400 for {endpoint} with {payload!r}")
            self.assertIn("error", json.loads(body.decode("utf-8")))


class TestCLISubcommands(unittest.TestCase):
    """Verify CLI subcommands & argparse help strings work cleanly (including Python 3.14)."""

    def _run_cli(self, argv: list[str]) -> tuple[int, str, str]:
        out_buf = io.StringIO()
        err_buf = io.StringIO()
        with contextlib.redirect_stdout(out_buf), contextlib.redirect_stderr(err_buf):
            code = cli_main(argv)
        return code, out_buf.getvalue(), err_buf.getvalue()

    def test_help_strings_percent_safe(self) -> None:
        parser = build_parser()
        # Formatting help across root and all subcommands verifies Python 3.14 %-escaping
        root_help = parser.format_help()
        self.assertIn("crisprguard", root_help.lower())
        for sub_name in ("design", "offtarget", "mmej", "serve"):
            with self.assertRaises(SystemExit) as ctx:
                with contextlib.redirect_stdout(io.StringIO()):
                    parser.parse_args([sub_name, "--help"])
            self.assertEqual(ctx.exception.code, 0)

    def test_cli_design_offtarget_mmej_exit_zero(self) -> None:
        code, out, _ = self._run_cli(["design", "--list-samples"])
        self.assertEqual(code, 0)
        self.assertIn("tp53_exon5", out)

        code, out, _ = self._run_cli(["design", "--sample", "tp53_exon5", "--enzyme", "spcas9", "--top", "5"])
        self.assertEqual(code, 0)
        self.assertIn("CrisprGuard", out)
        self.assertIn("sg01", out)

        code, out, _ = self._run_cli(
            ["design", "--sample", "brca1_exon11", "--enzyme", "cas12a", "--format", "json", "--top", "3"]
        )
        self.assertEqual(code, 0)
        parsed = json.loads(out)
        self.assertEqual(parsed["enzyme"], "cas12a")
        self.assertGreaterEqual(len(parsed["guides"]), 1)

        code, out, _ = self._run_cli(
            ["offtarget", "--guide", "GACCTGCCCTGTGCAGCTGT", "--sample", "tp53_exon5", "--format", "table"]
        )
        self.assertEqual(code, 0)
        self.assertIn("Specificity Score", out)

        code, out, _ = self._run_cli(
            ["mmej", "--sample", "egfp_reporter", "--cut-pos", "120", "--format", "json"]
        )
        self.assertEqual(code, 0)
        mmej_data = json.loads(out)
        self.assertIn("out_of_frame_probability", mmej_data)


if __name__ == "__main__":
    unittest.main()
