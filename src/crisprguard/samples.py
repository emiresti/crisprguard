"""Curated human & reporter CRISPR target exon sequences with mini reference backgrounds."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class SamplePreset:
    """Curated target gene exon + realistic homologous off-target background FASTA."""

    id: str
    name: str
    gene: str
    locus: str
    description_en: str
    description_tr: str
    target_fasta: str
    background_fasta: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


SAMPLES: dict[str, SamplePreset] = {
    "tp53_exon5": SamplePreset(
        id="tp53_exon5",
        name="Human TP53 Exon 5",
        gene="TP53",
        locus="chr17:7675053-7675324 (DNA-Binding Domain)",
        description_en=(
            "Tumor suppressor p53 Exon 5 spanning the R175H structural hotspot, paired with "
            "homologous TP63/TP73 paralog decoy loci for realistic SpCas9 & Cas12a CFD screening."
        ),
        description_tr=(
            "R175H yapisal sicak noktasini iceren tumor baskilayici p53 Ekzon 5 dizisi; gercekci "
            "SpCas9 ve Cas12a CFD hedef-disi taramasi icin TP63/TP73 paralog bolgeleriyle eslestirilmistir."
        ),
        target_fasta=(
            ">TP53_Exon5_NM_000546\n"
            "TTCTTTGCTGCCGTCTTCCAGTTGCTTTATCTGTTCACTTGTGCCCTGACTTTCAACTCT\n"
            "GTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCC\n"
            "AAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGC\n"
            "GCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCAC\n"
            "CATGAGCGCTGCTCAGATAGCGATGGTGAGCAGCTGGGGCTGGAGAGACGACAGGGCTGG\n"
        ),
        background_fasta=(
            ">TP63_DBD_homolog_chr3\n"
            "AAGTGCCTTGACCACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTCGACTCCACACCAC\n"
            "AGTCGCTGCCCAGTACACATGACGGAGGTTGTGAGACGCTGCCCCCAGCATGAATGCTGC\n"
            "TGCTCCTTTCTTTGCTGCCGTCTTCCAATTGCTTTATCTGTTAACTTGCACCCTGAATTT\n"
            ">TP73_paralog_chr1\n"
            "CCTGGTGAGGTTGTGAGGCGCTACCCCCACCATGAACGCTGCTCAGACAGCGATGGTGAG\n"
            "CTCGACGCCATGGCCATCTACAAGAAGTGACAGCACATGACGAAGGTTGTGAGGCACTGC\n"
            "GTTCAACTCTGTCTCCTTCCTCTTCCTACAGTCCTCCCCTGCCCTCAACGAGATGTTTTG\n"
            ">chr17_intergenic_decoy\n"
            "GCTAGCTGCCCTCAACAAGATGTTCTGCCAACTGGCCAGGACCTGCCCTATGCAGCTATG\n"
            "TGGGTTGATTCCACACTCCCGCCCGGCACCTGCGTCCGCACCATGGCCATCTACAAACAG\n"
            "CCAGTTGCTTTATCTGTTCACTTGTGCCCTGACTTTCATCTCTGTCTCCTTCCTCTTCAG\n"
        ),
    ),
    "brca1_exon11": SamplePreset(
        id="brca1_exon11",
        name="Human BRCA1 Exon 11",
        gene="BRCA1",
        locus="chr17q21.31 (RAD51-Interaction Domain)",
        description_en=(
            "AT-balanced DNA repair gene BRCA1 Exon 11 segment rich in both LbCas12a (TTTV) "
            "and SpCas9 (NGG) target sites, paired with BRCA1P1 pseudogene off-target loci."
        ),
        description_tr=(
            "Hem LbCas12a (TTTV) hem SpCas9 (NGG) bolgeleri bakimindan zengin DNA onarim geni "
            "BRCA1 Ekzon 11 segmenti; BRCA1P1 psodogen hedef-disi lokuslari ile eslestirilmistir."
        ),
        target_fasta=(
            ">BRCA1_Exon11_NM_007294\n"
            "GCTGCTTGTGAATTTTCTGAGACGGATGTAACAAATACTGAACATCATCAACCCAGTAAT\n"
            "AATGATTTGAACACCACTGAGAAGCGTGCAGCTGAGAGGCATCCAGAAAAGTATCAGGGT\n"
            "AGTTCTGTTTCAAACTTGCATGTGGAGCCATGTGGCACAAATACTCATGCCAGCTCATTA\n"
            "CAGCATGAGAACAGCAGTTTATTACTCACTAAAGACAGAATGAATGTAGAAAAGGCTGAA\n"
            "TTCTGTAATAAAAGCAAACAGCCTGGCTTAGCAAGGAGCCAACATAACAGATGGCTGCTG\n"
        ),
        background_fasta=(
            ">BRCA1P1_pseudogene_chr17\n"
            "CCAGTGCTTGTGAATTTCCTGAGACAGATGTAACAAACACTGAACATCGTCAACCCAGTA\n"
            "AATGATTTGAACACCGCTGAGAAGCGTGCAGCTGAAAAGCATCCAGAAAAGTATCAAGGG\n"
            "TTCTGTTTCAAACTTGCATGTGGAGCCATATGGCACAAAGACTCATGCCAGCTCACTACA\n"
            ">BARD1_adjacent_locus_chr2\n"
            "AGTTCAGCAGTTTATTACTCACTAAAGACAGAATGAACGTAGAAAAGGCTGAATTCTGTG\n"
            "ATGATAAAAGCAAACAGCCTGGCTTAGCGAGGAGCCAACATGACAGATGGCTGCTGCCAA\n"
            "CTGAAGCGTGCAGCTGAGAGGCAACCAGAAAAGTATCGGGGTAGTTCTGTTTCAAGCTTG\n"
            ">chr17q21_intergenic_decoy\n"
            "GGAATTTTCTGAGACGGATGTAACAAATATTGAGCATCATCAACCCAGCAATAATGATTT\n"
            "GAACACCACTGAAAAGCGTGCAGCTGGGAGGCATCCAGAGAAGTATCAGGGTAGTTCTAG\n"
        ),
    ),
    "pcsk9_exon1": SamplePreset(
        id="pcsk9_exon1",
        name="Human PCSK9 Exon 1",
        gene="PCSK9",
        locus="chr1p32.3 (Therapeutic LDL-C Knockout)",
        description_en=(
            "Clinical gene-editing target PCSK9 Exon 1 (signal peptide & prodomain) used in "
            "cardiovascular LDL-cholesterol lowering studies, with chr1/chr6 off-target decoys."
        ),
        description_tr=(
            "Kardiyovaskuler LDL-kolesterol dusurucu klinik gen duzenleme calismalarinda kullanilan "
            "PCSK9 Ekzon 1 (sinyal peptidi ve prodomain) hedefi ve chr1/chr6 hedef-disi lokuslari."
        ),
        target_fasta=(
            ">PCSK9_Exon1_NM_174936\n"
            "CTTTCGCTTTCAGCAAGGAGATGGGCACCGTCAGCTCCAGGCGGTCCTGGTGGCCGCTGC\n"
            "CACTGCTGCTGCTGCTGCTGCTGCTCCTGGGTCCCGCGGGCGCCCGTGCGCAGGAGGACG\n"
            "AGGACGGCGACTACGAGGAGCTGGTGCTAGCCTTGCGTTCCGAGGAGGACGGCCTGGCCG\n"
            "AAGCACCCGAGCACGGAACCACAGCCACCTTCCACCGCTGCGCCAAGGATCCGTGGAGGT\n"
            "TGCCTGGCACCTACGTGGTGGTGCTGAAGGAGGAGACCCACCTCTCGGTAAGTCAGGTTT\n"
            "GAAACC\n"
        ),
        background_fasta=(
            ">chr1p32_paralog_region\n"
            "GGAGCTTTCAGCAAGGAGATAGGCACCGTCAGCTCCAGGAGGTCCTGGTGGCCACTGCCA\n"
            "GGACGAGGACGGCGACTACGAGGAGCTGGTACTAGCCTTGCGTTCCGAGGAGGATGGCCT\n"
            "CCCGAGCACGGAACCACAGCCACCTTCCATCGCTGCGCCAAGGATCCGTAGAGGTTGCCT\n"
            ">LDLR_distal_decoy_chr19\n"
            "TCACACCTACGTGGTGGTGCTGAAGGAGGAGACACACCTCTCGGTAAGTCAGGTTTGAAA\n"
            "CAGGCCTGGCCGAAGCACCCGAGCACAGAACCACAGCCACCTTCCACCACTGCGCCAAGA\n"
            "CCACTGCTGCTGCTGCTACTGCTGCTCCTGGGTCCCGTGGGCGCCCGTGCGCAGGAGGAC\n"
            ">chr6_intergenic_offtarget\n"
            "AAGGAGATGGGCACCGTCAGCTCAAGGCGGTCCTGGTAGCCGCTGCCACTGCTGCTACTG\n"
            "TGCGTTCCGAGGAGGACGGCCTGGCCAAAGCACCCGAGCACGGAACCACAACCACCTTCC\n"
        ),
    ),
    "egfp_reporter": SamplePreset(
        id="egfp_reporter",
        name="eGFP Reporter Knockout",
        gene="eGFP",
        locus="Chromophore Domain (T65-Y66-G67)",
        description_en=(
            "Classic wet-lab CRISPR knockout validation target spanning the eGFP fluorophore "
            "(T65-Y66-G67) with EYFP/ECFP homologous variant decoys to test mismatch specificity."
        ),
        description_tr=(
            "Islak laboratuvar CRISPR nakavt dogrulamasi icin eGFP florofor bolgesini (T65-Y66-G67) "
            "kapsayan klasik hedef; EYFP/ECFP homolog varyantlariyla ozgulluk testi sunar."
        ),
        target_fasta=(
            ">eGFP_Chromophore_CDS\n"
            "GTTTAGTGAACCGTCAGAATTCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGG\n"
            "GGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTC\n"
            "CGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCAC\n"
            "CGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTG\n"
            "CTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGA\n"
            "AGGCTACGTCCAGTTTGAAAC\n"
        ),
        background_fasta=(
            ">EYFP_variant_locus\n"
            "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTTGAGCTGGAC\n"
            "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCAGCGAGGGCGAGGGCGATGCCACCTAC\n"
            "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCATGGCCCACC\n"
            "CTCGTGACCACCTTCGGCTACGGCCTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAA\n"
            ">ECFP_homolog_locus\n"
            "TTTAGTGAACCGTCAGAATTCGCCATCATGGTGAGCAAGGGCGAGGAGCTGTTCACTGGG\n"
            "CCACAAGTTCAGCGTGTCTGGCGAGGGCGAGGGCGATGCCACCTACAGCAAGCTGACCCT\n"
            "ACCACCCTGACCTGGGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCAT\n"
            ">AAVS1_human_background\n"
            "CCCCTCCACCCCACAGTGGGGCCACTAGGGACAGGATTGGTGACAGAAAAGCCCCATCCT\n"
            "TAGGCCTCCTCCTTCCTAGTCTCCTGATATTGGGTCTAACCCCCACCTCCTGTTAGGCAG\n"
        ),
    ),
}


def list_samples() -> list[SamplePreset]:
    """Return all curated sample presets in display order."""
    return list(SAMPLES.values())


def get_sample(sample_id: str) -> SamplePreset:
    """Retrieve a sample preset by ID (case-insensitive)."""
    key = sample_id.strip().lower()
    if key not in SAMPLES:
        available = ", ".join(SAMPLES.keys())
        raise ValueError(f"Unknown sample '{sample_id}'. Available samples: {available}")
    return SAMPLES[key]
