## Data
Biological Source: Whole-genome sequencing (WGS) of the HCC1395 triple-negative breast cancer cell line paired with its matched normal lymphoblastoid cell line (HCC1395BL) from the exact same donor.

truth_somatic_snv.vcf.gz: Contains 39,560 experimentally validated high-confidence somatic single-nucleotide variants (SNVs) confirmed across multiple sequencing platforms and callers. Any variant in this file gets label 1 (Somatic).

hc_regions.bed: Defines the ~2.5 billion base pairs of the human genome where sequencing coverage and mappability are high enough to reliably call variants. Variants outside these regions (e.g., centromeres/telomeres) have uncertain truth labels, so we filter to hc_regions.bed first.

raw_mutect2.vcf.gz: The raw, unfiltered output of MuTect2 on one Illumina WGS replicate. It contains tens of thousands of candidate calls—including true somatic mutations, germline leakage, and sequencing/mapping artifacts. Any candidate not in the truth VCF gets label 0 (Non-Somatic).
