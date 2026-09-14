# Code to regenerate Figures from Tarnopol et al. (2026) "Horizontal gene transfer rivals gene duplication as a source of anti-parasitoid immune innovation in the Drosophilidae" 

Please see the respective folders for R code and input data to generate the figures in this manuscript. Full resolution versions of all main text figures are also included in these folders.

Figures 1 and 4 use alltips_new.csv as input data. The columns in this file are as follows:

The following data are included in alltips_new.csv: 
- node: Node # for each tip
- label: Species name corresponding with tip label encoded in the tree
- species: Reformatted species name
- PPO1: Copy # of PPO1 including pseudogenes
- PPO1_qc: Data quality for PPO1 annotations
    - intact: all copies appear functional
    - inc: incomplete data (i.e., from poor assembly quality or only detected via transcriptome reads)
    - pseud_part: ≥ 1 copy pseudogenized or segregating polymorphism for pseudogenized allele
    - pseud_all: all copies pseudogenized
- PPO2: Copy # of PPO2 including pseudogenes
- PPO2_qc: Data quality for PPO2 annotations
    - intact: all copies appear functional
    - inc: incomplete data (i.e., from poor assembly quality or only detected via transcriptome reads)
    - pseud_part: ≥ 1 copy pseudogenized or segregating polymorphism for pseudogenized allele
    - pseud_all: all copies pseudogenized
- PPO3: Copy # of PPO3 including pseudogenes
- PPO3_qc: Data quality for PPO3 annotations
    - intact: all copies appear functional
    - inc: incomplete data (i.e., from poor assembly quality or only detected via transcriptome reads)
    - pseud_part: ≥ 1 copy pseudogenized or segregating polymorphism for pseudogenized allele
    - pseud_all: all copies pseudogenized
- PPO4: Copy # of PPO4 including pseudogenes
- PPO4_qc: Data quality for PPO4 annotations
    - intact: all copies appear functional
    - inc: incomplete data (i.e., from poor assembly quality or only detected via transcriptome reads)
    - pseud_part: ≥ 1 copy pseudogenized or segregating polymorphism for pseudogenized allele
    - pseud_all: all copies pseudogenized
- ppo_pseudogene: Identity of pseudogenized PPO copies
- ppo_notes: Notes on manual PPO curation 
- cdtb_count: Copy # of cdtB including pseudogenes
- cdtb_exons: Number of exons in each copy of cdtB
- cdtb_qc: Data quality for cdtB annotations
    - intact: all copies appear functional
    - inc: incomplete data (i.e., from poor assembly quality or only detected via transcriptome reads)
    - pseud_part: ≥ 1 copy pseudogenized or segregating polymorphism for pseudogenized allele
    - pseud_all: all copies pseudogenized
- cdtb_presence: Binary code for presence/absence of cdtB to help with subsequent analyses:
    - 0 = cdtB absent
    - 1 = cdtB present
- cdtb_locus: Chromosomal location ID for each cdtB insertion. Each number indicates a different chromosomal location
- notes_defensetype: Notes on data collected for anti-parasitoid defense mechanism 
- mechanism: Anti-parasitoid defense mechanism elucidated from the literature
- source: First author and year of manuscript describing anti-parasitoid defense mechanism
- DOI: DOI URL for anti-parasitoid defense mechanism source
- code
- evidence_for_em:

