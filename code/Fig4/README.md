# Code to generate Figure 4 from Tarnopol et al. (2026) "Horizontal gene transfer rivals gene duplication as a source of anti-parasitoid immune innovation in the Drosophilidae"

The pgls_input.csv file contains the following data: 

- label: Tip label retained in the pared tree
- bd.ms_0: Net diversification rate calculated under no extinction (e = 0)
- bd.ms_0.5: Net diversification rate calculated under moderate extinction (e = 0.5)
- bd.ms_0.9: Net diversification rate calculated under high extinction (e = 0.9)
- cdtb_presence: Character encoding of cdtB presence/absence
    - y = cdtB present
    - n = cdtB absent
- cdtb_prop: Proportion of species in representative clade that are cdtB positive
- cdtb_species: Binary encoding of cdtB presence/absence
    - 0 = cdtB absent
    - 1 = cdtB present
