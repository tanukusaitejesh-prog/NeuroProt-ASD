**Table 4. Zero-shot Spearman correlation with measured function (RNU4-2 saturation genome editing).**

| score                       |   n |    rho |     lo |    hi |
|:----------------------------|----:|-------:|-------:|------:|
| snrnavep_v2                 | 485 |  0.405 |  0.255 | 0.535 |
| nb_frac_prot                | 485 |  0.437 |  0.293 | 0.561 |
| frac_states_protein_contact | 485 |  0.347 |  0.198 | 0.482 |
| constraint_1_minus_oe       | 485 |  0.206 |  0.035 | 0.344 |
| pop_density_pos             | 485 |  0.149 | -0.006 | 0.287 |
| pop_density_w2              | 485 |  0.258 |  0.084 | 0.41  |
| cadd_phred                  | 411 |  0.072 | -0.07  | 0.208 |
| phylop447                   | 485 | -0.004 | -0.129 | 0.128 |

No model saw any functional measurement. This readout is free of the clinical-ascertainment circularity that affects population-derived scores, because the assay does not consult allele frequency.

| structure    | population            |   delta_rho |    lo |    hi |     p |
|:-------------|:----------------------|------------:|------:|------:|------:|
| snrnavep_v2  | constraint_1_minus_oe |       0.199 | 0.057 | 0.341 | 0.006 |
| nb_frac_prot | constraint_1_minus_oe |       0.231 | 0.093 | 0.363 | 0.002 |
| snrnavep_v2  | pop_density_pos       |       0.257 | 0.119 | 0.387 | 0     |
| nb_frac_prot | pop_density_pos       |       0.288 | 0.171 | 0.41  | 0     |
| snrnavep_v2  | pop_density_w2        |       0.147 | 0.012 | 0.281 | 0.028 |
| nb_frac_prot | pop_density_w2        |       0.178 | 0.06  | 0.3   | 0.002 |
| snrnavep_v2  | cadd_phred            |       0.252 | 0.053 | 0.439 | 0.01  |
| nb_frac_prot | cadd_phred            |       0.334 | 0.144 | 0.525 | 0     |
| snrnavep_v2  | phylop447             |       0.409 | 0.244 | 0.557 | 0     |
| nb_frac_prot | phylop447             |       0.44  | 0.271 | 0.604 | 0     |
