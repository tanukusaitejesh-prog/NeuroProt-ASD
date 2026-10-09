**Table 2. Ablation of every component and of each conformational state.**

| Ablation                                             | Dominant   | Recessive   | dDominant                                                                                      |
|:-----------------------------------------------------|:-----------|:------------|:-----------------------------------------------------------------------------------------------|
| **None (full MuSE, 11 structures)**                  | **0.806**  | **0.636**   | -                                                                                              |
| Remove the contact-graph smoothing term              | 0.788      | 0.635       | -0.018                                                                                         |
| Single structure, chosen without the test gene       | 0.597      | 0.654       | -0.209                                                                                         |
| Single structure, chosen WITH the test gene (oracle) | 0.795      | 0.675       | -0.011                                                                                         |
| Only 5Z56 (Bact-mature)                              | 0.795      | 0.644       | -0.011                                                                                         |
| Only 6FF7 (Bact)                                     | 0.792      | 0.634       | -0.015                                                                                         |
| Only 5XJC (C*)                                       | 0.752      | 0.608       | -0.054                                                                                         |
| Only 7DVQ (minor-Bact)                               | 0.747      | 0.675       | -0.060                                                                                         |
| Only 6QDV (P)                                        | 0.687      | 0.599       | -0.120                                                                                         |
| Only 6QW6 (tri-snRNP)                                | 0.655      | 0.637       | -0.152                                                                                         |
| Only 6QX9 (pre-B)                                    | 0.608      | 0.604       | -0.198                                                                                         |
| Only 6Y5Q (U2-snRNP)                                 | 0.587      | 0.625       | -0.219                                                                                         |
| Only 8Y6O (minor-preB)                               | 0.582      | 0.654       | -0.224                                                                                         |
| Only 3JCR (tri-snRNP)                                | 0.581      | 0.629       | -0.226                                                                                         |
| Only 5O9Z (B)                                        | 0.551      | 0.607       | -0.255                                                                                         |
| ADD conformational-trajectory features (scripts/51)  | n.s.       | n.s.        | +0.030, but fails its degree-matched null (p = 1.00) and the gain is one gene with 7 positives |
| ADD protein-partner identity (scripts/52)            | n.s.       | n.s.        | -0.026 (does not transfer across genes)                                                        |

Mean AUROC across held-out units. The final two rows are components that were specified with a kill criterion in advance, tested, and rejected; they are reported because a negative result about an obvious extension is evidence about the problem, not an omission.
