**Table 1. MuSE against 13 baselines under the pre-registered held-out-unit protocol.**

| Baseline                                      | Dominant AUROC   | Dominant dAUROC [95% CI]   | Recessive AUROC   | Recessive dAUROC [95% CI]   |
|:----------------------------------------------|:-----------------|:---------------------------|:------------------|:----------------------------|
| **MuSE (this work)**                          | **0.810**        | -                          | **0.662**         | -                           |
| gnomAD constraint (1 - o/e, 10-nt window)     | 0.766            | +0.044 [-0.123, +0.112]    | 0.658             | +0.004 [-0.112, +0.111]     |
| gnomAD allele density (+/-2 nt)               | 0.686            | +0.124 [-0.137, +0.182]    | 0.597             | +0.065 [-0.038, +0.172]     |
| Trained secondary-structure model             | 0.665            | +0.145 [+0.068, +0.213]*   | 0.621             | +0.041 [-0.026, +0.116]     |
| CADD (PHRED)                                  | 0.655            | +0.121 [+0.021, +0.241]*   | 0.625             | +0.037 [-0.083, +0.162]     |
| Paralogue pathogenic density (MitoTIP-style)  | 0.615            | +0.195 [+0.136, +0.310]*   | 0.743             | -0.081 [-0.183, +0.013]     |
| gnomAD allele density (position)              | 0.609            | +0.201 [-0.064, +0.259]    | 0.588             | +0.073 [-0.021, +0.169]     |
| phyloP (447-way)                              | 0.593            | +0.217 [+0.079, +0.282]*   | 0.639             | +0.023 [-0.035, +0.081]     |
| Distance to nearest known pathogenic position | 0.554            | +0.256 [+0.185, +0.438]*   | 0.669             | -0.007 [-0.109, +0.086]     |
| Rfam element-type enrichment                  | 0.518            | +0.292 [+0.216, +0.363]*   | 0.618             | +0.044 [-0.057, +0.146]     |
| RNA-FM (foundation model, zero-shot)          | 0.460            | +0.350 [+0.143, +0.408]*   | 0.623             | +0.007 [-0.085, +0.095]     |
| ViennaRNA pairing probability                 | 0.434            | +0.376 [+0.299, +0.534]*   | 0.515             | +0.147 [+0.017, +0.284]*    |
| ViennaRNA fold ddG                            | 0.260            | +0.549 [+0.353, +0.646]*   | 0.519             | +0.143 [+0.009, +0.277]*    |
| ViennaRNA duplex ddG                          | 0.251            | +0.598 [+0.194, +0.672]*   | 0.636             | +0.056 [-0.077, +0.191]     |

Mean AUROC across held-out units, by disease mechanism. dAUROC is MuSE minus the baseline, averaged over units, with a paired position-block bootstrap 95% CI (2,000 resamples). * marks intervals excluding zero. Hierarchical (unit-then-position) intervals are in Supplementary Table S1 and are wider throughout, reflecting that only four units contribute to each context.
