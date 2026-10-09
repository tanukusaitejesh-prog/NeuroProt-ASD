**Table 3. Four evaluation regimes of increasing difficulty.**

| Regime                             | What                                                                                |    n | Positives   | AUROC                                   |
|:-----------------------------------|:------------------------------------------------------------------------------------|-----:|:------------|:----------------------------------------|
| Held-out gene, dominant mechanism  | Every variant of one gene withheld entirely                                         | 1211 | 72          | 0.766 [0.713, 0.811]                    |
| Held-out gene, recessive mechanism | Same, recessive positives                                                           |  566 | 177         | 0.591 [0.540, 0.639]                    |
| Held-out gene, pooled              | Both mechanisms, rank-pooled                                                        | 1777 | 249         | 0.636 [0.599, 0.672]                    |
| Zero-shot functional transfer      | RNU4-2 saturation genome editing; no functional measurement used in training        |  485 | continuous  | Spearman 0.405 [0.255, 0.535] (Table 4) |
| Prospective, timestamped           | Predictions deposited 2026-10-09 with SHA-256, evaluation protocol fixed in advance | 2777 | pending     | to be evaluated on future reports       |

Held-out-gene is strictly harder than the cluster split used by comparable methods papers: no variant, position or paralogue of the test gene is seen in training. The prospective set is deposited and hashed, so it cannot be revised after the fact.
