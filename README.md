# NeuroProt-ASD

**A Paralog-Safe Multi-Modal Benchmark for Autism Risk Gene Prioritization and In-Vivo Validation Confounder Analysis**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Pytest: 7/7 passing](https://img.shields.io/badge/pytest-7%2F7%20passing-brightgreen.svg)](tests/)

---

## 🔬 Overview

**NeuroProt-ASD** is a machine learning benchmarking framework designed to evaluate the true translational utility of protein language models (PLMs) in psychiatric genetics while eliminating two critical methodological vulnerabilities:
1. **Homology Leakage:** Cross-validation is strictly grouped across **430 distinct paralog families**, preventing sequence similarity between training and testing sets.
2. **Confounder Mediation:** Systematic parametric and non-parametric evaluation of in-vivo mammalian knockout validations against baseline evolutionary constraint ($p\mathrm{LI}$) and gene length.

### Core Biological Modalities Integrated:
1. **Molecular & Structural PLMs:** Real canonical Swiss-Prot proteome embeddings across three ESM-2 scales (8M, 35M, 150M parameters) and SaProt 3D Foldseek tokens.
2. **Spatiotemporal Fetal Transcriptomics:** BrainSpan human neocortex developmental windows.
3. **Functional Interactome Networks:** STRING protein-protein interactions and Krishnan et al. brain-specific functional networks.
4. **Statistical Genetics & Evolutionary Constraint:** TADA de novo mutation Bayes Factors, gnomAD $p\mathrm{LI}$ loss-of-function intolerance, and baseline mutation rates.

---

## 📊 Key Benchmark Findings

### 1. Paralog-Safe Modality Ablation (Table 1)
| Modality / Model Framework | Paralog-Safe AUROC | 95% Bootstrap CI |
|:---|:---:|:---:|
| **Confounders & Conservation Baseline** | **0.8905** | [0.8558 – 0.9219] |
| **Structural PLM-only (ESM-2 8M + SaProt 3Di)** | **0.6169** | [0.5583 – 0.6723] |
| **Spatiotemporal Transcriptomics (BrainSpan)** | **0.8048** | [0.7597 – 0.8476] |
| **Functional Interactome (STRING + Krishnan)** | **0.8775** | [0.8397 – 0.9126] |
| **Classical Integration (BrainSpan + Networks + Genetics)** | **0.9257** | [0.8918 – 0.9538] |
| **NeuroProt-ASD (Full Multi-Modal Fusion)** | **0.9193** | [0.8880 – 0.9509] |

### 2. The Foundation Model Scale Bottleneck
Scaling ESM-2 by nearly 20-fold (**8M: 0.6768**, **35M: 0.6248**, **150M: 0.6643**) yields a flat performance plateau without closing the gap with tissue transcriptomics (0.8048), establishing that the **tissue-context bottleneck is an invariant property of protein models** in complex psychiatric traits.

### 3. Fair Prior Art Re-evaluation (forecASD Parity)
When legacy state-of-the-art (**forecASD**, Brueggeman et al. 2020) is re-evaluated under our identical 430-family paralog-safe GroupKFold protocol:
- Published in-sample score of **0.9510** drops to **0.9042** (Level-1 BrainSpan+STRING RF) and **0.9228** (Full Stack Ensemble).
- **NeuroProt-ASD (0.9193)** exhibits statistical parity with forecASD Full Stack (paired bootstrap $p = 0.7210$, $\Delta = -0.0034$), while numerically exceeding Level-1 RF (+0.0151 AUROC).

### 4. Methodological Warning: In-Vivo Knockout Mediation
While novel candidates show significant univariable knockout enrichment in the International Mouse Phenotyping Consortium (IMPC; **47.29% hit rate vs 36.95% background, OR = 1.531, $p = 0.0109$**), multivariable logistic regression and an independent **1,000-iteration matched-null permutation test ($p = 0.5390$)** confirm this enrichment is **fully mediated by evolutionary constraint ($p\mathrm{LI}$)** rather than autism-specific biology.

---

## 🚀 Getting Started

### Installation
```bash
git clone https://github.com/tanukusaitejesh-prog/NeuroProt-ASD.git
cd NeuroProt-ASD
pip install -r requirements.txt  # or install torch, scikit-learn, transformers, pandas, scipy
```

### Reproduce Master Pipeline
```bash
# Runs full paralog-safe pipeline on CUDA / CPU
python run_neuroprot_pipeline.py
```

### Run Test Suite
```bash
python -m pytest tests/test_pipeline.py -v
```

---

## 📄 Full Manuscript Report
The complete manuscript including figures, case studies (*ARID2*, *ATN1*, *ITSN1*, *ADAM22*, *HDAC9*, *CACNA1E*), and full regression tables is available in [`FINAL_REPORT.md`](FINAL_REPORT.md).

---

## 📜 License
This project is open-source under the MIT License.
