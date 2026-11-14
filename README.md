# 🧬 ODS — Ordering DNA Sequences (`ODS-Ordering-DNA-Sequences`)

**Combinatorial Bioinformatics Optimization: Ordering DNA Sequences by Minimizing Pairwise Manhattan ($L_1$) Distance via Simulated Annealing**  
**Author:** Rinaldo Di Giorgio ([DiGiorgio Labs — `digiorgiolabs.com`](https://digiorgiolabs.com))  
**License:** Public Domain ([The Unlicense](LICENSE)) / MIT

---

## Overview

`ODS-Ordering-DNA-Sequences` applies **Simulated Annealing** to find a global ordering $\pi = (\pi_1, \pi_2, \dots, \pi_n)$ of genomic DNA sequences $S_1, \dots, S_n$ (represented as $k$-mer frequency vectors in $\mathbb{Z}^{4^k}$) that minimizes total pairwise **Manhattan ($L_1$) distance** along the sequence chain:

$$\min_{\pi} \sum_{i=1}^{n-1} \left\| \mathbf{v}(\pi_i) - \mathbf{v}(\pi_{i+1}) \right\|_1$$

### Quick Start

```bash
javac OrderingDnaSequences.java
java OrderingDnaSequences
```

---

## 🔓 DiGiorgio Labs Open-Source & Public Domain Catalog

* 🌐 **[View Full Interactive Open-Source Showcase on `digiorgiolabs.com/open-source` →](https://digiorgiolabs.com/open-source)**
* ⭐ **[Browse All Original Open-Source Repositories on GitHub →](https://github.com/digitalrinaldo?tab=repositories&type=source)**
* 📄 See [`OPEN_SOURCE_CATALOG.md`](OPEN_SOURCE_CATALOG.md) for the complete directory of DiGiorgio Labs open-source projects, patents, and engineering chronicles.
