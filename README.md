# Computational pipeline of de novo Zn proteases

**Anqi Chen, Institute for Protein Design, University of Washington, Oct. 23, 2025**

The repository for the computational design of de novo Zn proteases. The Jupyter notebook `zn_protease.ipynb` demonstrates one design round, including:

- Backbone generation using RoseTTAFold Diffusion 2 — Molecular Interfaces (RFD2-MI)
- Backbone benchmarking and filtering
- Sequence design using EnhancedMPNN(2) + FastRelax cycles
- Structure prediction using AlphaFold3 (AF3)

Python scripts for command generation, analysis, and filtering used in this procedure are included in this repository.

**Availability of tools**:

- RFD2-MI source code will be released after the publication of Bauer, 2025 (1).  
- EnhancedMPNN is available at [fused_mpnn](https://github.com/baker-laboratory/fused_mpnn).

---

## Hardware & Software

- **GPU inference**: NVIDIA A4000, A6000, L40, H100 on a shared computing cluster  
  - CUDA 13.0, NVIDIA driver 580.82.07  
- **CPU sequence design, analysis and filtering**: AMD EPYC 7702P 64-Core nodes  
- **Python environment**: See `requirements.txt` 

---

## Dependencies

Install Python packages using:

```bash
pip install -r requirements.txt
```

---

## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Reference
1. Bauer, M. S. et al. De novo design of phospho-tyrosine peptide binders. Preprint at https://doi.org/10.1101/2025.09.29.678898 (2025).
2. Xue, F. et al. Improving protein sequence design through designability preference optimization. Preprint at http://arxiv.org/abs/2506.00297 (2025).


