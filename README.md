# Computational pipeline of de novo Zn proteases

**Anqi Chen, Institute for Protein Design, University of Washington, Oct. 23, 2025**

The repository for the computational design of de novo Zn proteases. The Jupyter notebook `zn_protease.ipynb` demonstrates one design round, including:

- Backbone generation using RoseTTAFold Diffusion 2 — Molecular Interfaces (RFD2-MI)
- Backbone benchmarking and filtering
- Sequence design using EnhancedMPNN(2) + FastRelax cycles
- Structure prediction using AlphaFold3 (AF3)

Python scripts for command generation, analysis, and filtering used in this procedure are included in this repository. RFD2-MI source code will be released after the publication of Bauer, 2025 (1). EnhancedMPNN is available at [fused_mpnn](https://github.com/baker-laboratory/fused_mpnn).

---

## Hardware & Software

- **GPU**: Inference computation for RFD2-MI and AF3 were performed using NVIDIA A4000, A6000, L40, H100 on a shared computing cluster with CUDA 13.0, NVIDIA driver 580.82.07. Inference time for RFD2-MI is 1-2 mins for each backbone on an A4000 gpu. For AF3 is ~ 1 min for each sequence.   
- **CPU**: Sequence design, analysis and filtering were performed on AMD EPYC 7702P 64-Core cpu nodes. The cpu run time for the sequence design of one backbone in 10 cycles of EnhancedMPNN-FR is around 10 mins. Using the parallelization parameters specified in the Jupyter notebook `zn_protease.ipynb`, the cpu run time for each job in the analysis, and filtering steps are less than 10 mins. 
- **Python environment**: See `requirements.txt` 

---

## Dependencies

Install required Python packages using (~ 10 mins, the time for the installation of RFD2-MI, EnhancedMPNN, AF3 is not included):

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


