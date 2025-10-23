# Computational pipeline of de novo Zn proteases
Anqi Chen, Institute for Protein Design, University of Washington, Oct.23 2025

Python dependencies are listed in `requirements.txt`.

The jupyter notebook `zn_protease.ipynb` demonstrates all computational steps in one full design round, including backbone generation using RoseTTAfold Diffusion 2 -- Molecular Interfaces (RFD2-MI), backbone benchmarking and filtering, sequence design using EnhancedMPNN(2)-FastRelax cycles, structure prediction using AlphaFold3 (AF3). Python scripts for command generation, analysis, filtering used in this procedure are included in this repository. 

The source code of RFD2-MI will be made available after the publication of Bauer,2025(1). 
EnhancedMPNN is available through https://github.com/baker-laboratory/fused_mpnn

Inference calculations for RFD2-MI and AF3 were performed using NVIDIA A4000, A6000, L40, and H100 GPUs on a shared computing cluster. CUDA version 13.0 and NVIDIA driver version  580.82.07 were used for GPU computations.
Sequence design and filtering steps were executed on CPU nodes equipped with AMD EPYC 7702P 64-Core processors.


## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.


## Reference
1. Bauer, M. S. et al. De novo design of phospho-tyrosine peptide binders. Preprint at https://doi.org/10.1101/2025.09.29.678898 (2025).
2. Xue, F. et al. Improving protein sequence design through designability preference optimization. Preprint at http://arxiv.org/abs/2506.00297 (2025).


