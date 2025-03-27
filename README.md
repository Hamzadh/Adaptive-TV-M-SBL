# Adaptive-TV-M-SBL

## Overview
This repository contains the official implementation of our paper "Adaptive and Self-Tuning SBL with Total Variation Priors for Block-Sparse Signal Recovery." 

The code provides a complete MATLAB implementation of the proposed algorithm along with detailed implementation guidelines and examples.

## Requirements
- MATLAB 2020b or newer
- The implementation of the MSBL-DoL algorithm (used for comparison) is available at [https://github.com/adityasant/TV-SBL](https://github.com/adityasant/TV-SBL)
- All other algorithms referenced in the paper were implemented from scratch

## Implementation Details

### Initialization
For optimal performance, we initialize the algorithm with the following parameters:
- Auxiliary variable: $\boldsymbol{C}^{(0)} = \boldsymbol{0}_N$
- Dual variable matrix: $\vecgreek{\lambda}^{(0)} = \boldsymbol{0}_N$
- M-SBL algorithm parameter: $\vecgreek{gamma}^{(0)} = 0.1 \times \boldsymbol{1}_N$

For both our proposed algorithm and MSBL-DoL, $\vecgreek{\gamma}^{(0)}$ can be initialized as described above or derived from the SBL solution (either upon convergence or after a few iterations) for a warm start.

### Convergence Parameters
- Stopping criterion: $\epsilon=10^{-3}$
- Maximum ADMM iterations: $t_{\mathrm{max}}=10$
- Maximum EM iterations: $k_{\mathrm{max}}=100$
- Number of considered neighbors: $\{|\Omega_i|=3\}$

### Support Detection
Given an estimate $\vec{X}$ from any sparse recovery algorithm, we define the estimated support of the signal as:

$$\hat{\mathcal{A}} = \{i \mid \|\hat{\vec{x}}_i\|_2 > \epsilon_{\mathrm{thr}}, \; \forall i \in \mathcal{N} \}$$

where $\epsilon_{\mathrm{thr}}$ is a small pre-defined positive threshold. In our simulations, we use $\epsilon_{\mathrm{thr}} = 0.1$. This threshold is determined empirically to optimize performance based on the signal characteristics.

Alternatively, the support can be estimated by selecting the $K_{nz}$ largest elements of the signal, where $K_{nz}$ denotes the number of nonzero elements.

## Citation
If you use this code in your research, please cite our paper:(this will be updated once the paper is published)
```
@article{...}
```

## Canatct infomation

For any question, inquiry or feedback, Please get in touch via email: hamza.djelouat@oulu.fi or djelouat.hmz@gmail.com 

## License
[Include license information here]

## Contact
[Include contact information here]
