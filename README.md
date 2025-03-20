# Adaptive-TV-M-SBL
This repository contains the official implementation of our paper "Adaptive and Self-Tuning SBL with Total Variation Priors for Block-Sparse Signal Recovery."
The contetnts includes complete MATLAB implementation of the proposed algorithm
as well as description for implementaion details.

All algorithms used in this paper are implemented in MATLAB 2020b. The implementation of the MSBL-DoL algorithm, provided by the authors, is available at \url{https://github.com/adityasant/TV-SBL}. The remaining algorithms were implemented from scratch.   
# Practical Implementaion 
## Initialization
We initialize the auxiliary variable as \(\vec{C}^{(0)} = \vec{0}_N\) and the dual variable matrix as \(\lambdab^{(0)} = \vec{0}_N\). Additionally, for the M-SBL algorithm, we set \(\gammab^{(0)} = 0.1 \times \vec{1}_N\). In both the proposed algorithm and MSBL-DoL, \(\gammab^{(0)}\) can be initialized similarly or derived from the SBL solution—either upon convergence or after a few iterations—as a warm start.  The stopping criterion is set to \( \epsilon=10^{-3} \), with maximum iterations of \( t_{\mathrm{max}}=10 \) for ADMM and \( k_{\mathrm{max}}=100 \) for EM. The number of considered neighbors is \( \{|\Omega_i|=3\}_{i=1}^N \). 
## Support detection
Given an arbitrary estimate \(\vec{X}\) from any sparse recovery algorithm, the estimated support of the signal is defined as  
\[
\hat{\mathcal{A}} = \{i \mid \|\hat{\vec{x}}_i\|_2 > \epsilon_{\mathrm{thr}}, \; \forall i \in \mathcal{N} \},
\]  
where \(\epsilon_{\mathrm{thr}}\) is a small pre-defined positive threshold. In our simulations, we set \(\epsilon_{\mathrm{thr}} = 0.1\). This threshold is not fixed but determined empirically to optimize performance based on signal values. Alternatively, the support can be estimated by selecting the \(K_{nz}\) largest elements of the signal, where \(K_{nz}\) denotes the number of nonzero elements.  
