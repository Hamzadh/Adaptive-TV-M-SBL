# Adaptive-TV-M-SBL
This repository contains the official implementation of our paper "Adaptive and Self-Tuning SBL with Total Variation Priors for Block-Sparse Signal Recovery."
The contetnts includes complete MATLAB implementation of the proposed algorithm
as well as description for implementaion details.

All algorithms used in this paper are implemented in MATLAB 2020b. The implementation of the MSBL-DoL algorithm, provided by the authors, is available at \url{https://github.com/adityasant/TV-SBL}. The remaining algorithms were implemented from scratch.   
# Practical Implementaion 

We initialize the auxiliary variable as \(\vec{C}^{(0)} = \vec{0}_N\) and the dual variable matrix as \(\lambdab^{(0)} = \vec{0}_N\). Additionally, for the M-SBL algorithm, we set \(\gammab^{(0)} = 0.1 \times \vec{1}_N\). In both the proposed algorithm and MSBL-DoL, \(\gammab^{(0)}\) can be initialized similarly or derived from the SBL solution—either upon convergence or after a few iterations—as a warm start.  
