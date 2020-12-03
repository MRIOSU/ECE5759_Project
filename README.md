# Optimization methods in cardiac MRI （ECE5759 Final）
### In this project, I implemented and compared several optimization methods for parallel imaging and compressive sensing in cardiac MRI. The numerical results are as expected: the momentum term and over-relaxation term can accelerate the convergence significantly with negligible additional computation.

#### Project report: https://github.com/MRIOSU/ECE5759_Project/blob/main/ECE5759_Project_Chen.pdf
## Parallel imaging (SENSE based)
Modeled as Least square problem. The gradient method (GM), fast GM (FGM) and optimized GM (OGM) were compared.

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/unifrom_samp.png" width="500">

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/PI_results.png" width="500">

Figure 2: Reconstructed cardiac cine images (only systolic frame shown here), the error maps and the convergent speeds using different optimization methods (GM, FGM, OGM) and acceleration rates (R = 2, 4, 6) for parallel imaging (PI). The stepsize α= 1/L and the total number of iteration N =  150.  The convergence speed of FGM and OGM are much faster than GM, and OGM has the best convergence performance. For R= 4,6, GM didn’t converge within 150 iterations. For R= 6, FGM and OGM diverged due to the high acceleration rate (optimization problem is ill posed).
## Compressive sensing (SENSE based)
Modeled as Lasso problem (l1 penalty, sparsifying transform: Temporal Fourier Transform). ISTA/PGM, FISTA/FPGM and POGM were compared.

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/random_samp.png" width="500">

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/CS_resluts.png" width="500">

Figure 4: Reconstructed cardiac cine images (only systolic frame shown here), the error maps and the convergent speeds using different optimization methods (ISTA, FISTA, POGM) and acceleration rates (R = 4, 6, 8) for compressive sensing (CS). The stepsize α= 1/L and the total number of iteration N = 150. The convergence speed of FISTA and POGM are much faster than ISTA, and POGM has the best convergence performance. For R= 8, ISTA didn’t converge within 150 iterations.
