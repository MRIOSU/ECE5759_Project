# Optimization methods in cardiac MRI （ECE 5759 Final Project）
### In this project, I implemented and compared several optimization methodes for parallel imaging and compressive sensing in cardiac MRI. The numerical results are as expected: the momentum term and over-relaxation term can accelerate the convergence significantly with negligible additional computation.
## Parallel imaging (SENSE based)
Least square problem

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/unifrom_samp.png" width="500">

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/PI_results.png" width="500">

Figure 2: Reconstructed cardiac images (systolic frame), the error maps and the convergent speeds with different optimization methods (GM, FGM, OGM) and acceleration rates (R = 2, 4, 6) for the parallel imaging (PI). The stepsize α= 1/L and the total number of iteration N =  150.  The convergence speed of FGM and OGM are much faster than GM, and OGM has the best convergence performance. For R= 4,6, GM didn’t converge within 150 iterations. For R= 6, FGM and OGM diverged due to the high acceleration rate (optimization problem is ill posed).
## Compressive sensing (SENSE based)
Lasso problem (l1 penalty, sparsifying transform: Temporal Fourier Transform)

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/random_samp.png" width="500">

<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/CS_resluts.png" width="500">

Figure 4: Reconstructed cardiac images (systolic frame), the error maps and the convergent speedswith different optimization methods (ISTA, FISTA, POGM) and acceleration rates (R = 4, 6, 8) for compressive sensing (CS). The stepsize α= 1/L and the total number of iteration N= 150. The convergence speed of FISTA and POGM are much faster than ISTA, and POGM hasthe best convergence performance. For R= 8, ISTA didn’t converge within 150 iterations.
