# Optimization methods in cardiac MRI （ECE 5759 Final Project）
### In this project, I implemented and compared several optimization methodes for parallel imaging and compressive sensing. The numerical results are as expected: the momentum term and over-relaxation term can accelerate the convergence significantly with negligible additional computation.
## Parallel imaging (SENSE based)
Least square problem
<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/unifrom_samp.png" width="700">
<img src="https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/PI_results.png" width="700">
Figure 2: Reconstructed cardiac images (systolic frame), the error maps and the convergent speeds with different optimization methods (GM, FGM, OGM) and acceleration rates (R = 2, 4, 6) for parallel imaging (PI). The stepsize α= 1/L and the total number of iteration N=  150.  As expected, the convergence speed of FGM and OGM are much faster than GM, and OGM has the best performance. For R= 4,6, GM didn’t converge within 150 iterations. For R= 6, FGM and OGM diverged due to the high acceleration rate.
## Compressive sensing (SENSE based)
Lasso problem (l1 penalty, sparsifying transform: Temporal Fourier Transform)
![](https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/random_samp.png)
![](https://github.com/MRIOSU/ECE5759_Project/blob/main/Results/CS_resluts.png)
