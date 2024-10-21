**Note: If your work uses this algorithm or makes improvements based on it, please be sure to cite this paper. Thank you for your cooperation.**
___
**注意：如果您的工作用到了本算法，或者基于本算法进行了改进，请您务必引用本论文，谢谢配合。**
___
The core algorithm of paper "Regularized reconstruction based on joint smoothly clipped absolute deviation regularization and graph manifold learning for fluorescence molecular tomography"

Jun Zhang, Gege Zhang, Yi Chen, Kang Li, Fengjun Zhao, Huangjian Yi, Linzhi Su and Xin Cao※.

Physics in Medicine & Biology.(2023).

Environments

Matlab 2020

Instructions for use:

Based on the approximate simplified equations and boundary conditions of RTE, the finite element method is used to obtain AX=B. Among them, A is the system matrix and B is the observation vector.

Put A and B into the algorithm and solve to obtain X. Note that in order to obtain a more accurate X, some parameters in the algorithm may need to be adjusted.

This algorithm utilizes graph manifold learning to explore the latent features of the reconstruction source, enhancing the morphological similarity of the algorithm. Therefore, it is necessary to modify the algorithm based on the specific model when computing.Specifically, during the calculation of the graph Laplacian matrix, the structural information of the model is required, which involves determining whether different nodes belong to the same organ. This part requires specific processing based on the particular model.