# SpTRSV Split Execution Framework
The sparse triangular solve split execution framework for CPU-GPU heterogeneous systems can automatically determine the suitability of SpTRSV for a given input matrix for split-execution, find the appropriate split point, and execute SpTRSV in a split fashion using two SpTRSV algorithms while automatically managing any required inter-platform communication.  The model is implemented as a C++/CUDA library supporting multiple CPU-GPU algorithms.

For more details and performance evaluation and analysis results of the framework, refer to our paper titled ["A Split Execution Model for SpTRSV"](https://ieeexplore.ieee.org/document/9409717) [1]

## Pre-requisites
Following are the pre-requisites to use or test the framework:

### 1) Hardware Platform
The framework is designed to work on heterogeneous CPU-GPU systems with current support for Intel CPUs and NVIDIA GPUs (Compute capability 6.0 or later). In [1], the framework has been evaluated on two heterogeneous machines, (i) Intel Xeon Gold (6148) CPU with NVIDIA Tesla V100 GPU, and (ii) Intel Core I7 (8700K) CPU with NVIDIA GTX1080 Ti GPU. 

### 2) Software Tools
The framework requires following software tools to be available (installed) on the machine for compilation and execution.

#### i) OS:
The framework has been evaluated on CentOS (release 7.4.1708) and Ubuntu 18.04.1

#### i) libUFget Library:
The [libUFget library](https://zenodo.org/record/3894753#.YTO8YykzZuQ) provides a C interface to SuiteSparse Matrix Collection (previously University of Florida Sparse Matrix Collection [2]). The framework uses the library to automatically download matrices from the SuiteSparse Matrix Collection. Follow the installation instructions provided in INSTALL.md in the main directory of libUFget library (https://github.com/ParCoreLab/Split_SpTRSV/extras/libUFget-1.0.3/INSTALL.md). 

#### ii) CUDA Toolkit
For compiling and running the code on NVIDIA GPUs, framework utilizes [NVIDIAs CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit). For evaluation of the framework [1], CUDA version 10.1 is used.

#### iii) Intel Compiler and Intel MKL library
The framework has been evaluated with Intel C++ compiler (icpc) and Intel MKL library available in Intel Parallel Studio 2019. For instructions on how to download and install Intel Parallel Studio 2019, you may follow instructions on ["How do I get an older version of an Intel® Software Development Product?"](https://software.intel.com/content/www/us/en/develop/articles/older-version-product.html).

## Using the Framework
After installing the required software tools on a supported CPU-GPU system, framework can be used as follows:

### - Downloading the repository:
To download the repository, use the command line:

```bash
git clone https://github.com/ParCoreLab/Split_SpTRSV
```
### - Compiling the code: 
To compile the code, run the make command in the main directory.

```bash
make
```

### - Running the code:
The framework can be run in one of the three modes, (1) Mode 1: SpTRSV executes in split execution mode, where the execution may be split between CPU or GPU, or run on a single platform, as determined by the framework, (2) Mode 2: SpTRSV executes on the GPU platform using ELMC algorithm [2], (3) Mode 3: SpTRSV executes on the CPU platform using Intel MKL library.

#### 1) Automatic matrix download
To run the framework in one of these modes for a square matrix (with a given ID) automatically downloaded from the SuiteSparse Matrix Collection, the command line looks like the following:

```bash
./silu_test matrix_id mode
```

For example, to run the framework in split SpTRSV execution mode for matrix "FullChip" (ID 2380), use the command line:

```bash
./silu_test 2380 1
```

Replace 1 with 2 or 3 to run unified SpTRSV for the matrix on GPU or CPU, respectively. 

#### 2) Using a matrix stored in Matrix Market Format
If instead of downloading the matrix, it is required to use an available matrix in matrix market format (.mtx extension), one may use the following command line:

```bash
./silu_test matrix_id mode filename.mtx
```

Here, filename.mtx is the filename of the matrix file in the Matrix Market Format. In this case, matrix_id is "don't care" and may be any integer number (e.g. 10). "mode" is same as previously described for "Automatic matrix download" case.

## Supported SpTRSV Algorithms:
In [1], framework is evaluated with split execution using ELMC, MKL and HTS SpTRSV algorithms. However, implementation for other algorithms is also available. These include Sync Free algorithm [4], SLFC algorithm (a predecessor of ELMC) [5], SpTRSV using cuSPARSE v2 [6] and ELMR (a row-wise version of ELMC) [3].

In the original version of ths split execution framework, the HTS library has been separately evaluated using the code available in (https://github.com/ParCoreLab/Split_SpTRSV/extras/hts) for the proof of concept. It is planned to integrate the HTS code into the framework. Moreover, we also plan to re-design the code so as to define separate classes for each SpTRSV algorithm and polymorphically select appropriate algorithm at runtime. This will make the code cleaner and more easier to understand, use and modify.    

Note: Thanks to Martin Köhler for making necessary changes to libUFget as per our request and releasing a newer version of libUFget (1.0.3).  

## References

[1] N. Ahmad, B. Yilmaz and D. Unat, "A Split Execution Model for SpTRSV," in IEEE Transactions on Parallel and Distributed Systems, vol. 32, no. 11, pp. 2809-2822, 1 Nov. 2021, doi: 10.1109/TPDS.2021.3074501.

[2] T. A. Davis and Y. Hu, "The University of Florida sparse matrix collection", ACM Trans. Math. Softw., vol. 38, no. 1, Dec. 2011.

[3] R. Li and C. Zhang, "Efficient parallel implementations of sparse triangular solves for GPU architectures", Proc. SIAM Conf. Parallel Process. Sci. Comput., pp. 106-117, 2020.

[4] W. Liu, A. Li, J. D. Hogg, I. S. Duff and B. Vinter, "Fast synchronization-free algorithms for parallel sparse triangular solves with multiple right-hand sides", Concurrency Computation: Practice Experience, vol. 29, no. 21, 2017.

[5] R. Li, "On parallel solution of sparse triangular linear systems in CUDA", CoRR, 2007.

[6] N. Corporation, "NVIDIA cuSPARSE library", 2021.




