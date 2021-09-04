# SpTRSV Split Execution Framework
The sparse triangular solve split execution framework for CPU-GPU heterogeneous systems can automatically determine the suitability of SpTRSV for a given input matrix for split-execution, find the appropriate split point, and execute SpTRSV in a split fashion using two SpTRSV algorithms while automatically managing any required inter-platform communication.  The model is implemented as a C++/CUDA library supporting multiple CPU-GPU algorithms.

For more details and performance evaluation and analysis results of the framework, refer to our paper titled ["A Split Execution Model for SpTRSV"](https://ieeexplore.ieee.org/document/9409717)

## Pre-requisites
Following are the pre-requisites to use or test the framework:

### 1) Hardware Platform
The framework is designed to work on heterogeneous CPU-GPU systems with current support for Intel CPUs and NVIDIA GPUs. In [1], the framework has been evaluated on two heterogeneous machines, (i) Intel Xeon Gold (6148) CPU with NVIDIA Tesla V100 GPU, and (ii) Intel Core I7 (8700K) CPU with NVIDIA GTX1080 Ti GPU. 

### 2) Software Tools
The framework requires following software tools to be available (installed) on the machine for compilation and execution.

#### i) OS:
The framework has been evaluated on CentOS (release 7.4.1708) and Ubuntu 18.04.1

#### i) libUFget Library:
The [libUFget library](https://zenodo.org/record/897632#.YTOngykzZuQ) provides a C interface to SuiteSparse Matrix Collection (previously University of Florida Sparse Matrix Collection [2]). The framework uses the library to automatically download matrices from the SuiteSparse Matrix Collection. Follow the installation instructions provided in INSTALL.md in the main directory of libUFget library (https://github.com/ParCoreLab/Split_SpTRSV/extras/libUFget-1.0.3/INSTALL.md). 

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


## References

[1] N. Ahmad, B. Yilmaz and D. Unat, "A Split Execution Model for SpTRSV," in IEEE Transactions on Parallel and Distributed Systems, vol. 32, no. 11, pp. 2809-2822, 1 Nov. 2021, doi: 10.1109/TPDS.2021.3074501.

[2] T. A. Davis and Y. Hu, "The University of Florida sparse matrix collection", ACM Trans. Math. Softw., vol. 38, no. 1, Dec. 2011.


