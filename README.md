# SpTRSV Split Execution Framework
The sparse triangular solve split execution framework for CPU-GPU can automatically determine the suitability of an SpTRSV for split-execution, find the appropriate split point, and execute SpTRSV in a split fashion using two SpTRSV algorithms while automatically managing any required inter-platform communication.  The model is implemented as a C++/CUDA library supporting multiple CPU-GPU algorithms.

For more details and performance evaluation and analysis results of the framework, refer to our paper titled ["A Split Execution Model for SpTRSV"](https://ieeexplore.ieee.org/document/9409717)

## Pre-requisites
Following are the pre-requisites to use or test the framework:

### 1) Hardware Platform
The framework is designed to work on heterogeneous CPU-GPU systems with current support for Intel CPUs and NVIDIA GPUs. In [1], the framework has been evaluated on two heterogeneous machines, (i) Intel Xeon Gold (6148) CPU with NVIDIA Tesla V100 GPU, and (ii) Intel Core I7 (8700K) CPU with NVIDIA GTX1080 Ti GPU. 

### 2) Software Tools
The framework requires following software tools to be available (installed) on the machine for compilation and execution.

#### i) libUFget Library: 



## References

[1] N. Ahmad, B. Yilmaz and D. Unat, "A Split Execution Model for SpTRSV," in IEEE Transactions on Parallel and Distributed Systems, vol. 32, no. 11, pp. 2809-2822, 1 Nov. 2021, doi: 10.1109/TPDS.2021.3074501.


