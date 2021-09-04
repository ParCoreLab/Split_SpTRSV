The UF Sparse Collection C interface                                 {#mainpage} 
====================================
 
The UF Sparse Collection ( http://www.cise.ufl.edu/research/sparse/ ) is a huge
collection of sparse matrices for academic research and benchmarking of
algorithms such a linear equation solvers or eigenvalue solvers. The collection 
is normally used via a web interface or a MATLAB script which accesses the
collection and download matrices by their ID or their name. 

The libUFget library provides a C interface to the UF Sparse Collection which
allows to access the collection form C codes and to download the matrices as in
the MATLAB interface. Furthermore, it converts the database into an SQLITE
database which than can be query via the SQL language to select a set of
matrices according to the requirements of the user. Using an iterator it is
possible to perform benchmarks for all matrices matched by the query. 

Additionally the library provides function to read the matrices into the
coordinate storage format as well as function to convert them into Compressed
Sparse Row (CSR) and Compressed Sparse Column (CSC) storage. All functions
support 32 bit integers as well as 64 bit integers to allow even to use the
largest matrices currently available in the collection and to support 64 bit
integer enabled codes directly without further data manipulation. 

Installation
------------
See [INSTALL.md](INSTALL.md) for details. 

Usage 
-----
The library requires one header file 
~~~~~~~~~~~~~~~~~~~~~
#include <libufget.h> 
~~~~~~~~~~~~~~~~~~~~~
and to linkage of your application against *libufget.so* via 

    gcc yourapp.c -o yourapp -lufget 

The library is split into three categories: 
* [Interface to the collection](@ref collection) 
* [Query the collection](@ref query)
* [Download, Read and Convert Sparse Matrices](@ref matrix) 

The ufget-update Tool
---------------------
The library also provides a tool called *ufget-update* which manages the local
matrix cache from the command line. The tool allows to update the collection
database via:

    ufget-update update

It can download the whole collection or single matrices by running 

    ufget-update download-all 
    ufget-update download GROUPNAME MATRIXNAME

Finally it can clean up the collection, which removed everything from the cache
except of the database file by 

    ufget-update clean 

The tool is installed to the standard binary directory of the given installation
prefix on `make install`. 

Environment 
-----------
The following environment variables control the behavior of the library: 
* If **UF_COLLECTION_VERBOSE** is not zero it makes the library verbose and
  printing debug messages.
* If **UF_COLLECTION_CLEANUP** is not zero the library cleans up its cache at the
  end of the program. All downloaded matrices will be remove from the disk 
  cache and only the data base is retained. 
* If **UF_COLLECTION_BASEURL** is set the value is used as an alternative base url
  for the collection. This is useful if the UF collection moves or one hosts its
  own one. By default http://www.cise.ufl.edu/research/sparse/ is used as base
  url. 
* If **UF_COLLECTION_CACHE_DIR** is set the library uses the given path as cache
  directory. Otherwise the default path $HOME/.ufget/ is used.
						 
Further information about these environment variable can be found in the
documentation of [uf_collection_init](@ref uf_collection_init). 

Examples
--------
Documented examples can be found in the *examples/* directory. 

License
-------
The whole library is provided under the conditions of the GPLv2 or higher if you
want. 
