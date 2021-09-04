Standalone version of Hybrid Triangular Solver, from
    https://github.com/trilinos/Trilinos/tree/master/packages/shylu/shylu_node/hts

Edit Makefile to include your system details. Run make. The output is
./hts_test, the unit tester. Run like this:
    OMP_NUM_THREADS=16 ./hts_test

The unit test code also serves as an example of how to interface with HTS.
