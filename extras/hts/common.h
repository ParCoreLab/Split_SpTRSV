// common.h
// Created: 03-12-2019
// Author: Najeeb Ahmad

#include <iostream>
#include <string>
#include <cstring>
#include <time.h>
#include <math.h>

#include <sys/time.h>

#ifndef val_type
#define val_type double
#endif

#ifndef ind_type
#define ind_type int
#endif

#ifndef sz_type
#define sz_type int
#endif

#ifndef WARP_SIZE
#define WARP_SIZE   32
#endif

#ifndef WARP_PER_BLOCK
#define WARP_PER_BLOCK   16
#endif

#ifndef BLOCKDIM
#define BLOCKDIM    WARP_PER_BLOCK * WARP_SIZE
#endif

#define SUBSTITUTION_FORWARD  0
#define SUBSTITUTION_BACKWARD 1
#define SUBSTITUTION_BOTH     2
#define RANDOM_REF            3

using index_t = long long;
