#CXX      := $(CXX)
HCXX     := icpc
CXX      := nvcc
CXXFLAGS := -std=c++14 
CXXFLAGS += --compiler-options -mkl=sequential,-std=c++14,-Wall #-DCUSPARSE_ONLY #-DAUTO_TESTING #-DSAVE_MATRIX, #-DSTORE_BINARY, #,-DSAVE_MATRIX #-I${MKLROOT}/include
HCXXFLAGS := -std=c++14 -fPIC 
LDFLAGS  := -lufget -lcusparse #-lsqlite3 -larchive -lz -lbz2 -llzma -lmatio -lcurl -lssl
#LDFLAGS  += -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
LDFLAGS  +=  --linker-options -lpthread,-lm,-ldl  #--linker-options ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a,${MKLROOT}/lib/intel64/libmkl_sequential.a,${MKLROOT}/lib/intel64/libmkl_core.a,-lpthread,-lm,-ldl
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)
TARGET   := silu_test
INCLUDE  := -Iinclude/
SRC      :=                      \
   $(wildcard src/*.cpp)         \
   $(wildcard synk/*.cpp)        \

SRC_CUDA :=                      \
   $(wildcard src/*.cu)         \



OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
OBJECTS_CUDA := $(SRC_CUDA:%.cu=$(OBJ_DIR)/%.o)
OS := $(shell uname)
TARGET_SIZE=64
NVCCFLAGS   := -O3 -m$(TARGET_SIZE) -w -Xptxas -dlcm=cg
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_LDFLAGS += $(ALL_CCFLAGS)
CCBIN := -ccbin icpc 
################################################################################

# Gencode arguments
SMS ?= 61 70

ifeq ($(SMS),)
$(info >>> WARNING - no SM architectures have been specified - waiving sample <<<)
SAMPLE_ENABLED := 0
endif

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))

# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif
endif

LIBRARIES += -lcublas -lcusparse -lufget

ifeq ($(SAMPLE_ENABLED),0)
EXEC ?= @echo "[@]"
endif

# Generated via https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
#ifeq ($(OS), Darwin)
#	LDFLAGS += $(MKLROOT)/lib/libmkl_intel_ilp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl
#	CXXFLAGS += -DMKL_ILP64 -I$(MKLROOT)/include
#else ifeq ($(OS), Linux)
#	LDFLAGS += -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#	CXXFLAGS +=  -DMKL_ILP64 -I$(MKLROOT)/include
#endif

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp 
	@mkdir -p $(@D)
	$(CXX) $(CCBIN) $(CXXFLAGS) $(ALL_CCFLAGS) $(INCLUDE) -o $@ -c $<

$(OBJ_DIR)/%.o: %.cu
	@mkdir -p $(@D)
	$(CXX) $(CCBIN) $(CXXFLAGS)  $(GENCODE_FLAGS) $(ALL_CCFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS) $(OBJECTS_CUDA)
	@mkdir -p $(@D)
	$(CXX) $(CCBIN) $(CXXFLAGS) $(ALL_LDFLAGS) $(INCLUDE) $(OBJECTS) $(OBJECTS_CUDA) $(LDFLAGS) -o $(APP_DIR)/$(TARGET)

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*
