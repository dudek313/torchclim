# A framework for machine learning of climate model physics
This repository contains the code that allows to introduce ML models that were trained using PyTorch
into a climate model. The repository consists of a module that creates a shared object that is climate 
model (aka GCM) agnostic, and a reference implementation into CESM and the Community Atmospheric Model
(CAM).

The goal of the project is to create a framework that facilitates the entire development/training/testing/running of ML-based parametrizations in GCMs. It aims to reduce the development cycle allowing the ML model to run in a compiled form (i.e. LibTorch) while training and evaluating the ML-model in an interpreted environment. 



## Supported functionality


## Exporting PyTorch models into the framework


## Building the LibTorch plugin
The build process would produce the following by-products that should be incorporated to your GCM:
1. A shared library: "libtorch-plugin.so"  
2. A fortran module: "torch_plugin.mod"


Build steps:  
1. clone this repository:  
git clone git@github.com:dudek313/climate-model-physics-ml.git
2. cd climate-model-physics-ml/torch-wrapper
3. Download the LibTorch dependencies by running:  
./envs/nci-gadi/install-deps.sh  
4. Build:  
./build.sh  
5. Verify that "build/src/interface/libtorch-plugin.so" was created.



## Building and using the CESM/CAM reference implementation
Here we demonstrate the process on CESM 1.0.6 by creating an AMIP scenario.

```bash
cd [to/csem/scripts/directory]  
export CASE_PATH=[path/to/CSEM/scenario]  
./create_newcase -case $CASE_PATH -compset F_AMIP_CN -res f19_f19 -mach gadi  

cd $CASE_PATH/

git clone git@github.com:dudek313/climate-model-physics-ml.git
#build the torch wrapper (see "Building the LibTorch plugin")

#edit the "env_mach_specific" file
#add the path to the plugin shared object to LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH `pwd`/climate-model-physics-ml/torch-wrapper/build/src/interface/:$LD_LIBRARY_PATH

#once all the env_* files are set
./configure -case
cp climate-model-physics-ml/cesm-module/src.cam/* SourceMods/src.cam/

#edit Macros.gadi
#add include fortran mod path
#e.g.
MOD_NETCDF    := $(NETCDF_PATH)/include/Intel /[absolute path]/climate-model-physics-ml/torch-wrapper/build/src/interface/

#edit Macros.gadi
#add include path
#e.g.
INCLDIR := -I. -I/[absolute path]/climate-model-physics-ml/torch-wrapper/build/src/interface/

#edit Tools/Makefile
#add libtorch-plugin.so to the linking stage
ULIBS += -L$(LIBROOT) -lcsm_share -lmct -lmpeu -lpio -L/[absolute path]/climate-model-physics-ml/torch-wrapper/build/src/interface/ -ltorch-plugin

#edit SourceMods/src.cam/machine_learning_model_config.F90
#modify the configuration options to suite your needs
vim SourceMods/src.cam/machine_learning_model_config.F90

```





## Contributing
### Wishlist and gaps
### Plugging the framework into a new GCM/physics module
