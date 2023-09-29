
# TorchClim: A deep-learning framework for climate model physics

## Overview
TorchClim is a framework that allows the introduction of ML/AI models that were trained using PyTorch into a climate model (aka GCM). It facilitates a fast turnover of the train-test-run workflow allowing for quick development of ML/AI-based parametrizations into parallel and distributed environments.

The framework consists of two components: first, a plugin that is in charge of the interactions with the underlying ML/AI framework. Currently, this plugin relies on LibTorch as the underlying implementation. This plugin is largely agnostic to the specifics of the GCM that is using it. Second, we provide a reference implementation of the framework into CESM version 1.0.6 where we replace moist and radiative parametrizations in the Community Atmospheric Model (CAM) version 4 with an ML surrogate. The reference implementation offers a range of innovative features, facilitating many aspects of the train-test-run workflow of ML models for GCMs. For example, it offers tools to extract data from CAM before and after the desired point of insertion of a surrogate parametrization. It also offers the ability to switch parametrizations at runtime based on time-space requirements, etc. 

This repository contains two main folders:

1. torch-wrapper\
   This is the low-level plugin in charge of loading and calling to the underlying ML/AI model.
2. cesm-module/src.cam\
   The reference implementation demonstrates how to load the plugin and call it from a GCM.


## Usage
Assuming that you have a climate model or a distributed or parallel application, written in Fortran of c/c++, and you want to introduce an ML/AI model to it. After the installation of the framework, you will need to train your ML/AI model using PyTorch. Once the training is completed, export your surrogate model. Currently, TorchClim supports export as a torch script (with future extensions to support other frameworks aside from PyTorch via the ONYX interface). At this stage you will need to compile the TorchClim plugin, pointing it to the exported surrogate model. At this stage, you will need to load the TorchClim plugin into your GCM. See the reference implementation provided here for an example of how to achieve that.


## Prerequisites
TorchClim requires Fortran and c/c++ compilers. The reference implementation was tested using the intel compiler version 2021.5.0. You will need to edit "torch-wrapper/env/load-env.sh" to load the compilers of your choice. 

You will also need to edit "torch-wrapper/env/install-deps.sh" to point to the libtorch implementation that you wish to use. Ideally, this should be as close as possible to the version that is used during training. You will also need to choose between a GPU or CPU-based version, depending on your local infrastructure and the mod in which you which to use the surrogate model.


## Building the TorchClim plugin
The build process would produce the following by-products that should be incorporated into the GCM:
1. A shared library: "libtorch-plugin.so"  
2. A Fortran module: "torch_plugin.mod"

Before building the plugin, ensure:
1. That the plugin can find your surrogate model:\
   Edit the "torch_wrap_predict" function in "torch-wrapper/src/interface/torch-wrap.cpp" (see script_path variable). 
2. Any additional interfaces that you which to expose to the GCM are created:\
   Choose an existing interface from "torch-wrapper/src/interface/torch-plugin.f90" or create a new one.


Build steps:  
1. clone this repository:\  
   git clone git@github.com:dudek313/torchclim.git
2. cd torchclim/torch-wrapper
3. Download the LibTorch dependencies by running:\
   ./env/install-deps.sh  
4. Build:\  
   ./build.sh  
5. Verify that "build/src/interface/libtorch-plugin.so" was created.



## Building the CESM/CAM reference implementation with TorchClim
Here we demonstrate the process on CESM 1.0.6 by creating an AMIP scenario and running it on the Australian NCI/Gadi infrastructure (nci.org.au/).

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
cp torchclim/cesm-module/src.cam/* SourceMods/src.cam/

#edit Macros.gadi
#add include Fortran mod path
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
#modify the configuration options to suit your needs
vim SourceMods/src.cam/machine_learning_model_config.F90

```

## Configuration options of the CESM/CAM reference implementation

The configuration of the reference implementation is done through the "torchclim/cesm-module/src.cam/machine_learning_model_config.F90" module.
The following options are implemented:

1. ml_model_skip_first_steps: true/false, steps_to_skip: integer (model steps)\
   These two parameters enable startup using the original CAM parametrization (true by default for the first 24 hours).
3. ml_model_enabled: true/false\
   When set to .false., the original CAM parametrization will be used.
5. ml_model_standalone: true/false\
   When set to .false. alongside with 'ml_model_enabled' set to .true., both the original and ML surrogate will run side by side.
7. radiation_from_ml: true/false\
   Whether to override the ML radiation variables with the original CAM radiative parametrization.
9. output_state_vars_from_ctrl_run: true/false\
   Whether or not to output state variables before and after the intended point of insertion of the surrogate model. This allows to use data from the original CAM parametrization for the initial training of the surrogate model. 
11. ml_only_in_lat_band: true/false\
    Together with max_lat and min_lat, limit the surrogate model to a latitude band.





