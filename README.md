# A framework for machine learning of climate model physics
This repository contains the code that allows to introduce ML models that were trained using PyTorch
into a climate model. The repository consists of a module that creates a shared object that is climate 
model (aka GCM) agnostic, and a reference implementation into CESM and the Community Atmospheric Model
(CAM).

The goal of the project is to create a framework that facilitates the entire development/training/testing/running of ML-based parametrizations in GCMs. It aims to reduce the development cycle allowing the ML model to run in a compiled form (i.e. LibTorch) while training and evaluating the ML-model in an interpreted environment. 



## Supported functionality


The framework supports the following functionality:


## Building and using the LibTorch plugin
### Building the plugin


## Building and using the CESM/CAM reference implementation
## Building the reference implementation


## Exporting PyTorch models into the framework


## Contributing
### Wishlist and gaps
### Plugging the framework into a new GCM/physics module
