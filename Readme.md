# Digital Analysis of Fingerprint

## About

* Project realised by Pierre-Louis Cauvin, Ewen Lallinec, Tanguy Terrien, Jules Treton
* This project is also available at https://gitlab.ensimag.fr/fingerprint-modeling/fingerprint-modeling/-/tree/master/fingerprint
## Table of contents

* Install
* Project architecture
* Files
    * Headers
    * Source
    * Tests
    * Data_in

* Build
* Run

## Install

This project is realised using cmake, to install it in the current directory just type: cmake <path/fingerprint>

## Project architecture

The project has the following architecture:
``` bash
fingerprint
├── CMakeLists.txt
├── Readme.md
├── data_in
├── data_out
├── docs
├── include
├── src
└── tests

```

## Files 
The files are separated in headers, source and tests. The data used for the tests are in data_in.
### Headers
There are 8 headers in include :
``` bash
include
├── GeometricalWarps.h
├── Image.h
├── LinearFiltering.h
├── RegistrationOptimization.h
├── Restoration.h
├── SaveManyImages.h
├── ShowManyImages.h
└── Simulation.h
```
### Source
There are 7 source files in src:
``` bash
src
├── GeometricalWarps.cpp
├── Image.cpp
├── LinearFiltering.cpp
├── MorphologicalFiltering.cpp
├── RegistrationOptimization.cpp
├── Restoration.cpp
└── Simulation.cpp
```

### Tests 
There are 6 files of tests creating 6 different executables.
``` bash
tests
├── SimulationGeometry.cpp
├── SimulationLinearFiltering.cpp
├── SimulationPressure.cpp
├── SimulationRegistration.cpp
├── SimulationWetnessDryness.cpp
└── TestRestoration.cpp
```

### Data_in
All the files needed for the tests files are in this repository.

## Build
To build you need to do a make command in the directory were you installed the project.
You should get this structure for the build repository.

```bash
build
├── CMakeFiles
├── data_in -> symlink/fingerprint/data_in
├── data_out
├── docs
│   └── html  
├── src
└── tests
    ├── Geometry
    ├── LinearFiltering
    ├── Makefile
    ├── Pressure
    ├── Registration
    ├── Restoration
    └── WetnessDryness
```
## Run
In the tests directory you can execute each file manually and separately.
Each one of those file registers figures in data_out.

**Warning** 
The restoration executable takes a long time (few minutes) to execute do not execute it without being prepared!
