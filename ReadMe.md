# Quandenser: QUANtification by Distillation for ENhanced Signals with Error Regulation

Quandenser condenses quantification data from label-free mass spectrometry experiments.

## Installation

An installer for Ubuntu can be found on the Release page. Furthermore, Java >=1.6 has to be installed.

If you prefer to compile from source, or are running on a different operating system, [click here](#installation-from-source).

## Interface

Usage:

```
  quandenser -b <batch_file_list> -f <output_directory>
```

`<batch_file_list>` is a flat text file with the absolute path to each of the mzML files (one per line). The mzML files should be peak picked on both MS1 and MS2 level.

If everything executed correctly, you will find a file with the quantified feature groups at `<output_directory>/Quandenser.feature_groups.tsv` and one or several files with consensus spectra in `<output_directory>/consensus_spectra/MaRaCluster.consensus.part<x>.ms2`. Different output formats for the consensus spectra files can be specified using the `-o` option and the consensus spectra can subsequently be searched by any search engine.

## Example

An example run including downstream analysis with [Triqler](https://github.com/statisticalbiotechnology/triqler) can be found here: https://app.box.com/s/kp4219dc22l3gq27014nms8oco594c2i

The folder contains a ReadMe file with instructions on how to run the example.

## Installation from source

First, clone the repository using `git clone --recursive https://github.com/statisticalbiotechnology/quandenser.git`. Note the `--recursive` flag which is needed to pull in the submodules.

Quandenser depends on the Proteowizard and Boost libraries, which are automatically installed by the builder scripts. 

To install Quandenser, you can use the provided installation script `./quickbuild.sh` (Unix) or `./quickbuild.bat`/`./quickbuild64.bat` (Windows 32-bit and 64-bit respectively), which calls the appropriate build script for your platform located at `admin/builders/<platform>_build.<ext>`. Note that these scripts are currently largely untested except for Ubuntu. By default, it will install the executables in the `/usr/bin` folder (needs superuser rights). If you do not have superuser rights or want to install the executable somewhere else, modify the script accordingly by setting the `-DCMAKE_INSTALL_PREFIX` flag to the desired location.
