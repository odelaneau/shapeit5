---
layout: default
title: Compile SHAPEIT5
parent: Build from source
grand_parent: Installation
permalink: /docs/installation/build_from_source/compile_shapeit5
---
# Compile SHAPEIT5
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Compile SHAPEIT5
Download the last version of the SHAPEIT5 code using:
<div class="code-example" markdown="1">
```bash
git clone --recurse-submodules https://github.com/odelaneau/shapeit5.git
```
</div>

Adding `--recurse-submodules` allows cloning automatically the xcftools submodule required by shapeit5. 

Navigate to the downloaded folder using `cd shapeit5`.

You'll find there a folder containing all the software packages are other utility folders:

- **docker**: all scripts needed to build a docker file comprising all binaries
- **docs**: documentation in html
- **ligate**: ligate multiple phased BCF/VCF files into a chromosome length file
- **phase_common**: phase common sites, typically SNP array data
- **phase_rare**: phase rare variants onto a scaffold of common variants
- **resources**: genetics maps [b37/b38] and coordinates for [5/20] cM chunks  
- **static_bins**: static binaries of all executables
- **switch**: compute switch error rate and genotyping error rate given simulated or trio data
- **tasks**: scripts used to phase large datasets, good base to start pipelining
- **test**: simulated data for first-step testing of the method
- **versions**: versioning
- **xcftools**: tool to handle XCF file format


Each software in the suite contains the same folder structure:

- `bin`: folder for the compiled binary.
- `obj`: folder with all binary objects.
- `src`: folder with source code.
- `makefile`: Makefile to compile the program.

In order to compile a specific tool, for example _phase\_rare_, go in directory of the software (cd `phase_rare`) and edit the specific makefile at lines so that the following variables are correctly set up (look at the paths already there for an example):

- `HTSSRC`: path to the root of the HTSlib library, the prefix for HTSLIB_INC and HTSLIB_LIB paths.
- `HTSLIB_INC`: path to the HTSlib header files.
- `HTSLIB_LIB`: path to the static HTSlib library (file `libhts.a`).
- `BOOST_INC`: path to the BOOST header files (usually `/usr/include`). 
- `BOOST_LIB_IO`: path to the static BOOST iostreams library (file `libboost_iostreams.a`). 
- `BOOST_LIB_PO`: path to the static BOOST program_options library (file `libboost_program_options.a`). 

If installed at the system level, static libraries (*.a files) can be located with this command:

<div class="code-example" markdown="1">
```bash
locate libboost_program_options.a libboost_iostreams.a libhts.a
```
</div>

Once all paths correctly set up, proceed with the compilation using `make`. The binary can be found in the `bin/` folder of each tool and will have a name similar to `SHAPEIT5_phase_common`. You will need to copy the modified makefile in each tool (folder) of SHAPEIT5.

