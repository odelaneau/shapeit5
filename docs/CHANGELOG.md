---
title: CHANGELOG
layout: default
---

# CHANGELOG

All notable changes to this project are documented in this file.

{: .highlight }
Version 1.0.0 is now available [here] (https://github.com/odelaneau/shapeit5/releases) as docker image (`shapeit5_v1.0.0.tar.gz`) or static binaries. 

## v1.0.1

- Second release. Version used to generate the haplotypes of the 200k WGS samples of the UK Biobank.
- Support for pedigree file has been added to phase_common and phase_rare. Kids are scaffolded when possible according to Mendel logic and parental genomes.
- Support for new file formats has been added (Binaary and Sparse XCF files).
- IBD2 detection algorithm has been improved in phase_common.
- Multiple minor changes in default parameters.

## v1.0.0

- First release. Version used for the preprint
