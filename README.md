# Segmented HAPlotype Estimation and Imputation Tools version 5 (SHAPEIT5)

SHAPEIT5 estimates haplotypes in large datasets, with a special focus on rare variants.

- **docker**: all script needed to build a docker file comprising all binaries
- **docs**: documentation in html
- **ligate**: ligate multiple phased BCF/VCF files into a chromosome length file
- **resources**: genetics maps in b37 and b38 amd coordinates for WGS chunks
- **phase_common**: phase common sites, typically SNP array data
- **phase_rare**: phase rare variants onto a scaffold of common variants
- **static_bins**: static binaries of all executables
- **switch**: compute switch error rate and genotyping error rate given simulated or trio data
- **tasks**: scripts used to phase large datasets, good base to start pipelining
- **test**: simulated data for first-step testing of the method
- **versions**: versioning
- **xcftools**: tools to handle XCF files [experimental]

## Citation
If you use SHAPEIT5 in your research work, please cite the following paper:

Hofmeister RJ, Ribeiro DM, Rubinacci S., Delaneau O. [Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank](https://www.biorxiv.org/content/10.1101/2022.10.19.512867v1). BiorXiv (2022)

## Documentation

Documentation, installation instructions and tutorials can be found at:

[https://odelaneau.github.io/shapeit5/](https://odelaneau.github.io/shapeit5/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details