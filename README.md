# Segmented HAPlotype Estimation and Imputation Tools version 5 (SHAPEIT5)

[![Build](https://github.com/odelaneau/shapeit5/actions/workflows/build.yml/badge.svg)](https://github.com/odelaneau/shapeit5/actions) [![Build_Pages](https://github.com/odelaneau/shapeit5/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/odelaneau/shapeit5/actions)

SHAPEIT5 estimates haplotypes in large datasets, with a special focus on rare variants.

- **docker**: all script needed to build a docker file comprising all binaries
- **docs**: documentation in html
- **ligate**: ligate multiple phased BCF/VCF files into a chromosome length file
- **resources**: genetics maps in b37 and b38 amd coordinates for WGS chunks
- **phase_common**: phase common sites, typically SNP array data
- **phase_rare**: phase rare variants onto a scaffold of common variants
- **static_bins**: static binaries of all executables
- **simulate**: simulate simple haplotype datasets
- **switch**: compute switch error rate and genotyping error rate given simulated or trio data
- **tasks**: scripts used to phase large datasets, good base to start pipelining
- **test**: simulated data for first-step testing of the method
- **versions**: versioning
- **xcftools**: tools to handle XCF files [experimental]

## Download

It is strongly advised to prioritize the use of the latest released version (see released section on the right hand side) instead of cloning the latest version of the source code, as doing so is at your own risk.

## Citation
If you use SHAPEIT5 in your research work, please cite the following paper:

Hofmeister RJ, Ribeiro DM, Rubinacci S., Delaneau O. [Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank](https://www.nature.com/articles/s41588-023-01415-w). Nature Genetics (2023)

## Documentation

Documentation, installation instructions and tutorials can be found at:

[https://odelaneau.github.io/shapeit5/](https://odelaneau.github.io/shapeit5/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Guard against large files (50MB) and Git LFS

This repository includes a pre-commit hook that blocks committing files larger than 50 MB by default. The hook lives in `.githooks/pre-commit` and is versioned with the code.

Enable it locally once per clone:

```
bash scripts/setup-git-hooks.sh
```

You can temporarily change the limit for a single commit (discouraged):

```
MAX_GIT_FILE_SIZE_MB=100 git commit -m "Allow up to 100MB"
```

For large binary assets (e.g., BCF/VCF), use Git LFS to avoid oversized commits and bloating clone sizes. The repository provides default patterns in `.gitattributes` for common genomics formats. To initialize LFS locally:

```
git lfs install
# patterns are already in .gitattributes; adjust as needed
```

If you see the pre-commit hook block a commit due to file size, consider whether the file should be ignored, generated, or tracked with Git LFS.
