# Segmented HAPlotype Estimation and Imputation Tools version 5 (SHAPEIT5)

SHAPEIT5 estimates haplotypes in large datasets, with a special focus on rare variants.

- **docker**: all script needed to build a docker file for shapeit5
- **docs**: documentation in html
- **ligate**: ligate multiple phased BF/VCF files into a chromosome length file
- **maps**: genetics maps in b37 and b38
- **phase_common**: phase common sites, typically SNP array data
- **phase_rare**: phase rare variants onto a scaffold of common variants
- **static_bins**: static binaries of all execeutables
- **switch**: compute switch error rate and genotyping error rate given simulated or trio data
- **tasks**: scripts used to phase large datasets, good base to start pipelining
- **test**: simulated data for first-step testing of the method
- **versions**: versioning

Delaneau O., et al. Accurate, scalable and integrative haplotype estimation. Nature Communications volume 10, Article number: 5436 (2019). 
https://www.nature.com/articles/s41467-019-13225-y

## Documentation

https://odelaneau.github.io/shapeit5/

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
