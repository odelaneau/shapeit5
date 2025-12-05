---
layout: default
title: Home
nav_order: 1
description: "SHAPEIT5 is a tool for haplotype phasing of high coverage sequencing data."
permalink: /

---

![alt text](https://github.com/odelaneau/shapeit5/blob/main/docs/assets/images/branding/shapeit5_logo.png?raw=true)

<!---
# SHAPEIT5
{: .fs-9 .fw-500 }
-->

<!---
**S**egmented **HAP**lotype **E**stimation and **I**mputation **T**ools version **5**
{: .fs-5 }
-->

---

## About

SHAPEIT5 is a software package to estimate haplotypes in large genotype datasets (WGS and SNP array). 

## News

{: .new }
> **Version `5.1.0` is now available!**
> See [the CHANGELOG](https://github.com/odelaneau/shapeit5/blob/main/docs/CHANGELOG.md) for details.


## Citation

If you use SHAPEIT5 in your research work, please cite the following paper:

Hofmeister RJ, Ribeiro DM, Rubinacci S., Delaneau O. [Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. <br>Nature Genetics (2023) doi: https://doi.org/10.1038/s41588-023-01415-w ](https://www.nature.com/articles/s41588-023-01415-w)

---

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 .mx-auto }
[View source code on GitHub](https://github.com/odelaneau/shapeit5){: .btn .fs-5 .mb-4 .mb-md-0 }


## Description

SHAPEIT5 is composed of the following tools:

- **phase_common**. Tool to phase common sites, typically SNP array data, or the first step of WES/WGS data phasing pipeline. This tool replaces SHAPEIT4.
- **ligate**. Ligate multiple phased BCF/VCF files into a single whole chromosome file. Typically run to ligate multiple chunks of phased common variants.
- **phase_rare**. Tool to phase rare variants onto a scaffold of common variants (output of phase_common / ligate).
- **switch**. Program to compute switch error rate and genotyping error rate given simulated or trio data.
- **xcftools**. Program to convert between the various file formats used by shapeit5 [BCF from/to XCF].

[phase_common]({{site.baseurl}}{% link docs/documentation/phase_common.md %}){: .btn .btn-blue }
[ligate]({{site.baseurl}}{% link docs/documentation/ligate.md %}){: .btn .btn-blue }
[phase_rare]({{site.baseurl}}{% link docs/documentation/phase_rare.md %}){: .btn .btn-blue }
[switch]({{site.baseurl}}{% link docs/documentation/switch.md %}){: .btn .btn-blue  }
[xcftools]({{site.baseurl}}{% link docs/documentation/switch.md %}){: .btn .btn-blue  }


---

## Getting started

- [See documentation]({{site.baseurl}}{% link docs/documentation/documentation.md %})

---

## About the project

SHAPEIT5 is developed by [Olivier Delaneau's group](https://odelaneau.github.io/lap-page).

### License

SHAPEIT5 is distributed with [MIT license](https://github.com/odelaneau/shapeit5/blob/main/LICENSE).

### Organisations

<div class="d-flex justify-content-around">
  <div class="p-5"><a href="https://www.unil.ch/index.html"><img src="assets/images/lausanne_logo.jpg" align="right" alt="unil" style="height:50px"></a></div>
  <div class="p-5"><a href="https://www.sib.swiss/"><img src="assets/images/sib_logo.jpg" align="right" alt="sib" style="height:50px"></a></div>
  <div class="p-5"><a href="https://www.snf.ch/en/Pages/default.aspx"><img src="assets/images/snf.gif" align="right" alt="snf" style="height:50px"></a></div>
</div>

### Contributing

SHAPEIT5 is an open source project and we very much welcome new contributors. When contributing to our repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change.
#### Thank you to the contributors of SHAPEIT5!

<ul class="list-style-none">
{% for contributor in site.github.contributors %}
  <li class="d-inline-block mr-1">
     <a href="{{ contributor.html_url }}"><img src="{{ contributor.avatar_url }}" width="32" height="32" alt="{{ contributor.login }}"/></a>
  </li>
{% endfor %}
</ul>

We thank the [Just the Docs](https://github.com/just-the-docs/just-the-docs) developers, who made this awesome theme for Jekyll.
