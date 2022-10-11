---
layout: default
title: Home
nav_order: 1
description: "SHAPEIT5 is a tool for haplotype phasing of high coverage sequencing data."
permalink: /
---
{: .warning }
Website under construction: content not available yet!

<img src="assets/images/branding/shapeit_logo.png" align="right" alt="Shapeit5" style="height:150px">

# SHAPEIT5
{: .fs-9 .fw-500 }

A tool for haplotype phasing of high-coverage sequencing data
{: .fs-5 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 .mx-auto }
[View source code on GitHub](https://github.com/odelaneau/shapeit5){: .btn .fs-5 .mb-4 .mb-md-0 }

---


<div class="code-example" markdown="1">
## About

SHAPEIT5 is a fast and accurate method for estimation of haplotypes for SNP array and and high-coverage datasets containing hundreads of thousands of samples.

## Citation

If you use SHAPEIT5 in your research work, please cite the following papers:

[Robin Hofmeister, Diogo Ribeiro, Simone Rubinacci, Olivier Delaneau. XXXXX. BiorXiv (2022)](https://www.nature.com/articles/s41588-020-00756-0)

## Features

- **Feature1**. Description 1.
- **Feature2**. Description 2.

## News

{: .new }
> **Version `1.0.0` is available!**
> See [the CHANGELOG](https://github.com/odelaneau/shapeit5/blob/main/versions/CHANGELOG.md) for a detailed breakdown.

## SHAPEIT5

Phase common
Version 1.0.0

Phases common variants using an improved version of the SHAPEIT4 algorithm

Phase rare
Version 1.0.0

Performs phasing of extremely rare variants down to singletons, using a scaffold of haplotypes built with Phase common

Documentation

---

## Getting started

{: .warning }
Website under construction: content not available yet!

### Dependencies


### Quick start: Use SHAPEIT5 on SNP array data


### Configure SHAPEIT5

- [See configuration options]({{ site.baseurl }}{% link docs/configuration.md %})

---

## About the project

SHAPEIT5 is developed by [Olivier Delaneau](https://odelaneau.github.io/lap-page).

### License

SHAPEIT5 is distributed with [MIT license](https://github.com/odelaneau/shapeit5/blob/main/LICENSE).

### Contributing

SHAPEIT5 is an open source project and we very much welcome new contributors. When contributing to our repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. Read more about becoming a contributor in [our GitHub repo](https://github.com/odelaneau/shapeit5#contributing).

#### Thank you to the contributors of SHAPEIT5!

<ul class="list-style-none">
{% for contributor in site.github.contributors %}
  <li class="d-inline-block mr-1">
     <a href="{{ contributor.html_url }}"><img src="{{ contributor.avatar_url }}" width="32" height="32" alt="{{ contributor.login }}"/></a>
  </li>
{% endfor %}
</ul>

We thank the [Just the Docs](https://github.com/just-the-docs/just-the-docs) developers, who made this awesome theme for Jekyll.
