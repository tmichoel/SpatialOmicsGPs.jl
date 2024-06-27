# SpatialOmicsGPs

A crude (for now) implementation of various Gaussian process models for [spatial transcriptomics][4] data. Currently planned to contain implementations of:

**SpatialDE** [[Svensson et al.]][1]:

- [x] Variance parameter estimation in SpatialDE model
- [x] Statistical significance of spatial vs. non-spatial covariance
- [ ] Expression histology

**Spatial Variance Component Analysis (SVCA)** [[Arnol et al.]][2]:

- [ ] Variance parameter estimation in SVCE model

The package is developed in parallel with another package, [FaSTLMMlight](https://github.com/tmichoel/FaSTLMMlight.jl) that provides the variance parameter estimation methods used here.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tmichoel.github.io/SpatialOmicsGPs.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tmichoel.github.io/SpatialOmicsGPs.jl/dev/)
[![Build Status](https://github.com/tmichoel/SpatialOmicsGPs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tmichoel/SpatialOmicsGPs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tmichoel/SpatialOmicsGPs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tmichoel/SpatialOmicsGPs.jl)


[1]: https://doi.org/10.1038%2Fnmeth.4636
[2]: https://doi.org/10.1016/j.celrep.2019.08.077
[3]: https://arxiv.org/abs/2210.10210
[4]: https://en.wikipedia.org/wiki/Spatial_transcriptomics
[5]: https://doi.org/10.1038/nmeth.1681