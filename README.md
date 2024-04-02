# SpatialTranscriptomicsGPs

A crude (for now) implementation of various Gaussian process models for [spatial transcriptomics][4] data. Currently planned to contain implementations of:

**FaST-LMM** [[Lippert et al.]][5]:

- [x] LMMs/GPs with full-rank kernel matrix
- [ ] LMMs/GPs with low-rank kernel matrix factorization

**SpatialDE** [[Svensson et al.]][1]:

- [ ] Variance parameter estimation in SpatialDE model
- [ ] Statistical significance of spatial vs. non-spatial covariance

**Spatial Variance Component Analysis (SVCA)** [[Arnol et al.]][2]:

- [ ] Variance parameter estimation in SVCE model

I'm also hoping to implement the method of [Greengard et al.][3] to see if it leads to further speed-up.

[![Build Status](https://github.com/tmichoel/SpatialTranscriptomicsGPs.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/tmichoel/SpatialTranscriptomicsGPs.jl/actions/workflows/CI.yml?query=branch%3Amaster)

[1]: https://doi.org/10.1038%2Fnmeth.4636
[2]: https://doi.org/10.1016/j.celrep.2019.08.077
[3]: https://arxiv.org/abs/2210.10210
[4]: https://en.wikipedia.org/wiki/Spatial_transcriptomics
[5]: https://doi.org/10.1038/nmeth.1681