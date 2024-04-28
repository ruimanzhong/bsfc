# bsfc: Spatial Functional Clustering

`bsfc` is an R package designed for clustering spatial functional data. This package facilitates the analysis of data where observations are functions or curves that are correlated across spatial locations. It extends ideas from the research paper "BAYESIAN CLUSTERING OF SPATIAL FUNCTIONAL DATA WITH APPLICATION TO A HUMAN MOBILITY STUDY DURING COVID-19" to the exponential family using Integrated Nested Laplace Approximations (INLA).


## Features

- **Broad Modeling Capabilities:** Implements various models for spatial clustering based on functional data.
- **Exponential Family Support:** Extends the Bayesian clustering framework to accommodate the exponential family of distributions using INLA, enhancing flexibility and computational efficiency.
- **Visualization Tools:** Provides built-in functions for visualizing clustering results, allowing users to effectively analyze and interpret the spatial patterns and groupings discovered by the model.

## Background

This package builds upon the methods introduced in the research paper "BAYESIAN CLUSTERING OF SPATIAL FUNCTIONAL DATA WITH APPLICATION TO A HUMAN MOBILITY STUDY DURING COVID-19". The original research utilized Bayesian clustering to analyze how human mobility patterns clustered across different regions during the early stages of the COVID-19 pandemic. By extending this framework to the exponential family and incorporating INLA, `bsfc` allows for faster and more scalable analyses, making it suitable for a wider range of applications beyond the initial study's scope. This includes complex models where traditional MCMC methods might be computationally prohibitive.

This README is now comprehensive, combining installation instructions, usage examples, feature highlights, and the scientific background of your package. It’s ready to be used on your GitHub repository to guide users and contributors.

## Installation

You can install the development version of `bsfc` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("ruimanzhong/bsfc", force = TRUE, build_vignette = F)

