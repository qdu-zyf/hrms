# Hierarchial Meta-Storms: An R package for community ecologists —— popular ordination methods, ecological null models & diversity analysis

![Version](https://img.shields.io/badge/Version-1.01%20-brightgreen)
![Release date](https://img.shields.io/badge/Released%20date-Nov.%2019%2C%202020-brightgreen)



# Contents

- [Introduction](#introduction)

- [System Requirement and dependency](#system-requirement-and-dependency)

- [Installation guide](#installation-guide)

- [Usage](#usage)

- [Example dataset](#example-dataset)

- [Tools in this package](#tools-in-this-package)



# Introduction

Hierarchical Meta-Storms (HMS) comprehensively calculates the dissimilarities of microbiome functional profiles by considering multi-level metabolic pathway hierarchy. It contains two core components: *i)* a dissimilarity algorithm that comprehensively calculates the distances among microbiome functional profiles by considering their multi-level metabolic pathway hierarchy among functional gene families, and *ii)* a PCoA implementation optimized by multi-thread parallel computing to rapidly parse out the beta-diversity pattern for thousands of samples. It takes the microbiome functional profiles of KO and their relative abundance as input, and computes and outputs their pairwise distance matrix and then the principle coordinates of PCoA. 

In addition, the standalone package is also developed by C++ ([Github: hierarchical-meta-storms](https://github.com/qdu-bioinfo/hierarchical-meta-storms.git)) for direct installation and use under Linux and MAC operating systems.

# System Requirement and dependency

## Hardware Requirements

Hierarchical Meta-Storms only requires a standard computer with sufficient RAM to support the operations defined by a user. For typical users, this would be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

  RAM: 8+ GB  
  CPU: 4+ cores, 3.3+ GHz/core

## Software Requirements

OpenMP library is the C/C++ parallel computing library. Most Linux releases have OpenMP already been installed in the system. In Mac OS X, to install the compiler that supports OpenMP, we recommend using the Homebrew package manager:

```
brew install gcc
```

The package depends C++ (>= 4.7), R (>= 2.10) and links to Rcpp, RcppArmadillo and RcppEigen package, you can install these packages in R environment by:

```
install.packages(c("Rcpp", "RcppArmadillo", "RcppEigen"))
```

# Installation guide

## Automatic Installation by 

At present, Hierarchical Meta-Storms provides a fully automatic installer for easy installation.

**a. Download the package**

```
git clone https://github.com/qdu-bioinfo/hrms.git	
```

**b. Install the package**

```
R CMD INSTALL hrms
```

The package should take less than 1 minute to install on a computer with the specifications recommended above.

The example dataset could be found at “example” folder. Check the “example/Readme” for details about the demo run.

# Usage

**a.  Compute the distance matrix**

```
CompFunc(abd_matrix, rev=0, dist_type=0, is_sim=0)
```

The method returns the pairwise distance or similarity matix.

**b. Implement the PCoA**

```
GetPcoa(dist_matrix, k=3)
```

The method returns the coordinates matrix of PCoA. 

# Example dataset

Here we provide a demo dataset (Synthetic Dataset 1) with functional KO profiles of 30 artificial microbiomes in “data” folder. In this package, “dataset1.ko.abd” is the relative abundance on KOs-level, and “dataset1.meta” is the group information of the samples.

To run the demo, you can either:

```
cd hrms/example
sh Readme
```

or type the following command:

```
Rscript MS-comp-func.R -i dataset1.ko.abd -o dataset1.dist
Rscript MS-get-pc.R -d dataset1.dist -m dataset1.meta -o dataset1.pcoa.pdf
```

Then the output file “dataset1.dist” is the pairwise distance of the 30 samples as well as  "dataset1.pcoa.pdf.pc" is the coordinates results of PCoA based on the distance matrix and "dataset1.pcoa.pdf" is PCoA results printing to pdf format. 

This demo run should take less than 1 minute on a recommended computer.

# Methods in this package

**a. CompFunc**

It calculates the hierarchical meta-storms distance matrix among microbiome functional profiles. Run:

```
??CompFunc
```

in R environment for detailed parameters.

**b. GetPcoa**

It calculates the PCoA based the distance matrix. Run:

```
??GetPcoa
```

in R environment for detailed parameters.





