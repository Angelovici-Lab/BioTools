# BioTools

<!-- badges: start -->
<!-- badges: end -->

BioTools is a repository that has many scripts to solve bio problems.

## Installation

You can install the BioTools from [Github](https://github.com/Angelovici-Lab/BioTools) with:

``` r
# Run this in your terminal
git clone https://github.com/Angelovici-Lab/BioTools.git
cd BioTools
cd scripts
```

## Usage

### 1) Hclust Script
``` r
Rscript Hclust.R [-h] -i I -o O [-k K]

optional arguments:
  -h, --help  show this help message and exit
  -i I        Input file path
  -o O        Output folder path
  -k K        Number of clusters (default value is 5)
```

#### Hclust Script Example

This is a basic example which shows you how to use Hclust script:

``` r
cd /path/to/BioTools/scripts
Rscript Hclust.R -i ../data/Hclust.csv -o ../output/07_30_2020
```
