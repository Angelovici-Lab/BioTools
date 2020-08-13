# BioTools

<!-- badges: start -->
<!-- badges: end -->

BioTools is a repository that has many scripts which can be utilized to solve biological problems.

## Installation

You can install the BioTools from [Github](https://github.com/Angelovici-Lab/BioTools) with:

``` r
# Run this inside R environment
install.packages("devtools", dependencies = TRUE)
devtools::install_github("Angelovici-Lab/BioTools")
```
```
# Run this in your terminal
git clone https://github.com/Angelovici-Lab/BioTools.git
cd BioTools/scripts
```

## Usage

### 1) Hclust Script
```
Rscript Hclust.R [-h] -i I -o O [-k K]

optional arguments:
  -h, --help  show this help message and exit
  -i I        Input file path
  -o O        Output folder path
  -k K        Number of clusters (optional; default value is 5)
```

#### Hclust Script Example

This is a basic example which shows you how to use Hclust script:

```
cd /path/to/BioTools/scripts
Rscript Hclust.R -i ../data/Hclust.csv -o ../output/07_30_2020
```

### 2) PCA Script
```
Rscript PCA.R [-h] -i I -o O [-s S]

optional arguments:
  -h, --help  show this help message and exit
  -i I        Input file path
  -o O        Output folder path
  -s S        Start column (optional; default value is 1)
```

#### PCA Script Example

This is a basic example which shows you how to use PCA script:

```
cd /path/to/BioTools/scripts
Rscript PCA.R -i 03_22_2020_Arabidopsis_1001_BAA.csv -o ../output/08_06_2020
```

## Package Update

To upgrade BioTools to the latest version, please remove the package and re-install the latest BioTools package:

``` r
# Run this inside R environment
remove.packages("BioTools")
devtools::install_github("Angelovici-Lab/BioTools")
```

``` r
# Run this in your terminal
cd /path/to/BioTools
git pull
```
