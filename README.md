# BioTools

<!-- badges: start -->
<!-- badges: end -->

BioTools is a repository that has many tools and scripts which can be utilized to solve biological problems.

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
cd BioTools/tools
```

## Usage

### 1) Hclust Script
```
Rscript Hclust.R [-h] -i I -o O [-k K] [--clustering_distance CLUSTERING_DISTANCE] [--clustering_method CLUSTERING_METHOD]

mandatory arguments:
  -i I                                          Input file path
  -o O                                          Output folder path

optional arguments:
  -h, --help                                    show this help message and exit
  -k K                                          Number of clusters (optional; default value is 5)
  --clustering_distance CLUSTERING_DISTANCE     Clustering distance
  --clustering_method CLUSTERING_METHOD         Clustering method
```

#### Hclust Script Example

This is a basic example which shows you how to use Hclust script:

```
cd /path/to/BioTools/tools
Rscript Hclust.R -i ../data/Hclust.csv -o ../output/07_24_2020
```

### 2) PCA Script
```
Rscript PCA.R [-h] -i I -o O -s S

mandatory arguments:
  -i I        Input file path
  -o O        Output folder path
  -s S        Start column index

optional arguments:
  -h, --help  show this help message and exit
```

#### PCA Script Example

This is a basic example which shows you how to use PCA script:

```
cd /path/to/BioTools/tools
Rscript PCA.R -i ../data/03_22_2020_Arabidopsis_1001_BAA.csv -o ../output/07_25_2020 -s 3
```

### 3) Heatmap Script
```
Rscript Heatmap.R [-h] -i I -o O

mandatory arguments:
  -i I        Input file path
  -o O        Output folder path

optional arguments:
  -h, --help  show this help message and exit
```

#### Heatmap Script Example

This is a basic example which shows you how to use Heatmap script:

```
cd /path/to/BioTools/tools
Rscript Heatmap.R -i ../data/Hclust.csv -o ../output/07_26_2020
```

### 4) T-test and FDR Script
```
Rscript Ttest_and_FDR.R [-h] -i I -o O [--cores CORES] [--fdr_threshold FDR_THRESHOLD]

mandatory arguments:
  -i I                                  Input file path
  -o O                                  Output folder path

optional arguments:
  -h, --help                            show this help message and exit
  --cores CORES                         Number of clusters (optional; default value is 1)
  --fdr_threshold FDR_THRESHOLD         Threshold to filter FDR values
```

#### T-test and FDR Script Example

This is a basic example which shows you how to use T-test and FDR script:

```
cd /path/to/BioTools/tools
Rscript Ttest_and_FDR.R -i ../data/Maize_Drought_proteins_result.csv -o ../output/01_27_2021 --cores 10
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
