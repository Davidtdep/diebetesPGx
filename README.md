# diabetesPGx

This repository contains the complete R code and datasets used in the article:  
**_Pharmacogenomics in Diabetes: What We Know, What We Overlook, and Who It Affects in Colombia_** (under submission).

## Purpose

This study explores the influence of genetic ancestry on pharmacogenomic profiles related to diabetes in Colombian populations. We assessed allele frequency differences between African- and European-ancestry groups, performed correlation analyses, and conducted a bibliometric evaluation of studies indexed in PharmGKB.

This repository includes both the full dataset and analysis scripts, enabling full reproducibility of all results.

## Repository Structure

- `main_analysis.R`: Main script with numbered and annotated sections for reproducibility.
- `data/`: Includes:
  - `data.tsv`: Full set of annotated pharmacogenomic variants from PharmGKB.
  - `data_significant.csv`: Filtered dataset used for statistical and graphical analyses.

## Required R Packages

The following R packages are required:

- `dplyr`
- `ggplot2`
- `ggExtra`
- `ggsci`
- `reshape2`
- `circlize`
- `openxlsx`
- `rentrez`
- `httr`
- `jsonlite`
- `xml2`
- `rcrossref`
- `fuzzyjoin`
- `stringr`

Install them using:

```r
pkgs <- c("dplyr", "ggplot2", "ggExtra", "ggsci", "reshape2", "circlize", 
          "openxlsx", "rentrez", "httr", "jsonlite", "xml2", "rcrossref", 
          "fuzzyjoin", "stringr")
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

## Analyses Included

The `main_analysis.R` script performs the following:

- **Data cleaning and filtering**: Selection of statistically significant SNPs with single-allele annotations.
- **Ancestry-specific allele frequency comparison**: Among African (PLQ, CHG) and European (ATQCES, ATQPGC, CLM) Colombian populations.
- **Correlation matrix and heatmap**: Pearson correlations across populations.
- **Scatterplots**: Visual comparison of allele frequencies between ancestry groups.
- **Chord diagram**: Distribution of pharmacogenomic study cases and controls by ancestry.
- **Descriptive plots**: Top drugs and genes associated with efficacy/toxicity outcomes.
- **Bibliometric metadata extraction**: Using Entrez, Unpaywall, Crossref, and SCImago APIs.
- **Descriptive bibliometric analysis**: Number of publications by WHO region, income group, journal quartile, and open access status.

## How to Use

1. Clone or download this repository.
2. Open `main_analysis.R` in RStudio.
3. Ensure the files in the `data/` folder are correctly referenced.
4. Install required R packages (see above).
5. Run the script from beginning to end or explore specific numbered sections.

## Reproducibility

> ✅ **This repository is fully reproducible**: All datasets and scripts are included, and no external data download is required.

## Citation

If you use or adapt this code, please cite the associated article (once published) and acknowledge this repository. A DOI will be added upon publication.

## License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute the code with attribution.
