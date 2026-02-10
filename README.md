# runDESeq2

### Easily run differential expression comparisons in batch with DESeq2

Bulk RNA-seq experiments often include multiple experimental groups and involve performing multiple pairwise comparisons between different subsets of samples.

The functions in this repository allow you to perfom one or multiple differential expression analyses using the DESeq2 package.
<br/><br/>

### Dependencies

These function require the tidyverse and DESeq2 R packages, which can be installed in R with:

```
install.packages("tidyverse")
```
and
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```
<br/>

## batch_runDESeq2 function (run multiple comparisons in batch)

This is a wrapper function that passes all the arguments to the *runDESeq2* function looping over all the desired comparisons.
<br/><br/>

**Comparisons sheet**: The parameters corresponding to each comparison are provided in an Excel file (XLSX) provided to the function with the **input_comparisons** parameter. Should contain the following columns:
- **custom_label**: the label for the comparison. It will be included in the results directory and files names.
- **target_var**: the name of the metadata variable of interest for the comparison. Must be a factor with two levels in the metadata. As per DESeq2's vignette, the first level corresponds to the reference and the second level to the treatment.
- **covariates** (optional): the name of the metadata variables that should be used as covariates. If more than one, separate them by commas.
- **filter_var** (optional): if before performing the comparisons you want to filter your data using another variable, indicate the metadata column name of that variable here.
- **filter_levels** (if using *filter_var*): the level of *filter_var* that should be subset before performing the analysis.
- **outdir**: the directory where all result files will be output.

This is how the XLSX file would look like to perform the comparison 
'Male gonads versus female gonads using chromosomes as a covariate and filtering for high-fat diet samples':

![Example XLSX file containing the comparison parameters. The custom_label column contains the string "gonad_HF", the target_var column the "gonads" string, the covariates column the "chromosomes" string, the filter_var column the "diet" string, the filter_levels column the "HF" string, short for high fat, and the outdir column the "results/DE/" string](/assets/example_xlsx.png "Example runDESeq2 input XLSX file")

To run multiple comparisons, simply provide multiple lines in this file with the corresponding paramters.
<br/><br/>

The rest of parameters of this function are:
- **counts**: a matrix with gene IDs as row names and samples as column names containing raw or gene-length-adjusted counts.
- **metadata**: a data frame with samples as row names and sample metadata variables as column names containing sample-level metadata.
- **geneAnno** (optional): a data frame containing gene-level metadata that will be added to the differential expression results.
- **gene_id_format** (optional): the column in *geneAnno* containing the gene IDs that match those in *counts* used to join both data frames.
- **pre_filtering**: if set to TRUE (default), genes with less than 10 counts in less than as many samples as the smallest experimental group involved in the comparison has will be removed (e.g. if the smallest group has 4 samples, genes with less than 4 samples with >10 samples will be excluded). Set to FALSE to avoid pre-filtering.
- **sink_to_log**: if set to TRUE (default), it will create a log file providing information about the differential expression analysis. Set to FALSE to show this information in the R console instead.
> Note: the samples in counts column names must match and be in the same order as those in the metadata row names.
<br/><br/>

## runDESeq2 function

This is the function that performs the differential expression analysis. Can be run separately from batch_runDESeq2 if you only wish to perform one comparison. It performs the following: 
- Creates a results directory/detects an existing one matching *custom_label* inside the output directory for the results from each comparisons.
- Creates a log file inside that directory if *sink_to_log* is set to TRUE (default).
- Filters the metadata and counts if *filter_var* and *filter_levels* were provided.
- Creates an XLSX file containing the final metadata used for the analysis.
- Verifies that all counts column names match all metadata row names.
- Creates the design formula in the format '~ covariates + target_var'
- Creates a DESeqDataSet object, performs pre-filtering if *pre_filtering* is set to TRUE (default), and performs count normalisation followed by differential expression analysis.
- Formats differential expression results (rounds the *baseMean*, *log2FoldChange* and *lfcSE* columns to two decimal places; formats the *pvalue* and *padj* coplumns to scientific notation; adds a column containing *custom_label*).
- Adds *geneAnno* to the results if provided, using *gene_id_format* as the column for joining.
- Generates XLSX files containing the differential expression results and normalised counts.
- Returns the DESeqDataSet object, the differential expression results and the normalised counts as a list.
<br/><br/>

### batch_runDESeq2 call example
```
de_results <- batch_runDESeq2(input_comparisons = comparisons_df,
  counts = counts_mtx,
  metadata = metadata_df,
  gene_id_format = "ensembl_gene_id",
  geneAnno = mm_gene_df)
```
