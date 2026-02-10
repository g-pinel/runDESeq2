# runDESeq2

## Easily run differential expression comparisons in batch with DESeq2

Bulk RNA-seq experiments often include multiple experimental groups and involve performing multiple pairwise comparisons between different subsets of samples.

The functions in this repository allow you to perfom one or multiple differential expression analyses using the DESeq2 package.

### batch_runDESeq2 function (run multiple comparisons in batch)

The parameters corresponding to each comparison are provided in an Excel file (XLSX) with the following columns:
- **custom_label**: the label for the comparison. It will be included in the results directory and files names.
- **target_var**: the name of the metadata variable of interest for the comparison. Must be a factor with two levels in the metadata. As per DESeq2's vignette, the first level corresponds to the reference and the second level to the treatment.
- **covariates** (optional): the name of the metadata variables that should be used as covariates. If more than one, separate them by commas.
- **filter_var** (optional): if before performing the comparisons you want to filter your data using another variable, indicate the metadata column name of that variable here.
- **filter_levels** (if using *filter_var*): the level of *filter_var* that should be subset before performing the analysis.
- **outdir**: the directory where all result files will be output.

This is how the XLSX file would look like to perform the comparison 
'I want to compare male gonads versus female gonads using chromosomes as a covariate and filtering for high-fat diet samples':

![Example XLSX file containing the comparison parameters. The custom_label column contains the string "gonad_HF", the target_var column the "gonads" string, the covariates column the "chromosomes" string, the filter_var column the "diet" string, the filter_levels column the "HF" string, short for high fat, and the outdir column the "results/DE/" string](/assets/example_xlsx.png "Example runDESeq2 input XLSX file")
