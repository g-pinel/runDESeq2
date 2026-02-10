library(tidyverse)
library(DESeq2)

runDESeq2 <- function(counts, metadata, target_var, covariates = NA, filter_var = NULL, filter_levels = NULL, outdir, gene_id_format = NULL, geneAnno = NULL, custom_label, pre_filtering = TRUE, sink_to_log = TRUE){
  
  # If needed, create results directory for the results
  # Create folder for results of differential expression analysis
  res_dirname <- custom_label
  res_dir <- file.path(outdir, res_dirname)
  
  # If results dir doesn't exist, create it
  if(!file.exists(res_dir)){
    message("Creating directory: ", res_dir)
    dir.create(res_dir, recursive = T)
  }else{
    message("Outputting files to directory: ", res_dir)
  }
  
  # Message output to file
  if(sink_to_log){
    message("Creating log file...")
    sink(file(paste0(res_dir, "/log.txt"), open = "wt"), 
         append=FALSE, 
         type = "message")
  }
  
  # Filter metadata
  if(!missing(filter_var)){
    metadata_filtered <- metadata %>% dplyr::filter(!!sym(filter_var) %in% filter_levels)
  }else{
    metadata_filtered <- metadata
  }
  
  message("Comparing variable: ", target_var, "\n")
  message("Variable levels: ", paste0(levels(metadata_filtered[,target_var]), collapse = ", "))
  if(all(!is.na(covariates))){message("Using covariates: ", paste(covariates, collapse = ", "))}
  
  # Get reference and treatment groups
  reference_group <- levels(metadata_filtered[,target_var])[1]
  treatment_group <- levels(metadata_filtered[,target_var])[2]
  
  # Get reference and treatment samples
  reference_samples <- metadata_filtered %>% dplyr::filter(!!sym(target_var) == reference_group) %>% rownames()
  treatment_samples <- metadata_filtered %>% dplyr::filter(!!sym(target_var) == treatment_group) %>% rownames()
  message("Using reference samples: ", paste0(reference_samples, collapse = ", "))
  message("Using treatment samples: ", paste0(treatment_samples, collapse = ", "))
  
  # Get filtered metadata for involved samples
  metadata_filtered <- metadata_filtered[rownames(metadata_filtered) %in%
                                           c(reference_samples, treatment_samples), ]
  
  # Save final used metadata
  message("Saving final metadata...")
  write.csv(metadata_filtered, paste0(res_dir, "/Metadata_", custom_label, ".csv"))
  
  # Get filtered counts
  counts_filtered <- counts[,rownames(metadata_filtered)]
  
  # Check needed for DESeq2
  if(all(colnames(counts_filtered) == rownames(metadata_filtered))){
    message("Sample names in count matrix and metadata match, continuing...")
  }else{
    stop("Mismatch between sample names in count matrix and metadata")
  }
  
  # Make formula
  if(all(is.na(covariates))){
    formula <- as.formula(paste("~", target_var))
  } else{
    formula <- as.formula(paste("~", paste(
      paste(covariates, collapse = "+"), 
      target_var,  
      sep = "+")))
  }
  message("Using formula: ", formula)
  
  # Construct DESeqDataSet
  message("Creating DESeqDataSet object...")
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts_filtered),
                                colData = metadata_filtered,
                                design = formula)
  
  # If indicated, run pre-filtering
  if(pre_filtering){
    message("Performing pre-filtering for genes with all samples <10 counts")
    keep <- rowSums(counts(dds) >= 10) >= min(c(length(reference_samples), length(treatment_samples))) # ncol(dds) equals number of samples in smallest group
    dds <- dds[keep,]
  }else{
    message("Parameter 'pre_filtering' set to FALSE, skipping pre-filtering")
  }
  
  # Run DE
  message("Running DE...")
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Format results
  res_df <- as.data.frame(res) %>% 
    rownames_to_column("gene_id") %>%
    dplyr::arrange(padj) %>% 
    dplyr::mutate(across(c(baseMean, log2FoldChange, lfcSE), ~round(.x, 2))) %>% 
    dplyr::mutate(across(c(pvalue, padj), ~format(.x, digits = 3))) %>% 
    dplyr::mutate(custom_label = custom_label) %>% 
    dplyr::mutate(across(c(pvalue, padj), as.numeric))
  
  if(!missing(geneAnno) & !missing(gene_id_format)){
    res_df <- res_df %>% 
      dplyr::left_join(geneAnno, by = c("gene_id" = gene_id_format))
  }else if(!missing(geneAnno) & missing(gene_id_format)){
    message("Please, provide a gene id format (gene_id_format parameter) matching a column in geneAnno to add gene-level annotations")
  }
  
  # Generate normalised counts
  message("Generating normalised counts file...")
  res_normCounts_df <- res_df %>%
    dplyr::left_join(counts(dds, normalized = TRUE) %>% 
                       as.data.frame() %>% 
                       rename_with(~str_c("Norm_", .), everything()) %>% 
                       rownames_to_column("gene_id"),
                     by = "gene_id") %>% 
    dplyr::select(-c("baseMean", "stat", "chromosome_name", "lfcSE", "pvalue"))
  
  # Stop output to log file
  if(sink_to_log){closeAllConnections()}
  
  write.csv(res_df, paste0(res_dir, "/DEG_", custom_label, ".csv"), row.names = F)
  write.csv(res_normCounts_df, paste0(res_dir, "/NormCounts_", custom_label, ".csv"), row.names = F)
  
  return(list(dds, res_df, res_normCounts_df))
}


batch_runDESeq2 <- function(input_comparisons, counts, metadata, gene_id_format = NULL, geneAnno = NULL, pre_filtering = TRUE, sink_to_log = TRUE){
  
  res <- list()
  
  for(i in 1:nrow(input_comparisons)){
    
    target_var <- input_comparisons[i, "target_var"] %>% as.character()
    covariates <- input_comparisons[i, "covariates"] %>% as.character()
    filter_var <- input_comparisons[i, "filter_var"] %>% as.character()
    filter_levels <- input_comparisons[i, "filter_levels"] %>% as.character()
    filter_levels <- str_split(filter_levels, pattern = ",") %>% unlist() %>% trimws()
    outdir <- input_comparisons[i, "outdir"] %>% as.character()
    custom_label <- input_comparisons[i, "custom_label"] %>% as.character()
    
    res[[i]] <- runDESeq2(counts = counts, 
                          metadata = metadata,
                          target_var = target_var, 
                          covariates = covariates, 
                          filter_var = filter_var, 
                          filter_levels = filter_levels, 
                          outdir = outdir,
                          gene_id_format = gene_id_format, 
                          geneAnno = geneAnno, 
                          custom_label = custom_label, 
                          pre_filtering = pre_filtering,
                          sink_to_log = sink_to_log)
  }
}
