removeNeighborSnps <- function(
                        snps, 
                        maxdist = 50000, 
                        biomart.config = biomartConfigs$hsapiens,
                        use.buffer = FALSE
                      ) {
  
  if(nrow(snps) < 1)
    return(snps)
  
  if(!is.null(snps$BP) & !is.null(snps$CHR)) {
    # dont need further information
  } else if(!is.null(snps$SNP)) {
    cat("Fetching SNP information from biomart...\n")
    snps.bm <- bm.snps(filter.val = snps$SNP, config = biomart.config, use.buffer = use.buffer)
    names(snps.bm)[ names(snps.bm) == "chrname" ] <- "CHR"
    names(snps.bm)[ names(snps.bm) == "chrom_start" ] <- "BP"
    # merge will sort the data frame, SNPs that do not match with ensembl are removed
    snps <- merge( snps, snps.bm, by.x = "SNP", by.y = "refsnp_id", suffixes = c(".original","") )
  } else {
    stop("Parameter 'snps' has to contain either column 'SNP' or 'BP' and 'CHR'")
  }
  
  snps <- snps[!is.na(snps$BP) & !is.na(snps$CHR), ]
  
  if(is.null(snps$P)) {
    # we require SNPs not to have identical BPs and CHR (would both be removed)
    snps <- snps[!duplicated(snps[, c("CHR", "BP")]), ]
    # for each SNP, remove when there are SNPs within maxdist AND having larger index
    return(snps[
              sapply(
                1:nrow(snps),
                function(idx) {
                  # snps in range
                  cand.idx <- which(
                      as.vector(snps$CHR) == as.vector(snps[idx, "CHR"]) & 
                      abs(as.numeric(snps$BP) - as.numeric(snps[idx, "BP"])) <= maxdist
                    )
                  if(length(cand.idx) < 1) 
                    return(TRUE)
                  if(any(max(cand.idx) > which(as.vector(snps$CHR) == as.vector(snps[idx, "CHR"]) & as.numeric(snps$BP) == as.numeric(snps[idx, "BP"]))))
                    return(FALSE)
                  else
                    return(TRUE)
                }
                )
              , ])
  } else {
    # we require SNPs not to have identical BPs and CHR (would both be removed)
    # remove duplicates (when P column exists, remove the larger P!)
    snps <- snps[order(snps$P), ]
    snps <- snps[!duplicated(snps[, c("CHR", "BP")]), ]
    # go through p-values descendingly, keep only those that are not in range with a higher p-value
    return(snps[
              sapply(
                1:nrow(snps),
                function(idx) {
                  # snps in range
                  cand.idx <- which(
                      as.vector(snps$CHR) == as.vector(snps[idx, "CHR"]) & 
                      abs(as.numeric(snps$BP) - as.numeric(snps[idx, "BP"])) <= maxdist
                    )
                  if(length(cand.idx) < 1) 
                    return(TRUE)
                  # among them, are there any with lower P (or equal P and lower BP)
                  if(
                      any(as.numeric(snps[cand.idx, "P"]) < as.numeric(snps[idx, "P"])) | 
                      any(as.numeric(snps[cand.idx, "P"]) == as.numeric(snps[idx, "P"]) & as.numeric(snps[cand.idx, "BP"]) < as.numeric(snps[idx, "BP"]))
                    )
                    return(FALSE)
                  else
                    return(TRUE)
                }
                )
              , ])
  }
}
