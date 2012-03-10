readPvalFiles <- function(pval.files, assert.positions.match = F, read.columns = c("CHR", "BP", "SNP", "P")) {

    cat("Reading pval files...\n")
    p.all <- NULL
    for( pval.file in pval.files )  {
        p.next <- read.table(pval.file, header=TRUE)
        if(!all(read.columns %in% colnames(p.next)))
          stop(paste("Invalid p-value file (GWAS result file) format, required columns ", paste(read.columns, collapse = "', '"), " were not found!\n", sep = "'"))
        p.next <- p.next[, read.columns]
        p.next <- na.omit(p.next)
        p.next$pheno = pval.file # pval dataset identifier, is name of pval file here
        p.all <- rbind(p.all,p.next)
    }
    
    if(assert.positions.match && "BP" %in% read.columns) {
      p.all.test <- merge(p.all[, c("SNP", "BP")], p.all[, c("SNP", "BP")], by.x = "SNP", by.y = "SNP")
      if(any(p.all.test$BP.x != p.all.test$BP.y))
        stop(paste(
          "Base positions do not match between pvalue files for SNP", 
          as.vector(p.all.test[which(p.all.test$BP.x != p.all.test$BP.y)[1], "SNP"]), 
          ". Use the bm.remapSnps function on both files first or do a manual remapping."
        ))
      rm(p.all.test)
    }

    return(p.all)
}

# returns data from a pedigree mapfile (linkage format, either 4 or 3 columns) without header. 
# Return value is a data frame with columns "CHR", "SNP", "BP"
readMapfile <- function(mapfile) {
  mapinfo <- read.table(mapfile, header = FALSE, row.names = NULL)
  if(ncol(mapinfo) == 4) mapinfo <- mapinfo[, c(1,2,4)]
  if(ncol(mapinfo) != 3) 
    stop("Mapfile format unknown")
  colnames(mapinfo) <- c("CHR", "SNP", "BP")
  return(mapinfo)
}
