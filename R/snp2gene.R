snp2gene <- function(
              snps, 
              by.ld = FALSE, 
              biomart.config = biomartConfigs$hsapiens, 
              use.buffer = FALSE,
              by.genename = TRUE,
              level = -1, 
              print.format = FALSE, 
              run.parallel = FALSE, 
              ... 
            ) {

  if(!all(names(list(...)) %in% c("map", "ped")))
    stop("Unknown parameter! Check argument list.\n")
  
  if(missing(snps)) 
    stop("Function parameter 'snps' has not been specified\n")

  if(!is.data.frame(snps) || nrow(snps) < 1) 
    stop("Function parameter snps is not a data frame or empty\n")

  if(is.null(snps$SNP))
    stop("Function parameter 'snps' does not contain 'SNP' column\n")

  snps.forbidden.columns <- c("CHR.original", "BP.original", "BP.mapped", "CHR.mapped", "chrname", "chromo", "direction", "genename", "geneid", "start", "end")
  if(any(colnames(snps) %in% snps.forbidden.columns))
    stop(paste("Function parameter 'snps' may not contain one of the columns", paste(snps.forbidden.columns, collapse = "', '"), "\n", sep = " '"))

  if( is.null(by.ld) || (by.ld != F && by.ld != T) )
    stop("Function parameter by.ld can only take the values TRUE or FALSE\n")
  
  if(is.null(by.genename))
    by.genename <- TRUE
  
  if(is.null(level))
    level <- -1
  if(!is.numeric(level))
    stop("Function parameter 'level' has to be numeric.")
  level <- round(level)

  if(is.null(use.buffer))
    use.buffer <- FALSE
  
  # lazy initialization - assigned on use (thus buffer data can be used without connection)
  delayedAssign("snpds", bm.init.snps(biomart.config))
  delayedAssign("geneds", bm.init.genes(biomart.config))

  cat("Processing input data...\n")
  snps.bm <- bm.snps(biomart.config, snpds, snps$SNP, use.buffer)
  names(snps.bm) <- c("refsnp_id", "CHR", "BP")   
  # merge will sort the data frame, SNPs that do not match with ensembl are removed
  snps <- merge(snps, snps.bm, by.x = "SNP", by.y = "refsnp_id", suffixes = c(".original",""))
  snps <- snps[!is.null(snps$CHR) & !is.null(snps$BP) & !is.na(snps$CHR) & !is.na(snps$BP), ]
  # there may be rare cases of SNPs with multiple positions in Biomart - keep the first
  snps <- snps[!duplicated(snps$SNP), ] 
  genes <- bm.genes(biomart.config, geneds, use.buffer)
  if(by.genename) {
    genes <- genes[!(genes$genename == ""), ]
    # ldgenes and closestGenes can cope with multi-chromosmes per genename, but we only want each genename once per chromo
    genes <- genes[!duplicated(genes[, c("genename", "chrname")]), ]
  }
  
  if(nrow(snps) < 1)
    stop("No SNPs left after biomart mapping. Aborting.\n")

  cat("Calculating closest genes...\n")
  if(by.ld) {

    ld <- ldGenes(snps, genes, run.parallel, ...)
    ld <- ld[!is.na(ld$start), ]
    res <- merge(snps, ld, by = "SNP")
    res$chrname <- NULL
    return(res)

  } else {
    
    if(level > 0) {
      
      cover.all <- list(closestGenes(snps, genes, "cover", level -1))
      up.all    <- list(closestGenes(snps, genes, "up", level -1))
      down.all  <- list(closestGenes(snps, genes, "down", level -1))
      
    } else if(level == 0) {
      
      cover.all <- closestGenes(snps, genes, "cover")
      up.all    <- closestGenes(snps, genes, "up")
      down.all  <- closestGenes(snps, genes, "down")
      
      # keep covering or, if not exists, the closer of the up / down genes
      res <- cover.all
      res <- vectorElements(res)
      check.idx <- is.na(cover.all$start)
      replace.up <- check.idx & (abs(up.all$end - up.all$BP) <= abs(down.all$start - down.all$BP))
      # can contain NA which is invalid in a selection vector. When down is NA, set TRUE to select up.
      replace.up[is.na(replace.up)] <- is.na(down.all$start)[is.na(replace.up)]
      replace.down <- check.idx & !replace.up
      res[replace.up, ] <- up.all[replace.up, ]
      res[replace.down, ] <- down.all[replace.down, ]
      
      # remove columns with both genename and geneid NA (see closest genes function)
      res <- res[!(is.na(res$geneid) & is.na(res$genename)), ]
      return(res[, !(colnames(res) %in% c("chromo", "chrname.1", "chromo.1", "chrname"))])
      
    } else {
      
      # annotate next updown genes including overlapping genes
      level <- 0 
      cover.all <- list(closestGenes(snps, genes, "cover", level))
      up.all    <- list(closestGenes(snps, genes, "up", level))
      down.all  <- list(closestGenes(snps, genes, "down", level))
    
      # loop as long as overlapping genes exist in any direction. 
      while(TRUE) {
        level <- level +1

        cat(paste("Calculating the", level+1, "-th closest genes (search for overlap)...\n"))
        cover <- closestGenes(snps, genes, "cover", level)
        up    <- closestGenes(snps, genes, "up", level)
        down  <- closestGenes(snps, genes, "down", level)

        # check for overlap
        # set to NA where genes do no overlap with level 0
        up[na.set(up$end < up.all[[1]]$start, TRUE), ] <- NA
        down[na.set(down$start > down.all[[1]]$end, TRUE), ] <- NA

        if(all(is.na(down)) && all(is.na(up)) && all(is.na(down))) {
          break
        } else {
          cover.all <- c(cover.all, list(cover))
          up.all    <- c(up.all, list(up))
          down.all  <- c(down.all, list(down))
        }
      }
    }

    cat("Formatting results...\n")
    if(print.format) {
      res <- data.frame(
                snps, 
                "CoveringGenes" = prettyCluster(cover.all), 
                "UpstreamGenes" = prettyCluster(up.all), 
                "DownstreamGenes" = prettyCluster(down.all)
             )
      # remove chromosome mapping columns
      return(res[, colnames(res) != "chrname"])
    } else {
      # combine all snp - gene assignments in one long table
      res.df <- list2df(c(cover.all, up.all, down.all))
      # remove columns with both genename and geneid NA (see closest genes function)
      res.df <- res.df[!(is.na(res.df$geneid) & is.na(res.df$genename)), ]
      return(res.df[, !(colnames(res.df) %in% c("chromo", "chromo.1", "chrname.1", "chrname"))])
    }
	
  }

}

snp2gene.LD <- function(snps, biomart.config = biomartConfigs$hsapiens, use.buffer = FALSE, by.genename = TRUE, run.parallel = FALSE, ...) 
  snp2gene(snps = snps, biomart.config = biomart.config, use.buffer = use.buffer, by.genename = by.genename, run.parallel = run.parallel, by.ld = TRUE, ...)

snp2gene.prox <- function(snps, level = -1, biomart.config = biomartConfigs$hsapiens, use.buffer = FALSE, by.genename = TRUE, print.format = FALSE)
  snp2gene(snps, print.format = print.format , biomart.config = biomart.config, use.buffer = use.buffer, by.genename = by.genename, level = level)

