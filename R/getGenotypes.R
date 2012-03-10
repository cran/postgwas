getGenotypes <- function(
                  snps, 
                  pop.id = 2, 
                  ped = NULL, 
                  map = NULL, 
#                  format.tab = NA, 
                  remove.homozygous = FALSE, 
                  as.geno.objects = TRUE) 
                {

    ###################################
    ######### PARAMETER CHECK #########
    ###################################

    if(missing(snps))
      stop("Parameter 'snps' has to be specified")

    snps <- as.vector(snps)

    if(is.null(snps) | length(snps) < 1) {
      warning("No SNPs specified to retrieve genotypes for")
      return(NULL)
    }
    
    if( !is.null(ped) & !is.null(map) ) {
      use.pedmap <- TRUE
    } else {
      if( !is.null(pop.id) & is.numeric(pop.id) ) {
        use.pedmap <- FALSE
      } else {
        stop("Specify either a numeric HapMap population ID or location of ped/map files")
      }
    }
    
    gts <- NULL

    ###################################
    ####### PED / MAP RETRIEVAL #######
    ###################################
    if( use.pedmap ) 
      {
        cat("Reading genotype files ...\n")
        mapinfo <- readMapfile(map)
        snp.colmask <- mapinfo$SNP %in% snps

        # read ped line by line, extract SNPs of interest
        pedfile <- file(ped, open = "r")
        ped.line <- readLines(pedfile, n = 1)

        while( length(ped.line) != 0 ) {
#           if(is.na(format.tab)) format.tab <- isPedTabFormatted(ped.line, nrow(mapinfo))
#           if(format.tab) {
#             ped.line.vec <- unlist(strsplit(ped.line, split = "\t", fixed = TRUE))
#             # remove leading non-genotype columns
#             gts.line <- ped.line.vec[(length(ped.line.vec) - length(snp.colmask) +1):length(ped.line.vec)]      
#             # selected target snps
#             gts.line.selected <- gts.line[snp.colmask]
#           } else {
            ped.line.vec <- unlist(strsplit(ped.line, split = "\\s"))
            # remove leading non-genotype columns
            gts.line <- ped.line.vec[(length(ped.line.vec) - 2*length(snp.colmask) +1):length(ped.line.vec)]
            # select target snps and combine the two alleles for each SNP
            gts.line.selected <- sapply(which(snp.colmask) *2, function(idx) paste(gts.line[idx -1], gts.line[idx]))
#          }
          # format for 'genotype' class (recode numeric nucleotides to character)
          gts.line.selected <- sub(" ", "", gts.line.selected, fixed = TRUE)
          gts.line.selected <- gsub("1", "A", gts.line.selected, fixed = TRUE)
          gts.line.selected <- gsub("2", "C", gts.line.selected, fixed = TRUE)
          gts.line.selected <- gsub("3", "G", gts.line.selected, fixed = TRUE)
          gts.line.selected <- gsub("4", "T", gts.line.selected, fixed = TRUE)
          # add line to whole gts matrix and read next line from file
          gts <- rbind(gts, gts.line.selected)
          ped.line <- readLines(pedfile, n = 1)
        }
        close(pedfile)
        colnames(gts) <- mapinfo$SNP[mapinfo$SNP %in% snps]
    } 
    else {
    ###################################
    ####### HAPMART RETRIEVAL #########
    ################################### 
      # hapmap queries do not work with r package 'biomart' because the virtualSchemaName is not properly set (always 'default' instead of e.g. 'rel27_NCBI_Build36')
      # construct the RESTful XML query address for martservice
      # number of SNPs might be large, divide into 250-element blocks
      snps.blocks <- makeBlocks(snps, block.size = 250)
      for( snps.block in snps.blocks ) {
        cat("Retrieving genotypes from HapMart...\n")
        query.url <- url(paste("http://hapmap.ncbi.nlm.nih.gov/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20%20virtualSchemaName%20=%20%22rel27_NCBI_Build36%22%20formatter%20=%20%22TSV%22%20header%20=%20%220%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%3CDataset%20name%20=%20%22hm27_variation%22%20interface%20=%20%22default%22%20%3E%3CFilter%20name%20=%20%22marker_name%22%20value%20=%20%22", paste(snps.block, collapse = ","), "%22/%3E%3CFilter%20name%20=%20%22pop_id%22%20value%20=%20%22", pop.id, "%22/%3E%22/%3E%3CAttribute%20name%20=%20%22marker1%22%20/%3E%3CAttribute%20name%20=%20%22genotype%22%20/%3E%3C/Dataset%3E%3C/Query%3E", sep = ""))
        gts <- c(gts, readLines(query.url))
        close(query.url)
      }
      # extract rs numbers and genotypes as strings
      gts <- strsplit(gts, split = "\t", fixed = TRUE) 
      gts <- unlist(gts)
      gts.rsnumbers <- gts[seq(1, length(gts), by = 2)]
      gts <- gts[seq(2, length(gts), by = 2)]
      # convert genotype strings to vectors (containing single genotypes per individ)
      gts <- sapply(gts, function(x) { strsplit(x, split = " ", fixed = TRUE) })
      # make genotype vectors uniform length (appending NAs by selecting over-length subset)
      gts.maxcount <- max( sapply(gts, length) )
      gts <- sapply(array(gts), function(x) { subset(x, rep(TRUE, gts.maxcount)) })
      dimnames(gts) <- list(NULL, gts.rsnumbers)
    }

    # replace everything except nucleotides with NA
    gts[grep("[^ACGT]", gts)] <- NA

    if( remove.homozygous ) {
      # tabulate all gts gives number of homo- and heterozygous - so we just extract those that have more than one 'table' entry
      hom.het.hom.count <- apply(gts, 2, table)
      # apply returns matrix if all are same zygosity state (same number of 'table' entries), convert to list in this case
      if( class(hom.het.hom.count) == "matrix" ) 
        hom.het.hom.count <- as.data.frame(hom.het.hom.count)
      gts.zygosity <- lapply(hom.het.hom.count, length)
      gts <- gts[, unlist(gts.zygosity) != 1]
    }

  if( as.geno.objects ) {
    if(is.null(gts))
      return(NULL)
    if( remove.homozygous ) 
      if( class(hom.het.hom.count) != "matrix" )
        gts <- as.matrix(gts)
    # remove all SNPs that have only NA genotypes (converting to class genotype fails in this case)
    gts.cols.nacount <- apply(is.na(gts),2,sum)
    gts <- gts[, !(gts.cols.nacount == nrow(gts))]
    cat("Creating genotype objects...\n")
    geno <- lapply(as.data.frame(gts), function(x) { genotype(x, sep = "") })
    if(length(geno) == 1) {
      # genotype conversion by data frame does not retain names for a single element
      names(geno) <- mapinfo$SNP[mapinfo$SNP %in% snps]
    }
    return(geno)
  } else {
    return(gts)
  }
}

#isPedTabFormatted <- function(ped.line, snpcount) {
#  # guess "tab" or "space"
#  if(length(unlist(strsplit(ped.line, split = "\t", fixed = TRUE))) >= snpcount) return(TRUE)
#  else if(length(unlist(strsplit(ped.line, split = "\\s"))) >= 2 * snpcount) return(FALSE)
#  else stop("Pedfile format unknown - too few columns detected")
#}
