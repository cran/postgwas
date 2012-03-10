biomartConfigs <- list(
hsapiens = list(
  snp = list(
    mart = "snp", 
    dataset = "hsapiens_snp", 
    filter = list(
      id = "snp_filter"), 
    attr = list(
      id = "refsnp_id", 
      chr = "chr_name", 
      bp = "chrom_start")
  ), 
  gene = list(
    mart = "ensembl",
    dataset = "hsapiens_gene_ensembl", 
    filter = list(
      id = "entrezgene", 
      name = "hgnc_symbol",
      chr = "chromosome_name", 
      startbp = "start", 
      endbp = "end"), 
    attr = list(
      id = "entrezgene", 
      name = "hgnc_symbol", 
      chr = "chromosome_name", 
      startbp = "start_position", 
      endbp = "end_position", 
      strand = "strand",
      goterms = "name_1006",
      domains = "superfamily")
  ),
  exon = list(
    filter = list(
      geneid = "entrezgene", 
      chr = "chromosome_name"), 
    attr = list(
      chr = "chromosome_name",
      startbp = "exon_chrom_start", 
      endbp = "exon_chrom_end",
      gene.startbp = "start_position",
      gene.endbp = "end_position")
  )
), 
scerevisiae = list(
  snp = list(
    mart = "snp", 
    dataset = "scerevisiae_snp", 
    filter = list(
      id = "snp_filter"), 
    attr = list(
      id = "refsnp_id", 
      chr = "chr_name", 
      bp = "chrom_start")
  ), 
  gene = list(
    mart = "ensembl",
    dataset = "scerevisiae_gene_ensembl", 
    filter = list(
      id = "entrezgene", 
      name = "wikigene_name",
      chr = "chromosome_name", 
      startbp = "start", 
      endbp = "end"), 
    attr = list(
      id = "entrezgene", 
      name = "wikigene_name", 
      chr = "chromosome_name", 
      startbp = "start_position", 
      endbp = "end_position", 
      strand = "strand",
      goterms = "name_1006",
      domains = "superfamily")
  ),
  exon = list(
    filter = list(
      geneid = "entrezgene",
      chr = "chromosome_name"), 
    attr = list(
      chr = "chromosome_name", 
      startbp = "exon_chrom_start", 
      endbp = "exon_chrom_end",
      gene.startbp = "start_position",
      gene.endbp = "end_position")
  )
)
)

biomartConfigs$oanatinus <- biomartConfigs$scerevisiae
biomartConfigs$tguttata <- biomartConfigs$scerevisiae
biomartConfigs$fcatus <- biomartConfigs$scerevisiae
biomartConfigs$rnorvegicus <- biomartConfigs$scerevisiae
biomartConfigs$ggallus <- biomartConfigs$scerevisiae
biomartConfigs$ecaballus <- biomartConfigs$scerevisiae
biomartConfigs$pabelii <- biomartConfigs$scerevisiae
biomartConfigs$drerio <- biomartConfigs$scerevisiae
biomartConfigs$tnigroviridis <- biomartConfigs$scerevisiae
biomartConfigs$mdomestica <- biomartConfigs$scerevisiae
biomartConfigs$dmelanogaster <- biomartConfigs$scerevisiae
biomartConfigs$ptroglodytes <- biomartConfigs$scerevisiae
biomartConfigs$sscrofa <- biomartConfigs$scerevisiae
biomartConfigs$mmusculus <- biomartConfigs$scerevisiae
biomartConfigs$btaurus <- biomartConfigs$scerevisiae
biomartConfigs$cfamiliaris <- biomartConfigs$scerevisiae

biomartConfigs$oanatinus$snp$dataset <- "oanatinus_snp"
biomartConfigs$tguttata$snp$dataset <- "tguttata_snp"
biomartConfigs$fcatus$snp$dataset <- "fcatus_snp"
biomartConfigs$rnorvegicus$snp$dataset <- "rnorvegicus_snp"
biomartConfigs$ggallus$snp$dataset <- "ggallus_snp"
biomartConfigs$ecaballus$snp$dataset <- "ecaballus_snp"
biomartConfigs$pabelii$snp$dataset <- "pabelii_snp"
biomartConfigs$drerio$snp$dataset <- "drerio_snp"
biomartConfigs$tnigroviridis$snp$dataset <- "tnigroviridis_snp"
biomartConfigs$mdomestica$snp$dataset <- "mdomestica_snp"
biomartConfigs$dmelanogaster$snp$dataset <- "dmelanogaster_snp"
biomartConfigs$ptroglodytes$snp$dataset <- "ptroglodytes_snp"
biomartConfigs$sscrofa$snp$dataset <- "sscrofa_snp"
biomartConfigs$mmusculus$snp$dataset <- "mmusculus_snp"
biomartConfigs$btaurus$snp$dataset <- "btaurus_snp"
biomartConfigs$cfamiliaris$snp$dataset <- "cfamiliaris_snp"

biomartConfigs$oanatinus$gene$dataset <- "oanatinus_gene_ensembl"
biomartConfigs$tguttata$gene$dataset <- "tguttata_gene_ensembl"
biomartConfigs$hsapiens$gene$dataset <- "hsapiens_gene_ensembl"
biomartConfigs$fcatus$gene$dataset <- "fcatus_gene_ensembl"
biomartConfigs$rnorvegicus$gene$dataset <- "rnorvegicus_gene_ensembl"
biomartConfigs$ggallus$gene$dataset <- "ggallus_gene_ensembl"
biomartConfigs$ecaballus$gene$dataset <- "ecaballus_gene_ensembl"
biomartConfigs$pabelii$gene$dataset <- "pabelii_gene_ensembl"
biomartConfigs$drerio$gene$dataset <- "drerio_gene_ensembl"
biomartConfigs$tnigroviridis$gene$dataset <- "tnigroviridis_gene_ensembl"
biomartConfigs$mdomestica$gene$dataset <- "mdomestica_gene_ensembl"
biomartConfigs$dmelanogaster$gene$dataset <- "dmelanogaster_gene_ensembl"
biomartConfigs$ptroglodytes$gene$dataset <- "ptroglodytes_gene_ensembl"
biomartConfigs$sscrofa$gene$dataset <- "sscrofa_gene_ensembl"
biomartConfigs$mmusculus$gene$dataset <- "mmusculus_gene_ensembl"
biomartConfigs$btaurus$gene$dataset <- "btaurus_gene_ensembl"
biomartConfigs$cfamiliaris$gene$dataset <- "cfamiliaris_gene_ensembl"

biomartConfigs$fcatus$gene$attr$name <- "external_gene_id"
biomartConfigs$rnorvegicus$gene$name <- "rgd_symbol"
biomartConfigs$mmusculus$gene$name <- "mgi_symbol"



prepareBuffer <- function() {
  if(!("postgwasBuffer" %in% search())) {
    attach(what = NULL, name = "postgwasBuffer")
  } else {
    if(environmentIsLocked(pos.to.env(which(search() == "postgwasBuffer")[1])))
      stop("Cannot access the environment for buffer data, 'postgwasBuffer' (is locked). Try setting use.buffer = FALSE or open a new R session.\n")
  }
  if(!exists("postgwas.buffer.genes", where = "postgwasBuffer"))
	  assign("postgwas.buffer.genes", NULL, pos = "postgwasBuffer")
  if(!exists("postgwas.buffer.snps", where = "postgwasBuffer"))
	  assign("postgwas.buffer.snps", NULL, pos = "postgwasBuffer")
  if(!exists("postgwas.buffer.exons.regionalplot", where = "postgwasBuffer"))
	  assign("postgwas.buffer.exons.regionalplot", NULL, pos = "postgwasBuffer")
  if(!exists("postgwas.buffer.genes.regionalplot", where = "postgwasBuffer"))
	  assign("postgwas.buffer.genes.regionalplot", NULL, pos = "postgwasBuffer")
}

# clearPostgwasBuffer <- function() {
#   cat("Clearing buffer variables if existing...\n")
#   while("postgwasBuffer" %in% search()) detach("postgwasBuffer")
#   while(length(getAnywhere("postgwas.buffer.snps")$where) > 0) 
#     rm("postgwas.buffer.snps", pos = getAnywhere("postgwas.buffer.snps")$where)
#   while(length(getAnywhere("postgwas.buffer.genes")$where) > 0)
#     rm("postgwas.buffer.genes", pos = getAnywhere("postgwas.buffer.gene")$where)
#   while(length(getAnywhere("postgwas.buffer.exons.regionalplot")$where) > 0)
#     rm("postgwas.buffer.exons.regionalplot", pos = getAnywhere("postgwas.buffer.exons.regionalplot")$where)
#   while(length(getAnywhere("postgwas.buffer.genes.regionalplot")$where) > 0)
#     rm("postgwas.buffer.genes.regionalplot", pos = getAnywhere("postgwas.buffer.genes.regionalplot")$where)
#   while(length(getAnywhere("postgwas.buffer.ld.regionalplot")$where) > 0) 
#     rm("postgwas.buffer.ld.regionalplot", pos = getAnywhere("postgwas.buffer.ld.regionalplot")$where)
# }

# deny retrieval and usage of exons when one of the exon config entries is NA
exons.avail <- function(biomart.config = biomartConfigs$hsapiens, use.buffer = FALSE) {
  # ok when valid buffer exists
  if(use.buffer && exists("postgwas.buffer.exons.regionalplot") && !is.null(postgwas.buffer.exons.regionalplot)) return(TRUE)
  # deny when attributes are not set
  if(is.null(biomart.config$exon) || any(is.na(unlist(biomart.config$exon)))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# args: a predefined species as character string, e.g. "hsapiens", or a list with all required biomart field names. Required fields are listed by calling names(biomartConfigs$hsapiens)
# value: a biomaRt dataset object
bm.init.snps <- function(species = NULL) {

  if(is.null(species))
    stop(paste("The following predefined configs are available (use as quoted string in the function call):", paste(names(biomartConfigs), collapse = ", ")))

  if(is.list(species)) { # user-defined biomart attribute list
    config <- species
  } else if(is.null(biomartConfigs[[species]])) {
      stop(paste("Biomart config not available for selected species. Available predefined configs:", paste(names(biomartConfigs), collapse = ", ")))
  } else {
    config <- biomartConfigs[[species]]
  }

  cat("Initializing biomart connection (SNPs)...\n")
  ds <- useMart(biomart = config$snp$mart, dataset = config$snp$dataset)

  # check attribs exist
  if(sum(listAttributes(ds)[, 1] %in% unlist(config$snp$attr)) < length(config$snp$attr)) 
    stop(paste("One of the biomart attributes [", paste(config$snp$attr, collapse = "; "), "] in dataset [", config$snp$dataset, "] does not exist"))
  # filters
  if(sum(listFilters(ds)[, 1] %in% unlist(config$snp$filter)) < length(config$snp$filter)) 
    stop(paste("One of the biomart filters [", paste(config$snp$filter, collapse = "; "), "] in dataset [", config$snp$dataset, "] does not exist"))

  return(ds)
}

# args: a predefined species as character string, e.g. "hsapiens", or a list with all required biomart field names. Required fields are listed by calling names(biomartConfigs$hsapiens)
# value: a biomaRt dataset object
bm.init.genes <- function(species = NULL) {

  if(is.null(species))
    stop(paste("The following predefined configs are available (use as quoted string in the function call):", paste(names(biomartConfigs), collapse = ", ")))

  if(is.list(species)) { # user-defined biomart attribute list
    config <- species
  } else if(is.null(biomartConfigs[[species]])) {
      stop(paste("Biomart config not available for selected species. Available predefined configs:", paste(names(biomartConfigs), collapse = ", ")))
  } else {
    config <- biomartConfigs[[species]]
  }

  cat("Initializing biomart connection (genes)...\n")
  ds <- useMart(biomart = config$gene$mart, dataset = config$gene$dataset)

  # check attribs exist
  if(sum(listAttributes(ds)[, 1] %in% unlist(config$gene$attr)) < length(config$gene$attr)) 
    stop(paste("One of the biomart attributes", paste(config$gene$attr, collapse = "; "), "in dataset", config$gene$dataset, "does not exist"))
  # exon attribs
  if(exons.avail(config, FALSE) && sum(listAttributes(ds)[, 1] %in% unlist(config$exon$attr)) < length(config$exon$attr)) 
    stop(paste("One of the biomart attributes", paste(config$exon$attr, collapse = "; "), "in dataset", config$gene$dataset, "does not exist"))
  # gene filters
  if(sum(listFilters(ds)[, 1] %in% unlist(config$gene$filter)) < length(config$gene$filter)) 
    stop(paste("One of the biomart filters", paste(config$gene$filter, collapse = "; "), "in dataset", config$gene$dataset, "does not exist"))

  return(ds)
}



# returns a data frame with the columns 
# "refsnp_id", "chrname", "chrom_start" 
# in the listed order. Columns are cast to vector. 
bm.snps <- function(
              config = biomartConfigs$hsapiens, 
              ds = bm.init.snps(config), 
              filter.val, 
              use.buffer = FALSE,
              remove.dupl = TRUE) {
  attrnames <- c(config$snp$attr$id, config$snp$attr$chr, config$snp$attr$bp)
  if(use.buffer && exists("postgwas.buffer.snps") && !is.null(postgwas.buffer.snps)) {
    cat("Using buffer data (SNPs)\n")
    if(!all(attrnames == colnames(postgwas.buffer.snps)))
      stop("Error in snp buffer data - clear buffer with rm(postgwas.buffer.snps, inherits = TRUE) or set use.buffer = FALSE")
    else
      snps <- postgwas.buffer.snps[postgwas.buffer.snps[, config$snp$attr$id] %in% filter.val, ]
  } else {
    snps <- getBM(
              attributes = attrnames, 
              filters = config$snp$filter$id, 
              values = filter.val, 
              mart = ds
            )
    if(use.buffer) {
      prepareBuffer()
      assign("postgwas.buffer.snps", snps, pos = "postgwasBuffer")
    }
  }
  colnames(snps) <- c("refsnp_id", "chrname", "chrom_start")
  
  # rare cases of multiple positions (chromosomes) per SNP in biomart can be removed (e.g. MHC assemblies)
  # keep the smallest chromsome name by character length
  # duplicate removes elements with larger index -> sort by character length
  if(remove.dupl) 
    snps <- snps[!(snps$refsnp_id %in% snps$refsnp_id[duplicated(snps$refsnp_id)]), ]

  snps$chrname <- as.vector(snps$chrname)
  return(snps)
}


# returns a data frame with the columns 
# "geneid", genename", "chrname", "start", "end"
# in the listed order.
# We can have strange chromosomes for genes like LRG genes, but this should not matter
bm.genes <- function(
                        config = biomartConfigs$hsapiens, 
                        ds = bm.init.genes(config), 
                        use.buffer = FALSE) {
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$chr, config$gene$attr$startbp, config$gene$attr$endbp)
  if(use.buffer && exists("postgwas.buffer.genes") && !is.null(postgwas.buffer.genes)) {
    cat("Using buffer data (genes)\n")
    genes <- postgwas.buffer.genes
    if(!all(attrnames %in% colnames(genes)))
      stop("Error in genes buffer data - clear buffer with rm(postgwas.buffer.genes, inherits = TRUE) or set use.buffer = FALSE")
  } else {
    genes <- getBM(attributes = attrnames, mart = ds)
    if(use.buffer) {
      prepareBuffer()
      assign("postgwas.buffer.genes", genes, pos = "postgwasBuffer")
    }
  }
  names(genes) <- c("geneid", "genename", "chrname", "start", "end")
  genes$geneid <- as.vector(genes$geneid)
  genes$genename <- as.vector(genes$genename)
  genes$chrname <- as.vector(genes$chrname)
  return(genes)
}


# returns a data frame with the columns 
# id, name, chr, startbp, endbp, strand 
# in the listed order.
bm.genes.regionalplot <- function(
                            config = biomartConfigs$hsapiens, 
                            ds = bm.init.genes(config), 
                            chr,      # position can be vectors for multiple regions
                            bp.start, 
                            bp.end, 
                            use.buffer = FALSE) {
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$chr, config$gene$attr$startbp, config$gene$attr$endbp, config$gene$attr$strand)
  if(use.buffer && exists("postgwas.buffer.genes.regionalplot") && !is.null(postgwas.buffer.genes.regionalplot)) {
    cat("Using buffer data (genes)\n")
    genes <- postgwas.buffer.genes.regionalplot
    if(!all(attrnames %in% colnames(genes)))
      stop("Error in genes buffer data - clear buffer with rm(postgwas.buffer.genes.regionalplot, inherits = TRUE) or set use.buffer = FALSE")
  } else {
    genes <- unique(list2df(lapply(
      1:length(chr), 
      function(idx) 
        getBM( 
          attributes = attrnames, 
          filters = c(config$gene$filter$chr, config$gene$filter$startbp, config$gene$filter$endbp),
          values = list(chr[idx], bp.start[idx], bp.end[idx]),
          mart = ds 
        )
    )))
    if(use.buffer) {
      prepareBuffer()
      assign("postgwas.buffer.genes.regionalplot", genes, pos = "postgwasBuffer")
    }
  }
  names(genes) <- c("id", "name", "chromo", "start", "end", "strand")
  genes$id <- as.vector(genes$id)
  genes$name <- as.vector(genes$name)
  genes$mean <- genes$start + abs((genes$start - genes$end)/2)
  return(genes)
}


# using the buffer, this function returns all exons in the buffer!
# filtering for chromosome is mandatory because highly variable regions like MHC can have multiple chromos for the same position
bm.exons <- function(
                   config = biomartConfigs$hsapiens, 
                   ds = bm.init.genes(config), 
                   filter.val = list(geneid = c(""), chromo = c("")), 
                   use.buffer = FALSE) {
  names(filter.val) <- c(config$exon$filter$geneid, config$exon$filter$chr)
  attrnames <- c(config$exon$attr$chr, config$exon$attr$startbp, config$exon$attr$endbp, config$exon$attr$gene.startbp, config$exon$attr$gene.endbp)
  if(use.buffer && exists("postgwas.buffer.exons.regionalplot") && !is.null(postgwas.buffer.exons.regionalplot)) {
    cat("Using buffer data (exons)\n")
    if(!all(attrnames == colnames(postgwas.buffer.exons.regionalplot))) {
      stop("Error in exon buffer data - clear buffer with rm(postgwas.buffer.exons.regionalplot, inherits = TRUE) or set use.buffer = FALSE")
    } else {
      exons <- postgwas.buffer.exons.regionalplot
    }
  } else {
    exons <- getBM(
      attributes = attrnames, 
      filters = c(config$exon$filter$geneid, config$exon$filter$chr), 
      values = filter.val, 
      mart = ds
    )
    if(use.buffer) {
      prepareBuffer()
      assign("postgwas.buffer.exons.regionalplot", exons, pos = "postgwasBuffer")
    }
  }
  names(exons) <- c("chromo", "start", "end", "genestart", "geneend")
  return(exons)
}

bm.genes.goterms <- function(
                       config = biomartConfigs$hsapiens, 
                       ds = bm.init.genes(config),
                       filter.name = biomartConfigs$hsapiens$gene$filter$name, 
                       filter.val
                     ) {
  
  if(is.null(config$gene$attr$goterms)) 
    stop("Biomart configuration is incomplete - the 'goterms' element has to be specified for the current operation\n")
  
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$goterms)
  goterms <- getBM( 
               attributes = attrnames, 
               filters = filter.name,
               values = filter.val,
               mart = ds 
             )
  colnames(goterms) <- c("geneid", "genename", "goterm.name")
  return(goterms)
}

bm.genes.domains <- function(
                       config = biomartConfigs$hsapiens, 
                       ds = bm.init.genes(config),
                       filter.name = biomartConfigs$hsapiens$gene$filter$name, 
                       filter.val
                     ) {
  
  if(is.null(config$gene$attr$domains)) 
    stop("Biomart configuration is incomplete - the 'domains' element has to be specified for the current operation\n")
  
  attrnames <- c(config$gene$attr$id, config$gene$attr$name, config$gene$attr$domains)
  domains <- getBM( 
               attributes = attrnames, 
               filters = filter.name,
               values = filter.val,
               mart = ds 
             )
  colnames(domains) <- c("geneid", "genename", "domain.name")
  domains <- domains[domains$domain.name != "", ]
  return(domains)
}


bm.id2name <- function(
                filter.ids = NULL,
                config = biomartConfigs$hsapiens, 
                ds = bm.init.genes(config), 
                use.buffer = FALSE
              ) {

  if(is.null(filter.ids))
    stop("Argument filter.ids has to be specified.\n")
  
  if(use.buffer && exists("postgwas.buffer.genes") && !is.null(postgwas.buffer.genes)) {
    if(!all(c(config$gene$attr$id, config$gene$attr$name) %in% colnames(postgwas.buffer.genes)))
      stop("Error in genes buffer data - clear buffer with rm(postgwas.buffer.genes, inherits = TRUE) or set use.buffer = FALSE")
    cat("Using buffer data (genes)\n")
  } else {
    use.buffer <- FALSE
  }
  
  if(is.data.frame(filter.ids) & config$gene$attr$id %in% colnames(filter.ids)) {
    if(use.buffer) {
      map <- merge(filter.ids, postgwas.buffer.genes, all.x = TRUE)
      map <- map[, c(colnames(filter.ids), config$gene$attr$name)]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$id,
        values = filter.ids[, config$gene$attr$id],
        mart = ds 
      )
      map <- merge(filter.ids, map, all.x = TRUE)
    }
  } else {
    if(use.buffer) {
      map <- postgwas.buffer.genes[
               postgwas.buffer.genes[, config$gene$attr$id] %in% as.vector(filter.ids), 
               c(config$gene$attr$id, config$gene$attr$name)
             ]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$id,
        values = filter.ids,
        mart = ds 
      )
    }
  }
  # make names unique, keep NAs
  map <- map[!duplicated(map[, config$gene$attr$name]) | is.na(map[, config$gene$attr$name]), ]
  return(map)

}


bm.name2id <- function(
                filter.names = NULL,
                config = biomartConfigs$hsapiens, 
                ds = bm.init.genes(config), 
                use.buffer = FALSE
              ) {
  if(is.null(filter.names))
    stop("Argument filter.names has to be specified.\n")
  
  if(use.buffer && exists("postgwas.buffer.genes") && !is.null(postgwas.buffer.genes)) {
    if(!all(c(config$gene$attr$id, config$gene$attr$name) %in% colnames(postgwas.buffer.genes)))
      stop("Error in genes buffer data - clear buffer with rm(postgwas.buffer.genes, inherits = TRUE) or set use.buffer = FALSE")
    cat("Using buffer data (genes)\n")
  } else {
    use.buffer <- FALSE
  }
  
  if(is.data.frame(filter.names) & config$gene$attr$name %in% colnames(filter.names)) {
    if(use.buffer) {
      map <- merge(filter.names, postgwas.buffer.genes, all.x = TRUE)
      map <- map[, c(colnames(filter.names), config$gene$attr$id)]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$name,
        values = filter.names[, config$gene$attr$name],
        mart = ds 
        )
      map <- merge(filter.names, map, all.x = TRUE)
    }
  } else {
    if(use.buffer) {
      map <- postgwas.buffer.genes[
        postgwas.buffer.genes[, config$gene$attr$name] %in% as.vector(filter.names), 
        c(config$gene$attr$id, config$gene$attr$name)
        ]
    } else {
      map <- getBM(
        attributes = c(config$gene$attr$id, config$gene$attr$name), 
        filters = config$gene$filter$name,
        values = filter.names,
        mart = ds 
      )
    }
  }
  # make IDs unique, keep NAs
  map <- map[!duplicated(map[, config$gene$attr$id]) | is.na(map[, config$gene$attr$id]), ]
  return(map)
}


bm.remapSnps <- function(
                  file, 
                  config = biomartConfigs$hsapiens, 
                  ds = bm.init.snps(config), 
                  use.buffer = FALSE, 
                  toFile = paste(file, ".remapped", sep = "")
                ) {
  
  f <- na.omit(read.table(file, header=TRUE))
  if(!("SNP" %in% colnames(f)))
    stop("Source file does not contain a column 'SNP' or is not readable by the read.table function.\n")
  
  cat("Retrieving SNP positions from biomart, this can take a while.\n")
  snps.bm <- bm.snps(
    filter.val = f$SNP, 
    config = config, 
    use.buffer = use.buffer
  )
  
  colnames(snps.bm)[colnames(snps.bm) %in% c("chrname", "chrom_start")] <- c("CHR", "BP")
  f.new <- merge(f, snps.bm, by.x = "SNP", by.y = "refsnp_id", suffixes = c(".original", ""))
  
  if(!is.null(toFile) && is.character(toFile) && length(toFile) == 1)
    write.table(f.new, toFile, row.names = F, col.names = T)
  
}
