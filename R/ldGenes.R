ldGenes <- function(snps, genes, run.parallel = FALSE, ...) {

  if(!is.numeric(snps$BP) | !is.numeric(genes$start) | !is.numeric(genes$end))
    stop("Base positions have to be numeric")
  
  if(is.null(list(...)$map)) {
    map <- NULL
    ped <- NULL
    map.dat <- NULL
  } else {
    map <- list(...)$map
    ped <- list(...)$ped
    map.dat <- readMapfile(map)
  }
  
  # get positions of SNPs according to genotype source
  # snps$BP matches genes positions, snps$BP.mapped matches genotype positions
  snps.mapped <- getSnpsByRS(as.vector(snps$SNP), map = map.dat)
  snps <- vectorElements(merge(snps, snps.mapped, by = "SNP", suffixes = c("", ".mapped")))
  genes <- genes[order(genes$start, na.last = NA), ]
  
  if(run.parallel) {
	require(parallel)
    myapply <- mclapply
  } else {
    myapply <- lapply
  }
  
  # get genes for each SNP (all in a 1 MB window + nearest)
  snps.genes <- lapply(
    1:nrow(snps), 
    function(idx) {
      bp <- snps[idx, "BP"]
      genes.chr <- genes[genes$chrname == snps[idx, "CHR"], ]
      nearest.bp <- min(abs(genes.chr$start - bp), abs(genes.chr$end - bp))
      return(unique(rbind(
        # window genes
        genes.chr[genes.chr$end >= (bp -500000) & genes.chr$start <= (bp +500000), ], 
        # nearest gene
        genes.chr[abs(genes.chr$start - bp) == nearest.bp | abs(genes.chr$end - bp) == nearest.bp, ] 
      )))
     }
  )
  names(snps.genes) <- snps$SNP
  
  # get SNPs in genes
  # snps.genes.snps: snp.genes.snps$rs-id$genename contains rs-ids of all snps in that gene (max. 100)
  snps.genes.snps <- myapply(
              1:length(snps$SNP), 
              function(idx) {
                snp <- as.vector(snps$SNP)[idx] 
                shift <- as.vector(snps$BP - snps$BP.mapped)[idx]
                genes <- snps.genes[[snp]]
                
                # get snps within window (when map.dat is NULL, uses HapMap retrieval)
                ldsnps.all <- getSnpsByWin(
                    genes$chrname[1], 
                    min(genes$start - shift), 
                    max(genes$end - shift), 
                    map = map.dat
                )
                # return only snps within genes
                return(mapply(
                    function(id, start, end) {
                      ldsnps.gene <- ldsnps.all[as.vector(ldsnps.all$BP) >= start & as.vector(ldsnps.all$BP) <= end, ]
                      ldsnps.gene <- as.vector(ldsnps.gene[order(ldsnps.gene$BP), "SNP"])
                      return(pruneVec(ldsnps.gene, 100))
                    }, 
                    as.vector(genes$genename),  # dummy, only used for mapply to return a named list
                    genes$start - shift, 
                    genes$end - shift,
                    SIMPLIFY = FALSE
                ))
              }
          )
  names(snps.genes.snps) <- snps$SNP
  
  # now make a list with names = snps and their corresponding genotypes
  cat("Loading Genotypes...\n")
  snps.all <- unique(c(
    unlist(snps.genes.snps),    # snps from genes
    names(snps.genes.snps)      # always include the query SNP itself
  ))
  if(run.parallel) {
    snps.blocks <- makeBlocks(snps.all, detectCores())
    genos.genes <- unlist(myapply(snps.blocks, getGenotypes, remove.homozygous = TRUE, as.geno.objects = TRUE, map = map, ped = ped), recursive = FALSE)
  } else {
    genos.genes <- getGenotypes(snps = snps.all, remove.homozygous = TRUE, as.geno.objects = TRUE, map = map, ped = ped)
  }

  # calc LD
  res <- myapply(
            # for all query SNPs
            snps$SNP, 
            function(snp) {
              genes.snps <- snps.genes.snps[[snp]]
              ld <- sapply(
                  names(genes.snps), 
                  function(gene) {
                    cat(paste("Calculating LD:", snp, "<-> ", gene, "(", length(genes.snps[[gene]]), "snps )", "\n"))
                    genos.gene <- genos.genes[names(genos.genes) %in% genes.snps[[gene]]]
                    geno.snp <- genos.genes[[snp]]
                    if(is.null(geno.snp) | length(genos.gene) < 2)
                      return(c(
                        ld.max = 0,
                        ld.mean = 0,
                        ld.sdev = 0
                      ))
                    lds <- unlist(sapply(genos.gene, LD, geno.snp)["R^2",]) 
                    return(c(
                            ld.max = round(max(lds), digits = 4),
                            ld.mean = round(mean(lds), digits = 4),
                            ld.sdev = round(sd(lds), digits = 4)
                        ))
                  }
              )

              ldgenes.res <- cbind(snps.genes[[snp]], t(ld))
              ldgenes.res$SNP <- snp
#               ldgenes.res$summary <- paste(names(genes.snps),
#                                         "( max =", 
#                                         round(ld["ld.max", ], digits = 4),
#                                         ", mean =",  
#                                         round(ld["ld.mean", ], digits = 4), 
#                                         ", quart =",  
#                                         round(ld["ld.quart", ], digits = 4), 
#                                         ")"
#                                       )
              return(ldgenes.res)
            }
        )

  return(list2df(res))
  
}

