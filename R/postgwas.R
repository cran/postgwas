postgwas <- function(
              gwas.resultfile, 
              suggestive.p = 1e-06, 
              genomewide.p = 5e-08, 
              biomart.config = biomartConfigs$hsapiens, 
			  GOpackagename = "org.Hs.eg.db", 
			  ped = NULL, 
              map = NULL
            ) {

  if(getRversion() >= '2.15.1') 
	  utils::globalVariables(c("postgwas.buffer.genes", "postgwas.buffer.snps", "postgwas.buffer.genes.regionalplot", "postgwas.buffer.snps.regionalplot", "postgwas.buffer.exons.regionalplot"))
  
  suppressWarnings(rm(postgwas.buffer.snps, postgwas.buffer.genes, postgwas.buffer.genes.regionalplot, postgwas.buffer.ld.regionalplot, pos = ".GlobalEnv"))
  suppressWarnings(rm(postgwas.buffer.snps, postgwas.buffer.genes, postgwas.buffer.genes.regionalplot, postgwas.buffer.ld.regionalplot, pos = "postgwasBuffer"))
  
  if(!is.null(GOpackagename)) {
    if(!all(c("GOSim", GOpackagename) %in% installed.packages()[, "Package"])) {
     cat(paste("Either package GOSim or", GOpackagename, "is not installed. Network will not be colorized by overrepresented GO terms.\n"))
      vertexcolor.GO.enrich <- NULL
    } else {
      vertexcolor.GO.enrich <- GOpackagename
    }
  }
  
  
  cat("\n\n****************** PREPARING SNPS ******************\n\n")
  
  snps <- readPvalFiles(gwas.resultfile, read.columns = c("SNP", "P"))
  if(any(duplicated(snps$SNP)))
    stop("SNPs are not unique in the GWAS result file (check for several models in the file, e.g. DOM/REC/GENO...).\n")
  
  # check for rs-IDs
  if(sum(grepl("rs\\d+", snps$SNP)) < nrow(snps) / 2)
    warning("Most SNPs are not named as dbSNP rs-ids. Make sure that the supplied IDs match the biomart config.\n")
  
  
  cat("\n\n****************** SNP TO GENE MAPPING ******************\n")
  
  cat("\nAnnotating genes (proximity)...\n")
  tryCatch(
    {
      s2g <- snp2gene.prox(
        snps[snps$P <= suggestive.p, ], 
        level = 0, 
        biomart.config = biomart.config, 
        use.buffer = TRUE
        );
      s2g <- s2g[, c("SNP", "CHR", "BP", "P", "geneid", "genename", "start", "end", "direction")]
    },
    error = function(e) cat(paste("snp2gene.prox: An error occured, proceeding with the next function.  ", e))
  )
  
  tryCatch(
    if(!is.null(ped) && !is.null(map)) {
      cat("\nAnnotating genes (LD)...\n")
      s2g.LD <- snp2gene.LD(
        snps[snps$P <= suggestive.p, ], 
        biomart.config = biomart.config, 
        use.buffer = TRUE, 
        ped = ped, 
        map = map
      )
      s2g.LD <- s2g.LD[s2g.LD$ld.max > 0.3, c("SNP", "CHR", "BP", "P", "geneid", "genename", "start", "end", "ld.max", "ld.mean", "ld.sdev")]
      s2g.LD$direction <- "LD"
      if(exists("s2g")) {
        s2g$ld.max <- NA
        s2g$ld.mean <- NA
        s2g$ld.sdev <- NA
        s2g <- rbind(s2g, s2g.LD)
      } else {
        s2g <- s2g.LD
      }
    },
    error = function(e) cat(paste("snp2gene.LD: An error occured, proceeding with the next function.  ", e))
  )
  
  if(exists("s2g")) 
    write.table(s2g, paste("postgwas_snp2gene_p", suggestive.p, ".csv", sep = ""), row.names = FALSE, sep = "\t")
  
  
  
  cat("\n\n****************** MANHATTAN PLOT ******************\n\n")
  
  tryCatch(
    manhattanplot(
      gwas.resultfile = gwas.resultfile, 
      highlight.logp = -log10(genomewide.p), 
      highlight.win = 25000, 
      biomart.config = biomart.config
    ), 
    error = function(e) cat(paste("manhattanplot: An error occured, proceeding with the next function.  ", e))
  )
  
  cat("\n\n****************** REGIONAL PLOTS ******************\n\n")
  
  tryCatch(
    {
      cat("Select representative SNPs (lowest P) within 50 kb windows for regionalplots...\n");
      snps <- removeNeighborSnps(snps[snps$P <= suggestive.p, c("SNP", "P")], maxdist = 200000, biomart.config = biomart.config, use.buffer = TRUE);
      if(nrow(snps) <= 0)
        stop("The dataset does not contain SNPs with suggestive significance.\n");
      
      if(!is.null(ped) && !is.null(map)) {
        ld.options <- list(
                        ped = ped,
                        map = map,
                        max.snps.per.window = ceiling(200 / log(nrow(snps) +1)), 
                        rsquare.min = 0.2, 
                        show.rsquare.text = FALSE
                      )
      } else  {
        ld.options <- NULL
      };
	  suppressWarnings(rm(postgwas.buffer.snps, pos = ".GlobalEnv"));
	  suppressWarnings(rm(postgwas.buffer.snps, pos = "postgwasBuffer"));
      cat(paste("Regionalplots of", nrow(snps), "loci...\n"));
      regionalplot(
        snps = snps, 
        gwas.resultfiles = gwas.resultfile, 
        ld.options = ld.options, 
        biomart.config = biomart.config, 
        use.buffer = TRUE
      )
    },
    error = function(e) cat(paste("regionalplot: An error occured, proceeding with the next function.  ", e))
  )

  
  if(exists("s2g")) {
  
    s2g <- s2g[, c("SNP", "genename", "P")]
    s2g <- s2g[order(s2g$P), ]
    s2g <- s2g[!duplicated(s2g$genename), ] # is sorted: keeps better p-values
    s2g <- vectorElements(s2g)
    
    if(nrow(s2g) <= 0)
      stop("No SNPs were mapped to genes.\n")
    
    cat("\n\n****************** NETWORK ANALYSIS ******************\n")
    cat("\nDownloading network data...\n")
    tryCatch(
      {
        network.dom <- getInteractions.domains(filter.ids = s2g$genename, max.occurence = NULL, toFile = "postgwas.interaction.download.domains");
        network.ppi <- getInteractions.ppi(filter.ids = s2g$genename, includeInteractors = FALSE, toFile = "postgwas.interaction.download.ppi");
        if(nrow(network.dom) > 0) {
          network.dom$label <- network.dom$domain.name
          network <- network.dom
        };
        if(nrow(network.ppi) > 0) {
          network.ppi$label <- "ppi";
          network <- network.ppi
        };
        if(nrow(network.dom) > 0 && nrow(network.ppi) > 0) {
          network <- rbind(network.ppi, network.dom[, c("genename.x", "genename.y", "label")]);
        };
        net <- gwas2network(
          gwas.mapped.genes = s2g, 
          network = network, 
          prune = "gwasonly", 
          vertexcolor.GO.enrich = vertexcolor.GO.enrich, 
          biomart.config = biomart.config, 
          use.buffer = TRUE
        )
      },
      error = function(e) cat(paste("gwas2network: An error occured.  ", e))
      )
    
  } 

  cat("\n\n****************** THANKS FOR USING POSTGWAS ******************\n")
  cat("Your output files have been produced.\n")
  cat("Try also out the individual postgwas functions - \n")
  cat("they offer additional features and customization possibilities.\n\n")
  
}
