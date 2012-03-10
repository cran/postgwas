data.regionalplot.ld <-function(regions.df, ld.options, cores) {
  
    cat("Extracting SNPs for LD computation: ")
    # for each region, remove snps so that we have at most max.snps.per.window per region (evenly selected)
    regions.df <- regions.df[order(regions.df$BP.mapped), ]
    # a vector of rs ids
    ldsnps <- unique(unlist(tapply(
               regions.df$SNP, 
               factor(regions.df$packet.snp), 
               pruneVec, 
               ld.options$max.snps.per.window
              )))
    genos <- getGenotypes(
          ldsnps, 
          ld.options$pop.id, 
          ped = ld.options$ped, 
          map = ld.options$map, 
          remove.homozygous = TRUE
        )

    # precomute LD (gives a list of r-square matrices for each region)
    return(tapply(
      regions.df$SNP, 
      factor(regions.df$packet.snp), 
      function(region.snps) {
        # select genos that belong to snps from that region
        genos.region <- genos[names(genos) %in% unique(as.vector(region.snps))]
        if(length(genos.region) < 2) return(NULL)
        cat("Precalculating ", length(genos.region)^2, "LD pairs for current region\n")
        if(cores > 1) {
          rsq <- LDparallel(as.data.frame(genos.region), cores)
        } else {
          rsq <- LD(as.data.frame(genos.region))$"R^2"
        }
      }
    ))
  
}