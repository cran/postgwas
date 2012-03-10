# g is a graph
# universe.idmap all genes defined by geneid and genename column
# geneid.col names one of these columns that are used in the graph as names (and match the appropriate universe.idmap column)
# value: list, graph = the colorized graph, legend.v = a data frame of vertices for the legend
gwas2network.vertexcolor.enrich <- function(
  g, 
  universe.idmap, 
  geneid.col, 
  custom.sel = F,
  toFile = NULL
  ) {
  
  # set base color (has to be set because recolorization depends on the previous color)
  V(g)$color <- "black"
  
  if(geneid.col == "genename") {
    # translate vertex names to entrez IDs
    genes.gwas <- merge(data.frame(genename = V(g)$name), universe.idmap, all.x = TRUE, sort = TRUE)
    # make names unique, keep NAs
    genes.gwas <- genes.gwas[!duplicated(genes.gwas$genename) | is.na(genes.gwas$genename), ]
  } else {
    genes.gwas <- data.frame(geneid = V(g)$name, P = V(g)$P)
  }
  
  signif.terms <- NULL
  legend.v <- NULL
  
  tryCatch(
    signif.terms <- gwasGOenrich(
      gwas.mapped.genes.entrez = genes.gwas, 
      gene.universe.entrez = universe.idmap, 
      # do not use kol.smirn                      p.default = 0.05, 
      alpha = NULL, 
      toFile = toFile
      ), 
    error = function(e) {
      cat(paste("Could not set vertexcolor by GO enrichment:  ", e))
    }
    )
  
  if(is.null(signif.terms) || nrow(signif.terms) <= 0) {
    
    cat("No significantly enriched GO-terms found.\n")
    
  } else {
    
    # list: names = color, values = index in signif.terms
    clr.idx <- list()
    
    if(custom.sel) {
      cat("The most significant GO terms:\n")
      nterms <- if(nrow(signif.terms) > 20) 20 else nrow(signif.terms)
      cat(paste(
        1:nterms, 
        signif.terms[1:nterms, 1], 
        signif(as.numeric(signif.terms[1:nterms, 2])), 
        signif.terms[1:nterms, 3], 
        sep = "\t", collapse = "\n"
        ))
      while({
        cat("\nChoose terms to colorize? [Y/N]\n"); 
        ans <- scan(n = 1, what = character(), quiet = TRUE); 
        length(ans) > 0 && ans %in% c("Y", "y")
      }) {
        cat("\nSpecify a color [e.g. red, blue, darkgreen, lightsalmon3, ...]\n")
        clr <- scan(n = 1, what = character(), quiet = TRUE)
        if(length(clr) < 1 || !(clr %in% colors())) {
          cat("\nInvalid color!\n")
          next
        }
        cat(paste("\nSpecify a term number [1-", nterms, "]\n", sep = ""))
        idx <- scan(n = 1, what = integer(), quiet = TRUE)
        if(length(idx) < 1 || idx > nterms || idx < 1) {
          cat("\nInvalid index!\n")
          next
        }
        clr.idx[[clr]] <- idx
      }
      
    } else {
      
      # colorize the most significant
      clr.idx[["red"]] <- 1
      if(nrow(signif.terms) > 1)
        clr.idx[["blue"]] <- 2
      if(nrow(signif.terms) > 2)
        clr.idx[["darkgreen"]] <- 3
      
    }
    
    if(length(clr.idx) > 0) {
      
      # names : color, value: genenames (vertices) to colorize
      clr.vnames <- lapply(
        clr.idx, 
        function(idx) {
          entrez <- strsplit(signif.terms[idx, 5], split = ",", fixed = TRUE)[[1]]
          if(geneid.col == "genename") {
            # map back to names
            entrez2name <- merge(data.frame(geneid = entrez), universe.idmap, all.x = TRUE, sort = TRUE)
            return(entrez2name[!duplicated(entrez2name$genename) | is.na(entrez2name$genename), "genename"])
          } else {
            return(entrez)
          }
        }
        )
      
      # define legend vertices
      legend.v <- data.frame(
        name = c(paste(signif.terms[unlist(clr.idx), 1], signif.terms[unlist(clr.idx), 2], signif.terms[unlist(clr.idx), 3], sep = "\n"), "multi-annotation"), 
        label = c(paste(signif.terms[unlist(clr.idx), 1], signif.terms[unlist(clr.idx), 2], signif.terms[unlist(clr.idx), 3], sep = "\n"), "multi-annotation"),
        weight = mean(V(g)$weight), 
        P = NA,
        SNP = NA, 
        degree = 0,
        marked = FALSE, 
        label.color = "black",
        color = c(names(clr.idx), "saddlebrown")
        )
      
      for(clr in names(clr.vnames)) {
        V(g)$color[V(g)$color != "black" & V(g)$name %in% clr.vnames[[clr]]] <- "saddlebrown"
        V(g)$color[V(g)$color == "black" & V(g)$name %in% clr.vnames[[clr]]] <- clr
      }
      
    } # end valid index exist
    
  } # end signif terms found
  
  return(list(graph = g, legend.v = legend.v))
}



# g is a graph
# geneid.col names one of these columns that are used in the graph as names
# value: new colorized graph
gwas2network.vertexcolor.regex <- function(
  g, 
  regex, 
  geneid.col, 
  bm.config = NULL, 
  ds = bm.init.genes(bm.config)
  ) {
  
  if(!all(names(regex) %in% colors()))
    stop("Vertex color regex list contains invalid color names. See colors() for all valid color names.\n")
  
  cat("Retrieving gene ontology annotation for vertex color...\n")
  go <- bm.genes.goterms(
    config = bm.config, 
    ds = ds, 
    filter.name = if(geneid.col == "genename") bm.config$gene$filter$name else bm.config$gene$filter$id, 
    filter.val = V(g)$name
    )
  
  # names : color, value: a list of genenames (vertices) to colorize
  clr.vnames <- lapply(
    regex, 
    function(expr) unique(go[grep(expr, as.vector(go$"goterm.name"), ignore.case = T), geneid.col])
    )
  
  V(g)$color <- "black"
  
  for(clr in names(clr.vnames)) {
    V(g)$color[V(g)$color != "black" & V(g)$name %in% clr.vnames[[clr]]] <- "saddlebrown"
    V(g)$color[V(g)$color == "black" & V(g)$name %in% clr.vnames[[clr]]] <- clr
  }
  
  return(g)
}


# returns a data frame 'go_id', 'P', 'Term', 'genecount', 'genefrac', 'genes'
# where the genes column contains all entrez genes for that term, separated by ','
gwasGOenrich <- function(
  gwas.mapped.genes.entrez, 
  gene.universe.entrez = bm.genes(use.buffer = TRUE), 
  p.default = NULL, # when NULL use Fisher test, otherwise Kolmog. SMirn. as described in GOSim
  alpha = 0.05, # when NULL, return all terms with p < 0.01, otherwise all exceeding bonferroni correction
  max.term.abundance = 0.2,
  method = "elim", # topgo method, classic for all terms
  toFile = nextFilename("gwasGOenrich", "csv")
  ) {
  
  if(!"GOSim" %in% names(installed.packages()[, 'Package'])) {
    stop("This functionality requires the 'GOSim' package to be installed.\n")
  } else {
    suppressPackageStartupMessages(require(GOSim, quietly = TRUE))
  }
  
  if(!is.null(p.default)) {
    # all = TRUE ensures that universe also contains query genes
    gene.universe.entrez <- merge(gene.universe.entrez, gwas.mapped.genes.entrez, all = T)[, c("geneid", "P")]
    gene.universe.entrez[is.na(gene.universe.entrez$P), "P"] <- p.default
    gene.universe.entrez <- na.omit(gene.universe.entrez)
    univ <- gene.universe.entrez$P
    names(univ) <- gene.universe.entrez$geneid
    
    gwas <- gwas.mapped.genes.entrez$P
    names(gwas) <- gwas.mapped.genes.entrez$geneid
    gwas <- gwas[!is.na(names(gwas))]
  } else {
    gwas <- as.vector(as.character(na.omit(gwas.mapped.genes.entrez$geneid)))
    univ <- c(gwas, as.vector(as.character(na.omit(gene.universe.entrez$geneid))))
  }
  
  if(length(gwas) < 1)
    return(NULL)
  
  capture.output(enrich.res <- GOenrichment(gwas, univ, method = method))
  
  if(length(enrich.res$p.values) < 1) 
    return(NULL)
  
  if(!is.null(alpha)) {
    num.tests <- length(unique(enrich.res$GOTerms$go_id))
    signif.terms <- enrich.res$p.values[as.vector(as.numeric(enrich.res$p.values)) < alpha/num.tests]
  } else {
    signif.terms <- enrich.res$p.values
  }
  signif.terms <- data.frame(go_id = names(signif.terms), P = signif.terms)
  signif.terms <- merge(signif.terms, enrich.res$GOTerms)[, 1:3]
  
  signif.terms <- as.data.frame(t(apply(
    signif.terms, 
    1, 
    function(row) 
      return(c(
        row, 
        genecount = length(enrich.res$genes[[row[1]]]),
        genes = paste(enrich.res$genes[[row[1]]], collapse = ",")
        ))
    )), stringsAsFactors = FALSE)
  
  signif.terms <- signif.terms[order(as.vector(as.numeric(signif.terms$P))), ]
  signif.terms$genefrac <- format.pval(as.numeric(signif.terms$genecount) / length(univ), digits = 1)
  signif.terms <- signif.terms[as.vector(as.numeric(signif.terms$genefrac)) < max.term.abundance, ]
  signif.terms <- signif.terms[as.vector(as.numeric(signif.terms$genecount)) > 1, ]
  
  if(!is.null(toFile) && is.character(toFile))
    write.table(signif.terms, toFile[1], row.names = FALSE, sep = "\t")
  
  return(signif.terms)
}
