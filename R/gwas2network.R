gwas2network <- function(
                  gwas.mapped.genes, 
                  network, 
                  prune = "shared", 
                  remove.superhubs = nrow(network) > 100,
                  p.default = 0.06, 
                  min.transparency = 0.25, 
                  max.transparency = 0.75, 
                  custom.layout = FALSE,
                  vertexcolor.GO.enrich = "org.Hs.eg.db",
                  vertexcolor.GO.regex = NULL,
                  max.communities = 500, 
                  edge.weight.fun = function(p1, p2, degree1, degree2, weight) {
                    return(((-log10(p1) +1) * (-log10(p2) +1))^2 * weight)
                  },
                  edge.weight.collapse.fun = function(weights, labels) {
                    return(mean(weights))
                  },
                  png.maxedges = 2500, 
                  file.verbosity = 2, 
                  biomart.config = biomartConfigs$hsapiens, 
                  use.buffer = FALSE, 
                  cores = 1
                ) {
  
  if(is.null(network) || nrow(network) == 0)
    stop("Parameter 'network' does not contain data.\n")
  
  if(ncol(network) < 2)
    stop("Parameter 'network' has to contain at least two columns.\n")
  
  if(any(c("net.weight.temp", "gwas2network.genename.x", "gwas2network.genename.y") %in% colnames(network))) 
    stop("Parameter 'network' may not contain the columns 'net.weight.temp', 'gwas2network.genename.x'', 'gwas2network.genename.y'.\n")
  
  # make sure weight column exists and does not contain NA (set to 1 instead)
  if(is.null(network$weight))
    network$weight <- 1
  network$weight[is.na(network$weight)] <- 1
  network <- vectorElements(network)
  
  if(!is.function(edge.weight.fun) || !is.function(edge.weight.collapse.fun))
    stop("Arguments 'edge.weight.fun' and 'edge.weight.collapse.fun' have to be functions.\n")
  
  if(!all(c("p1", "p2", "degree1", "degree2", "weight") %in% formalArgs(edge.weight.fun)))
    stop("Argument 'edge.weight.fun' has an incomplete argument list, must contain 'p1', 'p2', 'degree1', 'degree2', 'weight'.\n")
  
  if(!all(c("weights", "labels") %in% formalArgs(edge.weight.collapse.fun)))
    stop("Argument 'edge.weight.collapse.fun' has an incomplete argument list, must contain 'weights', 'labels'.\n")
  
  if(!is.null(cores) && cores > 1) {
    if(!'parallel' %in% installed.packages()[, 'Package'])
      stop("Function parameter cores > 1 requires the package 'parallel' to work correctly.\n")
    suppressPackageStartupMessages(stopifnot(require(parallel, quietly = TRUE)))
  }
  
  if(is.null(gwas.mapped.genes) || !is.data.frame(gwas.mapped.genes) || nrow(gwas.mapped.genes) == 0) {
    stop("Argument 'gwas.mapped.genes' has to be a nonempty data frame.\n")
  } else {
    gwas.mapped.genes <- vectorElements(gwas.mapped.genes)
    if(any(c("net.weight.temp", "gwas2network.genename.x", "gwas2network.genename.y") %in% colnames(gwas.mapped.genes))) 
      stop("Argument 'gwas.mapped.genes' may not contain the columns 'net.weight.temp', 'gwas2network.genename.x'', 'gwas2network.genename.y'.\n")
    if(is.null(gwas.mapped.genes$P)) 
      stop("Argument 'gwas.mapped.genes' has to contain a column 'P'.\n")
    if(!any(c("geneid", "genename") %in% colnames(gwas.mapped.genes))) 
      stop("Argument 'gwas.mapped.genes' has to contain either a column 'geneid' or 'genename'.\n")
    if(all(c("geneid", "genename") %in% colnames(gwas.mapped.genes))) {
      # determine the kind of identifier used in the network data
      if(sum(as.vector(na.omit(gwas.mapped.genes[, "geneid"])) %in% c(network[, 1], network[, 2]))
          > 
         sum(as.vector(na.omit(gwas.mapped.genes[, "genename"])) %in% c(network[, 1], network[, 2]))
         )
        geneid.col <- "geneid"
      else
        geneid.col <- "genename"
    } else {
      geneid.col <- if("geneid" %in% colnames(gwas.mapped.genes)) "geneid" else "genename"
    }
    cat(paste(
      "Using", 
      geneid.col, 
      "as primary identifier.", 
      round(sum(unique(gwas.mapped.genes[, geneid.col]) %in% c(network[, 1], network[, 2])) / length(unique(gwas.mapped.genes[, geneid.col])) * 100, 2),
      "% of query genes are contained in the network before processing.\n"
    ))
  }
  
  if(!is.null(vertexcolor.GO.regex)) {
    if(!is.list(vertexcolor.GO.regex))
      stop("Argument 'vertexcolor.GO.regex' is not a list.\n")
    if(!all(names(vertexcolor.GO.regex) %in% colors()))
      stop("Argument 'vertexcolor.GO.regex' contains invalid color names. See colors() for a list of valid color names.\n")
  }

  if(is.null(vertexcolor.GO.regex) && !is.null(vertexcolor.GO.enrich)) {
    if(!is.vector(vertexcolor.GO.enrich) || length(vertexcolor.GO.enrich) > 1 || typeof(vertexcolor.GO.enrich) != "character")
      stop("Argument 'vertexcolor.GO.enrich' has to be a character string defining an Annotation package name (see help pages).\n")
    if(!("GOSim" %in% installed.packages()[, "Package"]))
      stop("Vertex colorization by GO enrichment requires the GOSim package to be installed.\n")
    if(biomart.config$gene$attr$id != "entrezgene") 
      warning("Vertex colorization by significance normally requires entrez gene ids (specify properly in the biomart configuration)\n")
    require(GOSim)
    # load proper GO organism data
    success <- do.call(library, list(package = vertexcolor.GO.enrich, logical.return = TRUE))
    if(!success)
      stop("The GeneOnotology annotation package supplied in argument 'vertexcolor.GO.enrich' is not installed or cannot be loaded.\n")
    # reinit GOSim when intitalized with a different GO annotation package
    org.name.new <- get(grep("ORGANISM", ls(paste("package", vertexcolor.GO.enrich, sep = ":")), value = TRUE)[1])
    org.dat.new <- get(grep("GO", ls(paste("package", vertexcolor.GO.enrich, sep = ":")), value = TRUE)[1])
    if(sub("human", "Homo sapiens", GOSimEnv$organism) == org.name.new) {
      # everything is already preloaded
    } else {
      setEvidenceLevel(organism = org.name.new, gomap = org.dat.new)
      calcICs(path.package("GOSim"))
      setOntology("BP", DIR = path.package("GOSim"))
    }
  }
  
  if(!("igraph" %in% installed.packages()[, "Package"])) {
	if(interactive() && readline("Package 'igraph' is not installed. Try to install? [Y/N] ") %in% c("Y", "y")) {
		install.packages("igraph", repos = "http://cran.us.r-project.org")
	} else {
		stop("Package igraph is required for gwas2network but not installed.\n")		
	}
  } else {
	# there is a legacy package 'igraph0' that masks igraph and leads to errors
	# ensure that our igraph is first on search path 
	# and remember its current position for restore at function exit
    igraph.pos.orig <- which(search() == "package:igraph")
	suppressWarnings(detach("package:igraph", force = TRUE))
	suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = 2, quietly = TRUE)))
  }

  
  # we need an idmap for two purposes:
  # when ids are used in the graph, we want to display names in the plot
  # we need a universe of all genes for the organism, which is represented by the (complete) idmap
  cat("Retrieving all gene IDs and names from biomart for the organism (serves as gene universe for enrichment and ID to name mapping) ...\n")
  delayedAssign("ds", bm.init.genes(biomart.config))
  idmap <- bm.genes(
             config = biomart.config, 
             ds = ds, 
             use.buffer = use.buffer
           )[, 1:2]
  idmap.nona <- na.omit(idmap)
  
  
  g.dat <- data.gwas2network(
              gwas.mapped.genes = gwas.mapped.genes, 
              network = network, 
              geneid.col = geneid.col,
              prune = prune, 
              remove.superhubs = remove.superhubs,
              p.default = p.default, 
              edge.weight.fun = edge.weight.fun,
              edge.weight.collapse.fun = edge.weight.collapse.fun,
              toFile = file.verbosity > 1
            )
  e <- g.dat$e
  v <- g.dat$v

  g <- graph.data.frame(e, directed = FALSE, v)
  

  ######### vertex color and graph plot #########
  
  # GO enrichment
  if(!is.null(vertexcolor.GO.enrich) || !is.null(vertexcolor.GO.regex)) {
    cat("Setting vertex color...\n")
    if(!is.null(vertexcolor.GO.regex)) {
      g <- gwas2network.vertexcolor.regex(
        g = g, 
        regex = vertexcolor.GO.regex, 
        geneid.col = geneid.col,
        bm.config = biomart.config,
        ds = ds
        )
      legend.v <- NULL
    } else {
      clr.res <- gwas2network.vertexcolor.enrich(
        g = g, 
        universe.idmap = idmap, 
        geneid.col = geneid.col, 
        custom.sel = FALSE,
        toFile = if(file.verbosity > 2) "gwasGOenrich.csv" else NULL
        )
      g <- clr.res$graph
      legend.v <- clr.res$legend.v
    }
  } else {
    # no GO colors
    V(g)$color <- "black"
    legend.v <- NULL
  }
  
  
  ######### set transparency #########

  vw.max <- max(V(g)$weight)
  vw.min <- min(V(g)$weight)
  ew.max <- max(E(g)$weight)
  ew.min <- min(E(g)$weight)

  if(vw.min < 0)
    stop("Vertex weights < 0 found. Aborting.\n")

  # determine alpha (normalize weight to [0,1] and then expand to [min.transparence,max.transparency])
  if(vw.max == vw.min) {
    # avoid division by NULL
    V(g)$alpha <- max.transparency
  } else {
    V(g)$alpha <- ((V(g)$weight - vw.min) / (vw.max - vw.min)) * (max.transparency - min.transparency) + min.transparency
  }
  
  if(ew.max == ew.min) {
    E(g)$alpha <- max.transparency
  } else {
    E(g)$alpha <- ((E(g)$weight - ew.min) / (ew.max - ew.min)) * (max.transparency - min.transparency) + min.transparency
  }
  
  V(g)$label.color <- rgb(red = 0, green = 0, blue = 0, alpha = V(g)$alpha)
  E(g)$color       <- rgb(red = 0, green = 0, blue = 0, alpha = E(g)$alpha^1.7)
  # the vertex transparency has to be adjusted for each graph seperately, because the base vertex color can change after enrichment


  ######### define scale vertices #########
  # scale veritces have no valid gene names and may only be added after enrichment and so
  scale.vertices <- list(nv = 0)
  if(vw.max != vw.min) {
    max.scale <- max(2, ceiling(vw.max))
    alpha <- ((2:max.scale - min(vw.min, 2)) / (vw.max - min(vw.min, 2))) * (max.transparency - min.transparency) + min.transparency
    alpha[alpha > 1] <- 1
    alpha[alpha < 0] <- 0
    color <- rgb(0, 0, 0, alpha = alpha)
    scale.vertices <- list(
      nv = length(2:max.scale),
      name = format.pval(10^(-(2:max.scale))), 
      label = format.pval(10^(-(2:max.scale))),
      weight = 2:max.scale, 
      P = 10^(-(2:max.scale)), 
      SNP = format.pval(10^(-(2:max.scale))),
      marked = FALSE,
      label.color = color,
      color = color
    )
  }
  

  cat("Plotting the entire graph...\n")  
  g.whole <- g
  v.color          <- t(col2rgb(V(g.whole)$color)) / 255
  V(g.whole)$color <- rgb(red = v.color[, "red"], green = v.color[, "green"], blue = v.color[, "blue"], alpha = V(g.whole)$alpha)
  # translate labels to names when IDs were used
  if(geneid.col == "geneid")
    g.whole <- graph.id2name(g.whole, idmap.nona)
  g.whole <- addLegendVertices(g.whole, legend.v)
  # add scale vertices
  scale.vertices$graph <- g.whole
  g.whole <- do.call(add.vertices, scale.vertices)
  if(file.verbosity > 0) {
    gwas2network.plot(g.whole, file = "gwas2network.pdf", title = "gwas2network plot", v.weight.max = vw.max)
    if(ecount(g.whole) <= png.maxedges)
      gwas2network.plot(g.whole, device = png, file = "gwas2network.png", res = 200, v.weight.max = vw.max)
  }
  
  
  g.comps <- decompose.graph(g, mode = "strong", min.vertices = 2)  
  if(file.verbosity > 0) {
    cat("Plotting connected components / communities ...\n")
    g.comps <- g.comps[order(sapply(g.comps, vcount), decreasing = TRUE)]
  }
  
  for(g.comp in g.comps) {
    
    if(max.communities > 1) {
      
      cat("Detecting communities of current component...\n")
      
      membership <- spinglass.community(
                      graph = g.comp, 
                      weights=E(g.comp)$weight, 
                      spins = max.communities,
					  implementation = ifelse(ew.min < 0, "neg", "orig")
                    )$membership
      
      if(!is.null(cores) && cores > 1) {
        myapply <- mclapply
      } else {
        myapply <- lapply
      }
      
      # prepare modules (subgs is a ready-to plot list of graph objects)
      subgs <- myapply(
        names(table(factor(membership))), 
        function(mod.id, mc.cores) { # mc.cores isa dummy argument to catch the argument when normal apply is used
          subg <- induced.subgraph(g.comp, V(g.comp)[membership == mod.id])
          # spinglass can produce empty modules
          if(vcount(subg) < 1)
            next
          # community specific GO enrichment
          if(!is.null(vertexcolor.GO.enrich)) {
            V(subg)$color <- "black"
            cat(paste("Setting vertex color of community", as.numeric(mod.id) +1, "...\n"))
            clr.res <- gwas2network.vertexcolor.enrich(
              g = subg, 
              universe.idmap = idmap, 
              geneid.col = geneid.col, 
              custom.sel = FALSE,
              toFile = if(file.verbosity > 2) nextFilename("gwasGOenrichCommunities", "csv") else NULL
            )
            subg <- clr.res$graph
            legend.v <- clr.res$legend.v
            # set transparency 
            v.color       <- t(col2rgb(V(subg)$color)) / 255
            V(subg)$color <- rgb(red = v.color[, "red"], green = v.color[, "green"], blue = v.color[, "blue"], alpha = V(subg)$alpha)
          }
          
          # translate labels to names when IDs were used
          if(geneid.col == "geneid")
            subg <- graph.id2name(subg, idmap.nona)
          subg <- addLegendVertices(subg, legend.v)
          # add scale vertices
          scale.vertices$graph <- subg
          subg <- do.call(add.vertices, scale.vertices)
          return(subg)
        },
        mc.cores =cores
      )
      
      for(subg in subgs) {
        if(file.verbosity > 2) {
          # dump subgraph data
          v.dump <- merge(v, data.frame(name = V(subg)$name, color = V(subg)$color))
          write.table(v.dump, nextFilename("gwas2networkGraphCommunityVertices", "csv"), row.names = FALSE, sep = "\t")
          e.dump <- e[e$gwas2network.genename.x %in% as.vector(get.edgelist(subg)) & e$gwas2network.genename.y %in% as.vector(get.edgelist(subg)), ]
          write.table(e.dump, nextFilename("gwas2networkGraphCommunityEdges", "csv"), row.names = FALSE, sep = "\t")
        }
        
        # plot modules
        if(custom.layout && interactive() && ecount(subg) > 0) {
          cat("Drag vertices to the custom layout. When done, press Enter (before closing the window).\n")
          subg.tkp <- subg
          V(subg.tkp)$frame.color <- V(subg.tkp)$color
          V(subg.tkp)$label.color <- V(subg.tkp)$color
          V(subg.tkp)$color <- "white"
          tkp.id <- suppressWarnings(tkplot(subg.tkp))
          scan(n = 1, quiet = T)
          layout <- tkplot.getcoords(tkp.id)
          tkplot.close(tkp.id)
          layout[, 2] <- -layout[, 2]
          layout <- layout / 2.5
        } else {
          layout <- layout.fruchterman.reingold(subg, area = vcount(subg)^3)
        }
        if(file.verbosity > 0) {
          gwas2network.plot(subg, layout = layout, title = "gwas2network plot (communities)", file = nextFilename("gwas2networkCommunities", "pdf"), v.weight.max = vw.max)
          if(ecount(subg) <= png.maxedges)
            gwas2network.plot(subg, layout = layout, device = png, file = nextFilename("gwas2networkCommunities", "png"), res = 200, v.weight.max = vw.max)
        }
      }
    } else {
      # translate labels to names when IDs were used
      if(geneid.col == "geneid")
        g.comp <- graph.id2name(g.comp, idmap.nona)
      v.color         <- t(col2rgb(V(g.comp)$color)) / 255
      V(g.comp)$color <- rgb(red = v.color[, "red"], green = v.color[, "green"], blue = v.color[, "blue"], alpha = V(g.comp)$alpha)
      g.comp <- addLegendVertices(g.comp, legend.v)
      # add scale vertices
      scale.vertices$graph <- g.comp
      g.comp <- do.call(add.vertices, scale.vertices)
      if(file.verbosity > 0) {
        gwas2network.plot(g = g.comp, file = nextFilename("gwas2networkConnectedComponents", "pdf"), title = "gwas2network plot (connected components)", v.weight.max = vw.max)
        if(ecount(g.comp) <= png.maxedges) 
          gwas2network.plot(g = g.comp, file = nextFilename("gwas2networkConnectedComponents", "png"), device = png, res = 200, v.weight.max = vw.max)
      }
    }

  }

  
  ######### dump graph data, return #########
  
  v.dump <- merge(v, data.frame(name = V(g)$name, color = V(g)$color, size = v$weight * 2))
  if(file.verbosity > 1) {
    write.table(v.dump, "gwas2networkGraphVertices.csv", row.names = FALSE, sep = "\t")
    write.table(e, "gwas2networkGraphEdges.csv", row.names = FALSE, sep = "\t")
  }
  
  # restore search path
  suppressWarnings(detach("package:igraph", force = TRUE))
  suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = igraph.pos.orig, quietly = TRUE)))
  
  return(g.whole)
  
}




addLegendVertices <- function(
                       graph, 
                       v.df, 
                       weight = mean(V(graph)$weight)
                     ) {
  if(is.null(v.df))
    return(graph)
  else
    return(add.vertices(
      graph, 
      nv = nrow(v.df),
      name = as.vector(v.df$name), 
      label = as.vector(v.df$label),
      weight = weight, 
      P = 10^(-weight), 
      SNP = v.df$label,
      marked = FALSE,
      label.color = as.vector(v.df$color),
      color = rgb(red = 0, green = 0, blue = 0, alpha = 0)
    ))
}

graph.id2name <- function(graph, idmap.nona) {
    V(graph)$label <- sapply(
                       V(graph)$label, 
                       function(x) {
                         mapped <- idmap.nona[idmap.nona$geneid == x, "genename"]
                         if(length(mapped) == 0 || mapped == "") x else mapped[1]
                       }
                     )
  return(graph)
}
