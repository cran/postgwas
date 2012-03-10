manhattanplot <- function(
                    gwas.resultfile, 
                    highlight.logp = c(6, 7.3), 
                    highlight.win = 50000, 
                    highlight.color = NULL, 
                    highlight.text = "genes",
                    highlight.cex = 1,
                    max.y = NULL,
                    reduce.dataset = TRUE, 
                    plot.title = gwas.resultfile, 
                    biomart.config = biomartConfigs$hsapiens, 
                    use.buffer = FALSE, 
                    toFile = nextFilename("manhattanplot", "pdf")
                  ) {

  ######### prepare data #########

  if(is.null(gwas.resultfile) || !is.character(gwas.resultfile))
    stop("Argument gwas.resultfile has to be set to a path/filename (character string).\n")
  
  gwasdat <- readPvalFiles(gwas.resultfile)
  
  if( !("SNP" %in% names(gwasdat) && "CHR" %in% names(gwasdat) && "BP" %in% names(gwasdat) && "P" %in% names(gwasdat) ) )
    stop("Make sure your data frame contains columns CHR, BP, and P")

  if(is.null(highlight.logp) || length(highlight.logp) < 1)
    stop("Argument highlight.logp has to be set.\n")
  
  if(is.null(highlight.color)) {
    highlight.color <- rev(rainbow(length(highlight.logp)))
  } else if(!all(highlight.color %in% colors())) {
    stop("Highlighting colors not defined or unknown (see colors())")
  }

  # preprocess data
  gwasdat <- vectorElements(gwasdat[, c("SNP", "CHR", "BP", "P")])
  gwasdat$P  <- as.numeric(gwasdat$P)
  gwasdat$BP <- as.numeric(gwasdat$BP)
  gwasdat <- na.omit(gwasdat)
  gwasdat$logp <- -log10(gwasdat$P)
  gwasdat <- gwasdat[gwasdat$P > 0 & gwasdat$P < 1, ]
  
  if(reduce.dataset) {
    gwasdat <- gwasdat[gwasdat$logp > runif(nrow(gwasdat))^2 | gwasdat$logp < 1e-04, ]
    gwasdat <- gwasdat[gwasdat$logp > (runif(nrow(gwasdat))/1.1)^2 | gwasdat$logp < 1e-04, ]
  }
  
  # remove entries that exceed max.y
  if(!is.null(max.y) && is.numeric(max.y) && length(max.y) <= 1)
    gwasdat <- gwasdat[gwasdat$logp < max.y, ]
  
  # translate base position to coords (-> gwasdat$pos), making chromosmes sequential
  
  # use mapping to generate numeric chromosomes (order and number does not matter)
  chr.map  <- data.frame(CHR = sort(unique(gwasdat$CHR)), CHR.mapped = 1:length(unique(gwasdat$CHR)))
  gwasdat  <- merge(gwasdat, chr.map)
  chr.max.base <- tapply(gwasdat$BP, factor(gwasdat$CHR.mapped), max)  # for each chromosome the bp length (vector name = chrnames)
  chr.shift <- sapply(
                  names(chr.max.base), 
                  function(chr) sum(as.double(chr.max.base[as.numeric(names(chr.max.base)) <= as.numeric(chr)])) - chr.max.base[chr],
                  USE.NAMES = FALSE
                )
  for(i in names(chr.shift))
    gwasdat[gwasdat$CHR.mapped == i, "pos"] <- gwasdat[gwasdat$CHR.mapped == i, "BP"] + chr.shift[as.character(i)]

  # reorganize (and recycle) parameters when multiple p thresholds are given
  param <- cbind(highlight.logp, highlight.win, highlight.color, highlight.cex)
  if(nrow(param) > 1) {
    # highlight tresholds have to be sorted
    param <- param[order(highlight.logp), ]
    highlight.logp <- param[, "highlight.logp"]
    highlight.win <- param[, "highlight.win"]
    highlight.color <- param[, "highlight.color"]
    highlight.cex <- param[, "highlight.cex"]
  }

  
  ######### set highlighted SNP color and shape #########
  
  cat("Identifying highlighted regions...\n")
  clr <- data.frame(CHR = sort(unique(as.vector(gwasdat$CHR))))
  clr[seq(1, nrow(clr), 2), "color"] <- "grey10"
  clr[seq(2, nrow(clr), 2), "color"] <- "grey50"
  gwasdat <- merge(gwasdat, clr)
  
  gwasdat$x <- gwasdat$pos / max(gwasdat$pos)
  gwasdat$y <- gwasdat$logp / max(gwasdat$logp)
  
  # highl only contains the PEAK data point of highlighted regions
  highl <- mapply(
    function(logp, win, col, cex) {
      highl    <- gwasdat[gwasdat$logp > as.numeric(logp), ]
      highl    <- removeNeighborSnps(highl, as.numeric(win))
      highl$CHR.mapped <- as.numeric(as.vector(highl$CHR.mapped))
      highl$BP <- as.numeric(as.vector(highl$BP))
      highl$P <- as.numeric(as.vector(highl$P))
      highl$color <- rep(col, nrow(highl))
      highl$cex <- rep(cex, nrow(highl))
      return(highl)
    }, 
    highlight.logp, 
    highlight.win, 
    highlight.color,
    highlight.cex, 
    SIMPLIFY = FALSE
  )
  # for multiple p threshs: make intervals (currently highl only selected for upper bound, i.e. < p)
  if(length(highlight.logp) > 1)
    for(i in 1:(length(highlight.logp)-1))
      highl[[i]] <- highl[[i]][highl[[i]]$logp <= highlight.logp[i+1], ]

  # shape of all ordinaty data points: circle
  gwasdat$shape <- 1
  
  # all highlighted data points: triangle
  # also set color and cex of highlighted data points
  mapply(
    function(highl, win, color, cex) {
      if(nrow(highl) > 0) {
        for(idx in 1:(nrow(highl))) {
          gwasdat.highl.idx <- which(
                                 as.numeric(gwasdat$CHR.mapped) == as.numeric(highl[idx, "CHR.mapped"]) & 
                                 gwasdat$BP >= as.numeric(highl[idx, "BP"]) - win & 
                                 gwasdat$BP <= as.numeric(highl[idx, "BP"]) + win
                               )
          gwasdat[gwasdat.highl.idx, "color"] <<- color
          gwasdat[gwasdat.highl.idx, "shape"] <<- 2
          gwasdat[gwasdat.highl.idx, "cex"] <<- cex
        }
      }
    },
    highl, 
    as.numeric(highlight.win), 
    highlight.color,
    highlight.cex
  )
  
  
  ######### plot points and genes #########
  if(!is.null(toFile) && is.character(toFile) && length(toFile) == 1)
    pdf(toFile, width = 10, height = 5)
  
    
  # use npc: everything scaled to papersize percent
  cat("Scatterplot...\n")
  vp.plot <- viewport(width = 0.8, height = 0.8, name = "plot")
  pushViewport(vp.plot)
  
  gwasdat.plain <- gwasdat[gwasdat$shape == 1, ]
  grid.points(
    x = unit(gwasdat.plain$x, "npc"), 
    y = unit(gwasdat.plain$y, "npc"),
    pch = gwasdat.plain$shape,
    size = unit(0.0001, "npc"),
    gp = gpar(col = gwasdat.plain$color)
  )

  gwasdat.highl <- gwasdat[gwasdat$shape != 1, ]
  if(nrow(gwasdat.highl) > 0) {
    
    grid.points(
      x = unit(gwasdat.highl$x, "npc"), 
      y = unit(gwasdat.highl$y, "npc"),
      pch = 25,
      size = unit(0.005 * as.numeric(gwasdat.highl$cex)^1.7, "npc"),
      gp = gpar(col = gwasdat.highl$color, fill = gwasdat.highl$color)
    )
    
    # text (gene) annot
    for(highl.df in highl) {
      if(!is.null(highlight.text) && nrow(highl.df) > 0) {
        
        if(highlight.text == "SNP") {
          grid.text(
              label = highl.df$SNP, 
              x = unit(highl.df$x, "npc"), 
              y = unit(highl.df$y, "npc") + unit(0.013 * as.numeric(highl.df$cex), "npc"),
              just = "bottom",
              gp = gpar(col = highl.df$color, fill = highl.df$color, cex = 0.7 * as.numeric(highl.df$cex))
          )
        } else {
          
          cat("Geneplot...\n")
          colnames(highl.df)[colnames(highl.df) == "CHR.mapped"] <- "chrmp"
          highl.df <- snp2gene.prox(snps = highl.df, by.genename = TRUE, level = 1, biomart.config = biomart.config, use.buffer = use.buffer)
          
          if(nrow(highl.df) > 0) {
            
            highl.df.left <- highl.df[highl.df$direction == "down", ]
            if(nrow(highl.df.left) > 0)
              grid.text(
                  label = highl.df.left$genename, 
                  x = unit(highl.df.left$x, "npc") - unit(0.005 * as.numeric(highl.df.left$cex)^2, "npc"), 
                  y = unit(highl.df.left$y, "npc"),
                  just = "right",
                  gp = gpar(col = highl.df.left$color, fill = highl.df.left$color, cex = 0.7 * as.numeric(highl.df.left$cex))
              )
            
            highl.df.right <- highl.df[highl.df$direction == "up", ]
            if(nrow(highl.df.right) > 0)
              grid.text(
                  label = highl.df.right$genename, 
                  x = unit(highl.df.right$x, "npc") + unit(0.005 * as.numeric(highl.df.right$cex)^2, "npc"), 
                  y = unit(highl.df.right$y, "npc"),
                  just = "left", 
                  gp = gpar(col = highl.df.right$color, fill = highl.df.right$color, cex = 0.7 * as.numeric(highl.df.right$cex))
              )
            
            highl.df.top <- highl.df[highl.df$direction == "cover", ]
            if(nrow(highl.df.top) > 0)
              grid.text(
                  label = highl.df.top$genename, 
                  x = unit(highl.df.top$x, "npc"), 
                  y = unit(highl.df.top$y, "npc") + unit(0.013 * as.numeric(highl.df.top$cex), "npc"),
                  just = "bottom",
                  gp = gpar(col = highl.df.top$color, fill = highl.df.top$color, cex = 0.7 * as.numeric(highl.df.top$cex))
              )
          }
        }
      }
    }
  }
  
  ######### plot axes #########

  cat("Adding axes...\n")
  # y axis
  vp.yaxis <- viewport(x = unit(0.05, "npc"), width = 0.1, height = 0.8)
  upViewport()
  pushViewport(vp.yaxis)
  
  y.maxtick <- floor(max(gwasdat$logp))
  # at most 8 ticks
  ticks <- seq(1, y.maxtick, by = ceiling(y.maxtick/8))
  grid.text(
    label = ticks, 
    x = unit(0.5, "npc"), 
    y = ticks / max(gwasdat$logp),
    gp = gpar(col = "grey30", cex = 0.8)
  )
  
  grid.text(
    label = "-log10(p)", 
    x = unit(0.25, "npc"), 
    y = unit(0.5, "npc"),
    rot = 90, 
    gp = gpar(col = "grey30", fontface = "bold", cex = 0.8)
  )
  
  # x axis
  vp.axes <- viewport(y = unit(0.05, "npc"), width = 0.8, height = 0.1)
  upViewport()
  pushViewport(vp.axes)
  
  ticks.pos <- NULL
  for(i in chr.map$CHR.mapped) 
    ticks.pos = c(ticks.pos, mean(gwasdat[gwasdat$CHR.mapped == i, "x"]))
  
  grid.text(
    label = chr.map$CHR[odd(1:nrow(chr.map))], 
    x = ticks.pos[odd(1:nrow(chr.map))], 
    y = unit(0.5, "npc"),
    gp = gpar(col = "grey30", cex = 0.7)
  )
  
  grid.text(
    label = chr.map$CHR[even(1:nrow(chr.map))], 
    x = ticks.pos[even(1:nrow(chr.map))], 
    y = unit(0.6, "npc"),
    gp = gpar(col = "grey30", cex = 0.7)
    )
  
  grid.text(
    label = "Chromosome", 
    x = unit(0.5, "npc"), 
    y = unit(0.25, "npc"),
    gp = gpar(col = "grey30", fontface = "bold", cex = 0.8)
  )
  
  
  ######### plot title #########
  
  if(!is.null(plot.title)) {
    vp.title <- viewport(y = unit(0.975, "npc"), width = 0.1, height = 0.8)
    upViewport()
    pushViewport(vp.title)
    
    grid.text(
      label = plot.title, 
      x = unit(0.5, "npc"), 
      y = unit(0.5, "npc"),
      gp = gpar(col = "grey30", cex = 1)
    )
  }

  if(!is.null(toFile) && is.character(toFile) && length(toFile) == 1)
    dev.off()
  cat("Done.\n")
  return(NULL)
}
