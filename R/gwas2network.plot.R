gwas2network.plot <- function(
                   g = g,
                   layout = layout.fruchterman.reingold(g, area = vcount(g)^3),
                   vertex.to.edge.space = 1.3, 
                   device = pdf, 
                   ...
                 ) {

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
					 
  # there is a legacy package 'igraph0' that masks igraph and leads to errors
  # ensure that our igraph is first on search path 
  # and remember its current position for restore at function exit
  igraph.pos.orig <- which(search() == "package:igraph")
  detach("package:igraph", force = TRUE)
  suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = 2, quietly = TRUE)))
    
  layout[, 1] <- layout[, 1] - min(layout[, 1])
  layout[, 2] <- layout[, 2] - min(layout[, 2])
  
  papersize <- (max(layout))^0.45
  papersize <- if(papersize < 10) 10 else papersize
  papersize <- papersize * (if(!is.null(E(g)$label)) 1.25 else 1.1)  # more space when edge labels are plotted
  layout.normed <- layout / max(layout)
  layout.normed[is.na(layout.normed)] <- 0
  layout.normed <- as.data.frame(layout.normed) # we need that because subsetting a matrix to one row will yield a vector!

  args <- list(...)
  if(!is.null(list(...)$v.weight.max)) {
    vwm <- args$v.weight.max
    args$v.weight.max <- NULL
  } else {
    vwm <- max(V(g)$weight)    
  }
  V(g)$radius <- V(g)$weight / vwm
  
  
  if(any(c("width", "height", "units") %in% names(list(...))))
    stop("Arguments 'width', 'height' and 'units' cannot be specified in this function\n")
  if(identical(device, png) || identical(device, jpeg) || identical(device, tiff) || identical(device, bmp)) 
    do.call(device, c(width = papersize, height = papersize, units = "in", args))
  else
    do.call(device, c(width = papersize, height = papersize, args))


  # we create a viewport that leaves a border to the size of the full device
  frame <- viewport(width = unit(papersize - 3, "inches"), height = unit(papersize - 3, "inches"), name = "frame")
  pushViewport(frame)

  # plot edges first (are in background)
  if(ecount(g) > 0) {
    el <- na.omit(get.edgelist(g, names = F))   # vertex indices, eg V(g)$name[el], and layout index
    lapply(
      1:nrow(el), 
      function(idx) {
        # we define a start and end vertex for the edge
        x <- c(start = layout.normed[el[idx, 1], 1], end = layout.normed[el[idx, 2], 1])
        y <- c(start = layout.normed[el[idx, 1], 2], end = layout.normed[el[idx, 2], 2])
        
        # edges should be drawn from the vertex circle border
        # determine the point on the circle border where the edge should start
        # this is the intersection of edge and circle, where the edge is the hypothenusis or a triangle
    
        # radius of both circles (we add some space so that edge does not start directly at vertex border)
        r.start <- V(g)$radius[el[idx, 1]] * vertex.to.edge.space
        r.stop  <- V(g)$radius[el[idx, 2]] * vertex.to.edge.space
    
        # triangle
        adjacent <- abs(x["end"] - x["start"])
        opposite <- abs(y["end"] - y["start"])
        alpha    <- atan(opposite / adjacent)
        # for the triangle configuration with start vertex to the left and below end vertex
        x.adj <- cos(alpha)
        y.adj <- sin(alpha)

        if(x["start"] > x["end"]) 
          x.adj <- -x.adj
        if(y["start"] > y["end"])
          y.adj <- -y.adj
        
        x.just <- c(start = x.adj * r.start, end = -x.adj * r.stop)
        y.just <- c(start = y.adj * r.start, end = -y.adj * r.stop)
        grid.lines(
          x  = unit(x, "npc") + unit(x.just, "char"), 
          y  = unit(y, "npc") + unit(y.just, "char"),
          gp = gpar(
                 col = E(g)$color[idx]
               )
        )
        
        # draw labels
        if(!is.null(E(g)$label)) {
          grid.text(
            E(g)$label[idx], 
            x  = unit((x["start"] + x["end"]) / 2, "npc"), 
            y  = unit((y["start"] + y["end"]) / 2, "npc"),
            gp = gpar(
                   col = E(g)$color[idx], 
                   fontface = "italic", 
                   cex = 0.7  # make edge labels a bit smaller
                  )
          )
        }
      }
    )
  }
  
  # plot vertices
  if(vcount(g) > 0) {

    # category information is optional: all.cat can be empty, then category will be set to 'circle' globally
    # determine all categories used. We have to see that from column names, because
    # when for one category all vertices are a complete subset of another categories vertices, 
    # this will not be listed in category.main.
    all.cat <- list.vertex.attributes(g)[list.vertex.attributes(g) %in% c("circle", "square", "diamond")]
    
    if(length(all.cat) > 3)
      stop("Too many categories - cannot plot that.\n")
    
    if(length(all.cat) <= 0) {
      V(g)$category.main <- "circle"
    } else {
      # set unknown category values
      V(g)$category.main[is.na(V(g)$category.main)] <- all.cat[1]
      V(g)$category.main[!(V(g)$category.main %in% all.cat)] <- all.cat[1]
      for(curr.cat in all.cat)
        set.vertex.attribute(g, curr.cat, value = na.set(get.vertex.attribute(g, curr.cat)))
    }
    
    
    # plot the main category vertices
    v <- data.frame(
      cat = V(g)$category.main,
      color = V(g)$color,
      frame.color = "#00000000", # invisible frames, last two characters 00 = transparent, FF solid
      radius = V(g)$radius,
      marked = V(g)$marked
    )
    for(curr.cat in all.cat) 
      vertexplot(v = v[v$cat == curr.cat, ], layout = layout.normed[v$cat == curr.cat, ], type = curr.cat)
    
    
    # plot multi-category vertices on top of the main category vertices
    for(curr.cat in all.cat) {
      multi.idx <- get.vertex.attribute(g, curr.cat)
      # ignore NA multi category setting
      multi.idx[is.na(multi.idx)] <- FALSE
      multi.v <- data.frame(
        color = rep("#00000000", sum(multi.idx)), # invisible fill, last two characters 00 = transparent, FF solid
        frame.color = rep("white", sum(multi.idx)), 
#         sapply(
#           V(g)$color[multi.idx], 
#           function(clr) 
#             adjustcolor(clr, red.f = 0.9, blue.f = 0.9, green.f = 0.9, offset = c(0, 0, 0, 1))
#         ),
        radius = get.vertex.attribute(g, paste(curr.cat, "weight", sep = "."))[multi.idx] / vwm,
        marked = get.vertex.attribute(g, paste(curr.cat, "marked", sep = "."))[multi.idx]
      )
      vertexplot(v = multi.v, layout = layout.normed[multi.idx, ], type = curr.cat)
    }

    # plot vertex text (is always drawn once per gene, at the main category vertex)
    
    # set font (differs between multi category and single category vertices)
    V(g)$label.fontface <- "1" # plain
    for(curr.cat in all.cat)
      V(g)$label.fontface[get.vertex.attribute(g, curr.cat)] <- "4" # bold-italic

    grid.text(
      label = V(g)$label, 
      x = unit(layout.normed[, 1], "npc") + unit(V(g)$radius, "char"), 
      y = unit(layout.normed[, 2], "npc") + unit(V(g)$radius, "char"),
      just = "left",
      gp = gpar(
             col = V(g)$label.color, 
             fontface = as.numeric(as.vector(V(g)$label.fontface))
           )
    )
  }
 
  if(!identical(device, X11))
    dev.off()

  # restore search path
  suppressWarnings(detach("package:igraph", force = TRUE))
  suppressWarnings(suppressPackageStartupMessages(library(igraph, pos = igraph.pos.orig, quietly = TRUE)))

}




vertexplot <- function(v, layout, type = "circle") {

  if(!is.null(v) && nrow(v) > 0) {
    
    if(nrow(v) != nrow(layout))
      stop("Nonmatching size of vertex and layout data frames.\n")
    
    if(type == "circle") {
      grid.circle(
        x = layout[, 1], 
        y = layout[, 2], 
        r = unit(v$radius, "char"),
        gp = gpar(
          lex = 4 * v$radius,
          col = as.vector(v$frame.color),
          fill = as.vector(v$color)
        )
      )
    }

    if(type == "square") {
      grid.rect(
        x = layout[, 1], 
        y = layout[, 2], 
        width = unit(sqrt(pi * v$radius^2), "char"), 
        height = unit(sqrt(pi * v$radius^2), "char"), 
        gp = gpar(
          lex = 4 * v$radius,
          col = as.vector(v$frame.color),
          fill = as.vector(v$color)
        )
      )
    }

    if(type == "diamond") {
      grid.polygon(
        x = unit.c(   # x coords of the four corners of the diamond
              unit(layout[, 1], "npc"), 
              unit(layout[, 1], "npc") + unit(sqrt(pi * v$radius^2) /2, "char"), 
              unit(layout[, 1], "npc"), 
              unit(layout[, 1], "npc") - unit(sqrt(pi * v$radius^2) /2, "char")
            ), 
        y = unit.c(
              unit(layout[, 2], "npc") + unit(sqrt(pi * v$radius^2) /2, "char"), 
              unit(layout[, 2], "npc"), 
              unit(layout[, 2], "npc") - unit(sqrt(pi * v$radius^2) /2, "char"), 
              unit(layout[, 2], "npc")
            ),
        id = rep(1:nrow(v), 4), 
        gp = gpar(
          lex = 4 * v$radius,
          col = as.vector(v$frame.color),
          fill = as.vector(v$color)
          )
        )
    }

    # plot dual annot marks on vertex positions
    if(!is.null(v$marked) && sum(v$marked) > 0) {
      grid.points(
        x = layout[v$marked, 1], 
        y = layout[v$marked, 2], 
        size = unit(v$radius[v$marked], "char"),
        pch = 4,
        gp = gpar(
          lex = 2.5 * v$radius[v$marked], 
          col = "white"
        )
      )
    }
  }
  
}
