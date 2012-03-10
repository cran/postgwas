getInteractions.ppi <- function(
                         filter.ids = c(), 
                         additionalIdentifierTypes = c("OFFICIAL_SYMBOL", "ENTREZ_GENE", "ENSEMBL", "SWISSPROT"), 
                         taxId = 9606,
                         includeInteractors = FALSE, 
                         sourceDatabaseList = c("BioGRID", "INTACT", "MINT", "DIP"), 
                         toFile = "postgwas.interaction.download"
                   ) {

  if(!is.null(filter.ids) && is.vector(filter.ids) && length(filter.ids) > 0) {
    filter.ids <- unique(as.vector(filter.ids))
    if(is.null(additionalIdentifierTypes))
      stop("Argument 'additionalIdentifierTypes' has to be set\n")
    if(is.null(includeInteractors))
      stop("Argument 'includeInteractors' has to be set\n")
    filterquery <- paste(
                     "geneList=", 
                     paste(filter.ids, collapse = "|"), "&", 
                     "additionalIdentifierTypes=",
                     paste(additionalIdentifierTypes, collapse = "|"), "&",
                     "includeInteractorInteractions=FALSE", "&", 
                     "includeInteractors=", includeInteractors,
                     sep = ""
                   )
    # for many filter ids, we have to divide the queries (URL length cap)
    # when we do that, including interactor interactions is invalid (this feature is disabled anyways)
    if(length(filter.ids) > 500) {
      res <- list()
      filter.ids.part <- lapply(
                           seq(1, length(filter.ids), by = 500),  # the split points
                           function(split)
                             as.vector(na.omit(filter.ids[split:(split+499)]))
                         )
      res <- lapply(
               filter.ids.part, 
               function(part) 
                 getInteractions.ppi(
                   filter.ids = part,
                   taxId = taxId,
                   additionalIdentifierTypes = additionalIdentifierTypes, 
                   includeInteractors = includeInteractors, 
                   sourceDatabaseList = sourceDatabaseList, 
                   toFile = NULL
                 ) 
             )
      return(unique(list2df(res)))
    }
  } else {
    warning("Argument 'filter.ids' is not a vector with length > 0. Downloading all interactions.\n")
    filterquery = NULL
  }

  # determine number of results to retrieve - there is a cap of 10000 per query
  cat("Determining the number of interactions to download...\n")
  url.rows <- url(paste(
                  "http://webservice.thebiogrid.org/resources/interactions?", 
                  "selfInteractionsExcluded=true&format=count&", 
                  filterquery, "&", 
                  "taxId=",
                  taxId,
                  sep = ""
                ))

  rows <- as.numeric(readLines(url.rows, warn = FALSE))
  close(url.rows)
  
  if(rows == 0)
    return(data.frame(genename.x = "", genename.y = "")[-1, ])
  
  queries.data <- paste(
                    "http://webservice.thebiogrid.org/resources/interactions?",
                    "selfInteractionsExcluded=true&format=tab1&",
                    filterquery, "&",
                    "sourceDatabaseList=",
                    paste(sourceDatabaseList, collapse = "|"), 
                    "&taxId=",
                    taxId, 
                    "&start=",
                    gsub(" ", "", format(0:(rows/10000) * 10000, scientific = FALSE)), 
                    sep = ""
                  )
  queries.data[1] <- paste(queries.data[1], "includeHeader=true", sep = "&")
  urls.data <- lapply(queries.data, url)

  # paste together the repeated retrieves and download
  cat(paste("Downloading", rows, "interactions"))
  res <- lapply(
           urls.data, 
           function(url) {
             cat(".")
             dl <- list2df(strsplit(readLines(url, warn = FALSE), split = "\t", fixed = TRUE))
             close(url)
             return(dl)
           }
         )
  cat("\n")
  res <- as.data.frame(list2df(res))
  colnames(res) <- as.matrix(res)[1, ]
  res <- res[-1, ]  #remove header

  if(!is.null(toFile) && toFile != "")
    write.table(res, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)

  return(
    data.frame(
      genename.x = res[, "OFFICIAL_SYMBOL_A"], 
      genename.y = res[, "OFFICIAL_SYMBOL_B"], 
      stringsAsFactors = FALSE
    )
  )
}




# returns a data frame "genename.x", "genename.y", "geneid.x", "geneid.y", "domain.name"
getInteractions.domains <- function(
                             filter.ids = "", 
                             biomart.config = biomartConfigs$hsapiens,
                             filter.type = biomartConfigs$hsapiens$gene$filter$name,
                             ds = bm.init.genes(biomart.config),
                             min.occurence = NULL, 
                             max.occurence = length(filter.ids) / (logb(length(filter.ids), 3) +1),
                             toFile = "postgwas.interaction.download"
                           ) {
  
  if(is.null(filter.ids) || !is.vector(filter.ids))
    stop("Argument 'filter' has to be a non-NULL vector.\n")
  
  filter.ids <- as.vector(filter.ids)
  
  if(is.null(filter.type))
    stop("Argument 'filter' has to be set.\n")
  
  dom <- bm.genes.domains(
           config = biomart.config, 
           ds = ds, 
           filter.name = filter.type, 
           filter.val = filter.ids
         )

  dom.net <- merge(dom, dom, by = "domain.name")
  
  if(!is.null(min.occurence)) {
    keep <- table(dom.net$domain.name)
    keep <- names(keep)[keep > min.occurence]
    dom.net <- dom.net[dom.net$domain.name %in% keep, ]  
  }
  
  if(!is.null(max.occurence)) {
    keep <- table(dom.net$domain.name)
    keep <- names(keep)[keep < max.occurence]
    dom.net <- dom.net[dom.net$domain.name %in% keep, ]  
  }
  
  dom.net <- dom.net[, c(3,5,2,4,1)]
  
  if(!is.null(toFile) && toFile != "")
    write.table(dom.net, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(dom.net)
}




# returns a data frame "geneid.x", "geneid.y", "weight"
getInteractions.GO <- function(
                             filter.ids, 
                             GOpackagename = "org.Hs.eg.db", 
                             ontology = "BP", 
                             similarity.type = "max",
                             toFile = "postgwas.interaction.download"
                           ) {
  
  if(!"GOSim" %in% names(installed.packages()[, 'Package'])) {
    stop("This functionality requires the 'GOSim' package to be installed.\n")
  } else {
    suppressPackageStartupMessages(require(GOSim, quietly = TRUE))
  }

  if(is.null(filter.ids) || !is.vector(filter.ids) || length(filter.ids) <= 0)
    stop("Argument 'filter.ids' has to be set to a nonempty vector.\n")

  # load proper GO organism data
  success <- do.call(library, list(package = GOpackagename, logical.return = TRUE))
  if(!success)
	stop("The GeneOnotology annotation package supplied in argument 'GOpackagename' is not installed or cannot be loaded.\n")
  if(!exists("GOSimEnv")) GOSim:::initialize()
  # reinit GOSim when intitalized with a different GO annotation package
  org.name.new <- get(grep("ORGANISM", ls(paste("package", GOpackagename, sep = ":")), value = TRUE)[1])
  org.dat.new <- get(grep("GO", ls(paste("package", GOpackagename, sep = ":")), value = TRUE)[1])
  if(sub("human", "Homo sapiens", GOSimEnv$organism) == org.name.new) {
	# everything is already preloaded
  } else {
	setEvidenceLevel(organism = org.name.new, gomap = org.dat.new)
	calcICs(path.package("GOSim"))
	setOntology(ontology, DIR = path.package("GOSim"))
  }

  cat("Starting gene similyrity calculation with GOSim (this can take a while)...\n")
  simmatrix <- getGeneSim(as.vector(filter.ids), similarity = similarity.type)
  
  res <- list2df(lapply(
    colnames(simmatrix), 
    function(colname) 
      data.frame(
        geneid.x = colname, 
        geneid.y = rownames(simmatrix), 
        weight = simmatrix[, colname]
        )
    ))
  
  if(!is.null(toFile) && toFile != "")
    write.table(res, toFile, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  return(res)
  
}

