% File man/getInteractions.Rd
\name{getInteractions}
\alias{getInteractions.ppi}
\alias{getInteractions.domains}
\alias{getInteractions.GO}
\title{Retrieve Biological Interaction Network Data}

\description{
 A collection of functions to download or calculate network data for biological interaction (or relatedness, co-occurence etc).
}

\usage{
getInteractions.ppi(
  filter.ids = c(), 
  additionalIdentifierTypes = c("OFFICIAL_SYMBOL", "ENTREZ_GENE", "ENSEMBL", "SWISSPROT"), 
  taxId = 9606,
  includeInteractors = FALSE, 
  sourceDatabaseList = c("BioGRID", "INTACT", "MINT", "DIP"), 
  toFile = "postgwas.interaction.download"
)
getInteractions.domains(
  filter.ids = "", 
  biomart.config = biomartConfigs$hsapiens,
  filter.type = biomartConfigs$hsapiens$gene$filter$name,
  ds = bm.init.genes(biomart.config),
  min.occurence = NULL, 
  max.occurence = length(filter.ids) / (logb(length(filter.ids), 3) +1),
  toFile = "postgwas.interaction.download"
)
getInteractions.GO(
  filter.ids, 
  GOpackagename = "org.Hs.eg.db", 
  ontology = "BP", 
  similarity.type = "max",
  toFile = "postgwas.interaction.download"
)
}

\arguments{
  \item{filter.ids}{vector. A list of identifiers that that are used to retrieve network data. The type of IDs is specified by additional parameters (can also be symbols). See details.}
  \item{additionalIdentifierTypes}{character. Specifies which type of IDs to consider when filtering with \code{filter.ids}. Can be multiple types (default). Is equivalent to the equally named biogrid REST parameter \url{http://wiki.thebiogrid.org/doku.php/biogridrest}. }
  \item{taxId}{numeric(1). The NCBI taxonomy identifier for the organism to download data from.}
  \item{includeInteractors}{boolean(1). When TRUE, download also genes that interact with genes listed in the argument \code{filter.ids}. Is equivalent to the equally named biogrid REST parameter \url{wiki.thebiogrid.org/doku.php/biogridrest}. This is sufficient to identify a shared interactor network. The additional biogrid parameter includeInteractorInteractions is disabled.}
  \item{sourceDatabaseList}{character vector. Considers interactions from the named databases. Is equivalent to the equally named biogrid REST parameter \url{wiki.thebiogrid.org/doku.php/biogridrest}. All elements are concatenated and passed to the REST query.}
  \item{biomart.config}{list. Specifies field names for biomart data retrieval. See \code{\link{biomartConfigs}}.}
  \item{filter.type}{character. A biomart filter field name, that specifies the type of IDs listed in the argument \code{filter.ids}.}
  \item{ds}{object. Biomart dataset to use. Is initialized automatically, but can also be explicitely passed to reduce runtime and network load in multiple calls.}
  \item{min.occurence}{numeric(1). Specifies the minimum total number of occurences for a domain to be considered for interaction. This can occasionally be useful to reduce the specifity of modules and complexity of the network., i.e. only domains with common functions will form modules. When NULL, this is deactivated (default). }
  \item{max.occurence}{numeric(1). Specifies the maximal total number of occurences for a domain to be considered for interaction. This can be used to avoid very common domains to create giant modules while increasing the functional specifity of modules. When NULL, this is deactivated. Is by default a log fraction of the number of vertices.}
  \item{GOpackagename }{character(1). Has to be the package name of an installed GO annotation package (e.g. org.Hs.eg.db, see also \url{http://www.bioconductor.org/packages/2.10/data/annotation/}). GO term information from this package will be used for similarity calculation. }
  \item{ontology}{character(1). Can be "BP", "MF" and "CC". See \code{\link[GOSim]{setOntology}}.}
  \item{similarity.type}{character(1). The GO similarity measure between genes, as specified by the 'similarity' argument in \code{\link[GOSim]{getGeneSim}} in the package \code{GOSim}}.
  \item{toFile}{character(1). A file/pathname where the downloaded interaction data is saved in tab-delimited format. When NULL, does not write a resultfile.}
}

\details{
\describe{
  \item{getInteractions.ppi}{This function uses the RESTful download feature of biogrid. Biogrid comprises protein interaction data from a large number of respositories and a braod range of organisms. All function arguments are basically directly passed to the fields of the REST query. Visit the website for more details and which IDs can be used: \url{http://wiki.thebiogrid.org/doku.php/biogridrest}. Note that by default, very common domains are excluded with the max.occurence argument.}
  \item{getInteractions.domains}{Uses biomart to annotate (superfamily) domains to each gene. The edges are formed by common domains between proteins. Uses gene symbols by default.}
  \item{getInteractions.GO}{Calculates relations between the query genes based on similarity of their GO term architecture (using the GOSim package). This is time-intensive and takes seveal hours for ~1000 genes. Note that interactions at arbitrary strength are returned (consider 'weight' column), these should manually be filtered for a minimum similarity of e.g. weight > 0.75. GO term retrieval depends on bioconductor annotation packages. The human (org.Hs.eg.db) package is installed by default. When switching to a different organism by installing the appropriate package (\url{http://www.bioconductor.org/packages/2.10/data/annotation/}), the type of IDs used may change. The IDs used are then specified by the annotation package. See the example section of \link{gwas2network} how to overcome such an ID problem exemplified with yeast.} 
}  
}

\value{
  A data frame with the first two columns containing the interactors plus additional columns.
  \describe{
    \item{getInteractions.ppi:}{returns columns "genename.x", "genename.y" with names according to biogrid, which are currently gene symbols.}
    \item{getInteractions.GO:}{returns columns "geneid.x", "geneid.y", "weight", the latter containing the similarity measure calculated by GOSim, and IDs / names according to the annotation package configuration, currently entrez gene IDs for the human package org.Hs.eg.db.}
    \item{getInteractions.domains:}{returns columns "genename.x", "genename.y", "geneid.x", "geneid.y", "domain.name", the latter specifying the shared domain that defines the edge, and IDs and names according to the biomartConfig used.}
  }
}

\seealso{
  \code{\link{biomartConfigs}}
}
