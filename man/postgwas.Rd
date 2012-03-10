% File man/postgwas.Rd
\name{postgwas}
\alias{postgwas}
\title{Default GWAS Visualization, Annotation and Network}

\description{
  Run a default set of postgwas functions for a quick analysis of your GWAS dataset.
}

\usage{
postgwas(
  gwas.resultfile, 
  suggestive.p = 1e-06, 
  genomewide.p = 5e-08,
  biomart.config = biomartConfigs$hsapiens, 
  GOpackagename = "org.Hs.eg.db", 
  ped = NULL, 
  map = NULL
)
}

\arguments{
  \item{gwas.resultfile}{character(1). The path/filename to a GWAS resultfile containing the columns "SNP", "CHR", "BP", "P". It has to be readble by \code{\link{read.table}} with its default settings and header = TRUE. Make sure that each SNP occurs only once and is named according to dbSNP standards (rs-number). }
  \item{suggestive.p}{numeric(1). This p-value threshold selects SNPs that are considered for network analysis and regional plots.}
  \item{genomewide.p}{numeric(1). Currently, this only determines the height of the genomewide-line in the manhattanplot.}
  \item{biomart.config}{list. A configuration list for biomart data retrieval as described in \code{\link{biomartConfigs}}. Has to be changed when the species in the GWAS is not human.}
  \item{GOpackagename }{character(1). When GO enrichment in the network analysis should be used, this has to be the package name of an installed GO annotation package. It is preconfigured for human GO dataset, further species can be installed from \url{http://www.bioconductor.org/packages/2.10/data/annotation/} and the according package name passed to this argument when using non-human species. Can otherwise be NULL to avoid vertex colorization.}
  \item{ped}{character(1). Path/filename to a pedfile in linkage format as produced by plink. Should contain genotypes of the GWAS study. When specified, considers LD information for snp to gene assignment an produces and LD track in the regionalplots.}
  \item{map}{character(1). Path/filename to a mapfile in linkage format as produced by plink (as above)}
}

\details{
 Applies the most important functions of postgwas sequentially to the supplied gwas dataset. This will produce several output files (see value section). It further makes use of several buffer variables in the global workspace that contain downloaded data for re-use (postgwas.buffer.*). Preexisting buffers are cleared from the workspace. See the help pages of the individual functions to get to know how they use buffer variables.
}

\value{
  Does not have a proper return value, but creates files as a side effect. 
  The following files are created / replaced: postgwas.snp2gene.csv, gwas2network.pdf, gwas2network.png, gwasGOenrich.csv, gwas2networkParams.txt, gwas2networkGraphVertices.csv, gwas2networkGraphEdges.csv
  The following files are consecutively numbered: regionalplot*.pdf, gwas2networkCommunities*.pdf, gwas2networkCommunities*.png
}

\seealso{
  \code{\link{snp2gene}}, \code{\link{regionalplot}}, \code{\link{gwas2network}}, \code{\link{removeNeighborSnps}}
}

\examples{
### this example requires online access and the optional packages 'igraph' and 'GOSim' to be installed
\dontrun{
  gwas.resultfile <- system.file(
                       "extdata", 
                       "example.gwas1", 
                       package = "postgwas"
                      )
  pedmap.files <- system.file(
                    "extdata", 
                    c("example.map", "example.ped"), 
                    package = "postgwas"
                  )
  postgwas(
    gwas.resultfile, 
    suggestive.p = 5e-02, 
    map = pedmap.files[1], 
    ped = pedmap.files[2]
  )
}
}