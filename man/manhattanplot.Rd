% File man/manhattanplot.Rd
\name{manhattanplot}
\alias{manhattanplot}
\title{Generate a Manhattan Plot from GWAS Data}

\description{
Produces a manhattan plot with highlighting and annotation features.
}

\usage{
manhattanplot(
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
)
}

\arguments{
  \item{gwas.resultfile}{character(1). Contains the path/filename of a GWAS resultfile. The file must contain the columns "CHR", "BP", "SNP", "P" and an approriate header line and be readable by the \code{link{read.table}} function.}
  \item{highlight.logp}{integer. The value defines a threshold for highlighting SNPs with -log(p) < highlight.logp. See also the arguments below. Can be a vector of multiple thresholds. }
  \item{highlight.win}{integer. For a SNP above the highlighting threshold, all SNPS lying within highlight.win basepair distance are also highlighted (applied to both sides). This argument is further used to determine a representative (i.e. lead) SNP of multiple proximate SNPs that all hit the highlighting threshold. The gene annotation is only plotted at the lead SNP of a block. The length of this argument should match the length of highlight.logp or is recycled if necessary.}
  \item{highlight.color}{character. Can be set to a user-defined color as listed in colors(). The length of this argument should match the length of highlight.logp or is recycled if necessary. }
  \item{highlight.text}{character(1). When 'SNP', draws SNP identifiers above the peak SNP of each highlighted locus. When "genes" runs snp2gene to annotate the closest genes in each direction. When NULL, draws no annotation at all.}
  \item{highlight.cex}{numeric. This is applied as cex argument (see \code{\link{par}}) to the gene labels, and in similar form to the plotted point symbol size of highlighted SNPs. The length of this argument should match the length of highlight.logp or is recycled if necessary.}
  \item{max.y}{numeric. The maximum -log(p) value to display. This can be used when very small p's occur (remaining data is then displayed very compressed.)}
  \item{reduce.dataset}{boolean(1). When TRUE, randomly removes SNPs to reduce data load, where the probability for a SNP to be removed is the squared p-value. SNPs with p < 1e-04 are untouched.}
  \item{plot.title}{character(1). Title text of the plot.}
 \item{biomart.config}{list. A biomart config list as described in \code{\link{biomartConfigs}}. This is only cosidered when highlight.text = "genes" and directly passed to the snp2gene function for gene annotations.}
  \item{use.buffer}{boolean(1). As in \code{\link{snp2gene}}.}
  \item{toFile}{character. A file/pathname to write the plot to (in pdf format). When NULL, plots to the current device (usually the display).}
}

\details{
Plots a grayscale-separated manhattan plot, with SNPs exceeding a user-defined threshold being highlighted, colorized and annotated with surrounding genes (or custom text). 
If a single locus gets genes annotated at several SNPs belonging to the same locus, extend the highlight.win size to unify them. 
}

\value{
 NULL, plots to the current device as a side effect.
}

\seealso{
  \code{\link{snp2gene}}
}

\examples{

\dontshow{
# load snps annotation into biomart buffer (runs offline)
postgwas.buffer.snps <- read.table(system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE)
postgwas.buffer.genes <- read.table(system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE)
}
  
  # our file to plot
  gwas.resultfile = system.file("extdata", "example.gwas1", package = "postgwas")
  
  # basic plot
  manhattanplot(
    gwas.resultfile, 
    highlight.logp = c(3, 4), 
    use.buffer = TRUE, 
    plot.title = "truncated genomewide example dataset",
    toFile = NULL
  )
  dev.off()


\dontrun{
  # for yeast
  # we use our truncated synthetical dataset..
  gwas.resultfile = system.file("extdata", "example.gwas3", package = "postgwas")
  manhattanplot(
    gwas.resultfile, 
    highlight.win = 500, 
    highlight.logp = c(3, 4), 
    highlight.cex = 1.3, 
    biomart.config = biomartConfigs$scerevisiae, 
    use.buffer = TRUE
  )
}


  # we see that a locus appears twice (annotating the same genes)
  # increase window size to combine these loci
  # also, set alternative colors and size
  manhattanplot(
    gwas.resultfile, 
    highlight.logp = c(3, 4), 
    highlight.win = 130000,
    highlight.color = c("darkgreen", "red"),
    use.buffer = TRUE,
    plot.title = "truncated genomewide example dataset", 
    toFile = NULL
  )
}


