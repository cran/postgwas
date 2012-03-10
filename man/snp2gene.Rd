% File man/snp2gene.Rd
\name{snp2gene}
\alias{snp2gene}
\alias{snp2gene.prox}
\alias{snp2gene.LD}
\title{Annotate Genes to SNPs}

\description{
Adds columns with gene annotation to an input data frame of SNPs. Two modes are available: by proximity (genes covering the SNP and the next up- and downstream genes (including overlapping genes)) or, alternatively, genes that are in linkage equilibrium/disequilibrium with the query SNPs. Also adds biomart position annotation. 
}

\usage{
snp2gene.LD(
  snps, 
  biomart.config = biomartConfigs$hsapiens, 
  use.buffer = FALSE,
  by.genename = TRUE, 
  run.parallel = FALSE, 
  ...
) 
snp2gene.prox(
  snps, 
  level = -1, 
  biomart.config = biomartConfigs$hsapiens, 
  use.buffer = FALSE,
  by.genename = TRUE,
  print.format = FALSE
)
snp2gene(
  snps, 
  by.ld = FALSE, 
  biomart.config = biomartConfigs$hsapiens,
  use.buffer = FALSE, 
  by.genename = TRUE,
  level = -1, 
  print.format = FALSE, 
  run.parallel = FALSE, 
  ...
)
}

\arguments{
  \item{snps}{data frame. Has to contain a column named 'SNP' with rs-IDs (\url{http://www.ncbi.nlm.nih.gov/projects/SNP/}). May contain additional columns which are preserved. Will be made unique on the SNP column (duplicate SNPs discarded).}
  \item{biomart.config}{list. Contains values that are needed for biomaRt connection and data retrieval (i.e. dataset name, attribute and filter names for SNPs and genes). This is organism specific. Normally, a predefined list from \code{\link{biomartConfigs}} is specified, e.g. biomartConfigs$mmusculus for the mouse organism. Available predefined lists are shown with \code{names(biomartConfigs)}. See also \code{\link{biomartConfigs}}.}
  \item{by.ld}{boolean(1). When TRUE, calculates LD between the query SNP and all genes in a one-MB window plus the nearest gene up-and downstream. Genotype data can is retrieved from HapMap or local ped/map files, when these are supplied as additional parameters (see '...'). It is recommended to use ped/map data, which reduces webserver traffic and allows LD quantification in the actual cohort.}
  \item{level}{numeric(1). When 0, returns for each SNP the (single) closest gene. When > 0, returns the level-closest gene in each direction. When < 0, returns the closest gene in each direction, plus all genes that overlap with the closest in each direction.}
  \item{print.format}{boolean(1). When FALSE (default), a decomposed data frame is returned, with one row per SNP-gene association. Contains an additional column for the direction of the SNP-gene relation (upstream, downstream, ...). LD annotation always uses the decomposed return format. When TRUE, returns the input data frame (one row per SNP) with three additional columns for covering, up- and downstream genes.}
  \item{run.parallel}{bollean(1). CAUTION does only work on linux systems currently. Requires the package 'multicore'. Parallelizes LD calculation and genotype retrieval by forking the parent R process. Number of concurrent processes is set by invocation of parallel::detectCores(). Consumes x-fold memory according to the number of parallel processes. Parallelization might confuse the order of status messages.}
  \item{use.buffer}{boolean(1). When TRUE, buffers the downloaded data in two variables postgwas.buffer.snps and postgwas.buffer.genes. See also \link{postgwasBuffer}. When these variables already exists, data is not downloaded but used from there instead. This facilitates the possibility to re-use data, alter data or provide custom data. When setting use.buffer = T, make sure that existing buffer data is current!}
  \item{by.genename}{boolean(1). When TRUE, only genes that have a name are annotated. For genes with a single name but multiple IDs, a random representative ID is selected. When FALSE, uses all genes.}
  \item{...}{May be parameters used by a subfunction of snp2gene, namely "ped" and "map" for genotype retrieval, see documentation of the function \code{\link{getGenotypes}}. When using ped/mapfiles, all query SNPs have to occur in the mapfile.)}
}

\details{
  When gene annotation by proximity is requested (level < 0), we have to deal with cases of overlapping genes. This can be many (e.g. TRBV gene variants or dense gene clusters). 
  In such a case, all overlapping genes will be annotated (capped by cluster.maxsize) and form multiple rows in the decomposed format or, for the print format, occur as single string where gene names are separated by '|'. 
  LD calculation uses several measures to assess the amount of linkage disequilibrium between a SNP and a gene. In a preliminary step, the r-square correlation between the queried SNP and SNPs within the gene is calculated (HapMap CEU genotypes by default). Due to efficency constraints, only genes within a window of 1 MB are considered (includes partial overlap, when there are not gene within this range, the nearest is taken), and there is a cap of 100 SNPs per gene for LD calculation (when there are more, 100 SNPs are selected by position, evenly distributed over the gene). Output then contains 
  \itemize{
    \item{the best r-square value between the query SNP and all SNPs from the gene }
    \item{the mean r-square value between the query SNP and all SNPs from the gene }
  }
The whole annotation process makes use of the biomaRt database interface. All positions used and returned refer to the genome build version as provided by the according biomaRt retrival. The database and variables used can be specified by the user, e.g. it is possible to annotate accession numbers by altering the biomart.config, setting \code{myconfig$gene$attr$name <- 'uniprot_swissprot_accession'} 
}

\value{
  A data frame with all columns from the query data frame 'snps' plus columns 'BP' for the biomart position (renames an existing BP column to BP.original) and
  \itemize{
    \item{when by.ld = FALSE, columns 'geneid', 'genename', 'start', 'end' containing the gene id/name and position as in the biomart config and 'direction' giving the mode of annotation (cover, up- or downstream) }
    \item{when by.ld = TRUE, columns 'geneid', 'genename', 'start', 'end' containing the gene id/name and position as in the biomart config, and 'ld.mean' and 'ld.max' with the according rsquare value}
  }
  Each association of a SNP to a gene forms a row of the data frame. 
}

\seealso{
  \code{\link{biomartConfigs}}
}

\examples{

snps <- data.frame(SNP = c("rs172154", "rs759704"))
\dontshow{
## preload data for offline usage
postgwas.buffer.snps <- read.table(
  system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE
)
postgwas.buffer.genes <- read.table(
  system.file("extdata", "postgwas.buffer.genes", package = "postgwas"), header = TRUE
)

# offline LD annotation needs ped/map files e.g. from a genotyping project
pedmap.files <- system.file("extdata", c("example.map", "example.ped"), package = "postgwas")

## the buffer variables have been filled with proper data:
pedmap.files
head(postgwas.buffer.snps)
head(postgwas.buffer.genes)
}

\dontrun{
## simplest annotation, online
annot.prox <- snp2gene.prox(snps)
annot.ld <- snp2gene.LD(snps)

## extract a set of genes in LD e.g. for an enrichment analysis
geneset <- annot.ld[annot.ld$ld.max > 0.3, "genename"]

## the same with a different organism
snps <- data.frame(SNP = c("s01-10027", "s04-1469331", "s10-240843", "s15-474479"))
snp2gene.prox(snps, biomart.config = biomartConfigs$scerevisiae)

## how to use protein accession numbers
bconf <- biomartConfigs$hsapiens
bconf$gene$attr$name <- "uniprot_swissprot_accession"
snp2gene.prox(snps, biomart.config = bconf)

## this shows a list of further biomart attributes for human genes:
listAttributes(useMart("ensembl", "hsapiens_gene_ensembl"))
}

## annotation using buffer variables (runs offline when set properly)
snp2gene.LD(snps, use.buffer = TRUE, map = pedmap.files[1], ped = pedmap.files[2])
snp2gene.prox(snps, use.buffer = TRUE)

}

