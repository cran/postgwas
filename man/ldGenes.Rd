% File man/ldGenes.Rd
\name{ldGenes}
\alias{ldGenes}
\title{Evaluate LD between SNPs and their Neighboring Genes}

\description{
For all snps, get the next 4 genes up-and downstream together with their LD statistics. 
}

\usage{
ldGenes(snps, genes, run.parallel = FALSE, ...)
}

\arguments{
  \item{snps}{data frame. Has to contain columns 'SNP', 'CHR' and 'BP' for reference SNP IDs, chromosome and base position (both numeric!) }
  \item{genes}{data frame. Has to contain columns 'id' and numeric columns 'chromo', 'start' and 'end'. May contain additional columns (e.g. gene symbol) which are preserved in the return value. }
  \item{run.parallel}{boolean(1). CAUTION does only work on linux systems currently. Requires the package 'multicore'. Parallelizes LD calculation and genotype retrieval by forking the parent R process. Number of concurrent processes is set through parallel::detectCores(). Consumes x-fold memory according to the number of parallel processes. Parallelization confuses the order of status messages. }
  \item{...}{May be parameters used by a subfunction (e.g. for LD calculation, can specify ped- and map parameters for genotype retrieval, see documentation of the function 'getGenotypes') }
}

\details{
calculates the following LD measures between SNP and gene:
\itemize{
  \item{the snp from gene with best r^2}
  \item{the mean r^2 between query SNP and all snps within the gene and }
  \item{mean r^2 with the best quarter of snps from the gene.}
}
}

\value{
A data frame with the columns SNP, ld.max, ld.mean, ld.quart and all columns from the 'genes' parameter. For SNPs that lack a next gene at the requested level for any reason, row values are set to NA.
}

\seealso{
  \code{\link{getGenotypes}}
}
