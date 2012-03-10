% File man/LDparallel.Rd
\name{LDparallel}
\alias{LDparallel}
\alias{LDparallel.genotype}
\title{Calculate Pairwise Linkage Disequilibrium with Multiple Threads}

\description{
Performs calculation of pairwise linkage disequilibrium for a number of SNPs on multiple CPU cores. 
}

\usage{
LDparallel(g, num.processes, ...)
LDparallel.genotype(g1, g2, ...)
}

\arguments{
  \item{g}{list. Has to contain genotype objects. Each list element represents a SNP with its associated genotypes.}
  \item{num.processes}{integer. The number of parallel processes to use. Defaults to the number of cores.}
  \item{g1}{object(genotype). Contains all genotypes for a SNP.}
  \item{g2}{object(genotype). Second SNP, like g1.}
  \item{...}{Arguments passed to the 'genetics' LD function}
}

\details{
This function offers thread-parallel calculation of linkage disequilibrium for a list of genotypes. 
The returned matrix of r-square correlation values is symmetric since LD(snp1, snp2) == LD(snp2, snp1). 
Uses a number of cores as determined by the multicore:::detectCores() and distributes (non-symmetric) SNP pairs evenly on them. 
The LDparallel.genotype is basically a copy of the original function LD.genotype from the package genetics but does not return P and chi^2. 
}

\value{
A matrix containing r-square values (numeric datatype). Row and column names match the SNP names provided in parameter g. 
}

\references{
  \url{http://cran.r-project.org/web/packages/genetics/index.html}
}

\seealso{
  \code{\link{LD}}
}

\examples{

# filenames of ped/map files with genotype data
pedmap.files <- system.file(
                  "extdata", 
                  c("example.map", "example.ped"), 
                  package = "postgwas"
                )
genotypes <- getGenotypes(
               snps = c("rs13240", "rs2479", "rs8411", "rs35681"), 
               map = pedmap.files[1], 
               ped = pedmap.files[2]
              )
if(any(installed.packages()[, "Package"] == "multicore")) {
  LDparallel(genotypes)
} else {
  cat("Package multicore not installed, cannot run example.\n")
}
}


