% File man/bm.remapSNPs.Rd
\name{bm.remapSnps}
\alias{bm.remapSnps}
\title{Annotate Position to SNPs}

\description{
  Annotate current biomart base position information to a file containing SNPs.
}

\usage{
bm.remapSnps(
  file, 
  config = biomartConfigs$hsapiens, 
  ds = bm.init.snps(config), 
  use.buffer = FALSE, 
  toFile = paste(file, ".remapped", sep = "")
)
}

\arguments{
  \item{file}{character(1). Path/filename, has to be readable by read.table (default settings) and must contain a header line with column names. Has to contain a 'SNP' column that is used for annotation (values should match the selected biomart SNP identifier field). Rows that contain NA values are removed.}
  \item{config}{A list containing values that are needed for biomaRt connection and data retrieval. See \code{\link{biomartConfigs}}.}
  \item{ds}{Optionally, a biomart connection object for re-use scenarios to prevent connection overhead.}
  \item{use.buffer}{Uses the variable postgwas.buffer.snps for annotation when it exists and seems to contain valid biomart data. When it does not exists, writes the biomart data used to that variable.}
  \item{toFile}{character(1). Path/filename, destination to write the remapped data to. Disabled when NULL.}
}

\value{
  A data frame containing the original file data plus columns 'CHR' and 'BP' with the biomart position annotation. Existing 'BP and CHR' columns will be renamed with the suffix .original.
}

\seealso{
  \code{\link{biomartConfigs}}
}
