% File man/getGenotypes.Rd
\name{getGenotypes}
\alias{getGenotypes}
\title{Retrieve Genotypes}

\description{
Retrieves genotype data for the given SNPs. Can do online HapMap queries or read ped/map genotype files.
}

\usage{
getGenotypes(
  snps, 
  pop.id = 2, 
  ped = NULL, 
  map = NULL, 
  remove.homozygous = FALSE, 
  as.geno.objects = TRUE
)
}

\arguments{
  \item{snps}{character. A vector with SNP rs-ids to get genotypes for.}
  \item{pop.id}{numeric(1). This sets the Hapmap population ID for which genotypes should be retrieved (2 = CEU). When ped/map arguments are set, they will override the HapMart retrieval and use the given files. }
  \item{ped}{character(1). Path/filename to genotype files in PLINK format (can be omitted when pop.id is given). Reads various formats (numeric or character nucleotides, delimiters can be either space only or tab with spaces between alleles, ...). May contain genotypes for more SNPs than needed, and arbitrary additional columns ahead of genotype data (i.e. sex, family data, ...). }
  \item{map}{character(1). Reads linkage map format as well as --map3 (CHR, SNP, BP) format. }
  \item{remove.homozygous}{boolean(1). Returns only genotypes of SNPs that have more than one allele state. }
  \item{as.geno.objects}{boolean(1). Requires package 'genetics'. When TRUE, returns a data frame of genotype objects instead of genotype character codes. }
}

\value{
When return.geno.objects = TRUE, returns a data frame of genotype objects, one column for each SNP (see function 'genotype' in package 'genetics')
Otherwise, returns a matrix with genotypes, one column for each SNP. Genotype format is 2-character nucleotides, or NA, (e.g. 'AC').
May return NULL when genotypes are not available. May return an error String when HapMart is not available or the query is misconfigured. 
}

