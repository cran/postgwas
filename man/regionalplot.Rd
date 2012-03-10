% File man/regionalplot.Rd
\name{regionalplot}
\alias{regionalplot}
\title{Create Regional Association Plots of Multiple GWAS Datasets }

\description{
This is a plotting function for results of genomewide association studies (GWAS). It creates, for a list of SNPs (usually peak SNPs), regional plots for a given window size around these SNPs. These regional plots contain one or more p-value graphs (can plot multiple GWAS result datasets at once, e.g. for phenotypic comparisons), gene information, LD between SNPs from the GWAS and information about rare variants when available (e.g. resequencing studies). Most of the function parameters can be left at its defaults.
}

\usage{
regionalplot(
              snps, 
              gwas.resultfiles, 
              window.size = 1000000, 
              biomart.config = biomartConfigs$hsapiens,
              use.buffer = FALSE,
              plot.genes = TRUE, 
              draw.snpname = "auto",
              ld.options = list(
                             pop.id = 2, 
                             max.snps.per.window = 200, 
                             rsquare.min = 0.2, 
                             show.rsquare.text = FALSE
                           ),
              var.options = NULL,
              out.format = list(
                             file = "pdf", 
                             panels.per.page = "auto",
                             global.scale = "auto", 
                             paper.height = 11.7, 
                             paper.width = 8.3
                           ),
              max.logp = FALSE, 
              cores = 1,
              ytracks = ytracks.regionalplot, 
              panel.add = function(...) return(NULL)
)
}

\arguments{
  \item{snps}{data frame. Has to contain a column named "SNP" containing rs ids. For each SNP, a plot panel showing a region centered around the SNP is drawn.}
  \item{gwas.resultfiles}{character. The path/filename to one or more GWAS resultfiles containing the columns "SNP", "CHR", "BP", "P". They have to be readble by the \code{\link{read.table}} function with its default settings and header = TRUE. Make sure that each SNP occurs only once per file and is named according to dbSNP standards (rs-number). The 'CHR' column values should match ENSEMBL biomart conventions. When multiple files are supplied, they can contain different SNPs, but positions have to match.}
  \item{window.size}{integer(1). the total window size that defines the region to plot (i.e. window.size / 2 to each side of the query SNP)}
  \item{biomart.config}{list. Contains values that are needed for biomaRt connection and data retrieval (i.e. dataset name, attribute and filter names for SNPs and genes). This is organism specific. Normally, a predefined list from \code{\link{biomartConfigs}} is specified, e.g. biomartConfigs$mmusculus for the mouse organism. Available predefined lists are shown with \code{names(biomartConfigs)}. See also \code{\link{biomartConfigs}}.}
  \item{use.buffer}{boolean. When TRUE, several buffer variables in the global environment ("workspace") are used for SNP, gene, exon and LD data. When such variables (some or all) do not yet exist, data will be downloaded as usual and stored in these variables (e.g. for re-use in a second run or manual lookup). When they exist, the available data is used, skipping the usual download and calculation. A list of buffer variables and further explanation is given in the 'Details' section. Always remember to clear buffers when plotting different/new regions and use.buffer is TRUE.}
  \item{plot.genes}{boolean. Disables or enables gene/exon data track. When enabled (TRUE), genes are drawn and SNP and region positions will be mapped to match gene positions. Only genes with existing names are plotted.}
  \item{draw.snpname}{can be FALSE, "auto", or a custom data frame. When "auto", identifiers of all query SNP are drawn into the plots in a default manner. When FALSE, nothing is drawn. Can also be a data frame specifying a custom layout and text, see 'Details'.}
  \item{ld.options}{when NULL, LD track is omitted. Otherwise, takes a list with the following elements to draw an LD track in the plot:
    \describe{
      \item{pop.id}{a HapMap population identifier to retrieve genotyes for (can be NULL or omitted when ped and map arguments are given)}
      \item{ped}{character string containing a path/filename to a genotype pedigree file in linkage (PLINK) format. The file might contain more SNPs than actually needed. More details about the data format are described in the \code{\link{getGenotypes}} function. Can be omitted or NULL when pop.id is given.}
      \item{map}{path/filename to genotype mapfile in linkage format. Has to be in --map3 format (columns CHR, SNP, BP). Can be omitted or NULL when pop.id is given}
      \item{max.snps.per.window}{numeric, draw and calculate LD for at most max.snps.per.window SNPs in each region, evenly distributed. When regions overlap, ld pairs for the overlap will add up (you might want to use \code{\link{removeNeighborSnps}} to avoid these superfluous calculations). Such a restriction of LD calculation can be necessary in large regions for performance reasons. Default value is 200.}
      \item{rsquare.min}{the minimum LD value for an LD triangle to be visible in the plot. Remaining LD ist still calculated (and occasionally buffered). Default value is 0.2.}
      \item{show.rsquare.text}{annotates the r^2 value at the tip of each LD triangle. Default is FALSE.}
    }
  }
  \item{var.options}{when NULL, does not plot rare variants. Otherwise, a track with allele frequencies of additional variants e.g. from a resequencing project is drawn. Takes an argument list with the following elements: 
    \describe{
      \item{vcf}{a list with elements describing the input file(s) which have to be in tabix indexed variant call format (VCF with a data (.vcf.gz) and index (.tbi) file). Further, two files can be specified to draw a comparative histogram. Then, this argument is a list of two such lists. 
        \describe{
          \item{file}{character(1). This is full pathname of a tabix indexed vcf file to use for the variant display.}
           \item{remap.positions}{boolean. When FALSE, uses base positions from the vcf file and requires that they match with those used by the regionalplot function (which are either positions from the pvalue files when plot.genes = FALSE or the corresponding biomart positions otherwise). When TRUE, biomart positions are retrieved and used. This is done by extracting variants within a frame of 2 MB around the plotting region from the VCF file. When the shift between original and biomart positions exceeds this window, variant information might be incomplete. Unknown (e.g. de novo) SNPs are imputed using the offset between the nearest SNPs with annotated position.}
           \item{chrom.map}{a data frame. This is an option to correct different chromosome namings between the GWAS SNPs and variants from the vcf file. The data frame has to contain two columns named CHR and CHR.VCF, mapping chromosome names from the gwas data to chromosome names in the vcf file (use the \code{seqnamesTabix} function from the Rsamtools package to identify names used in the VCF file. Default is NULL (no mapping).}
        }
      }
      \item{vcf.id.prune}{a character vector. Each element of the vector is interpreted as a regular expression. Only variants that match one of these expressions in their ID field are displayed. Default is NULL (no pruning).}
      \item{vcf.info.prune}{a named character vector. Names have to be column names of the INFO field in the vcf file. Values are interpreted as a regular expression. Only variants that match one of these expressions in the appropriate INFO field are displayed. Default is NULL (no pruning).}
      \item{vcf.info.color}{a named character vector. Names have to be column names of the INFO field in the vcf file. Values are interpreted as regular expressions. For each of them, a color is chosen and matching variants in the appropriate INFO field will be displayed in the given color, on top of other variants. Colorization is applied in order of the expression list, thus variants that fall in several categories will be colored according to the last expression in the list. Default is NULL (no colorization).}
      \item{vcf.af.prune}{named numerical(1). The name has to be column names of the allele frequency field in the  INFO column of the vcf file (usually AF). Value is given in percent [0,50] for the minor allele. Only those variants will be displayed where the minor allele does not exceed this threshold. For two vcf files, the cap is applied to both individually. Default is AF=50 (display all)}
      \item{details}{when 1, draws a histogram plot of allele frequencies (or interval plot for two vcf files). When 2, adds thin vertical calibration lines for all colorized variants. When 3, adds variation ID and position for all colorized variants. Lists the original position in the first place and the remapped position second (for two vcf files, shows original position of the first file in first place). Default is 3.}
      \item{hist.alpha}{An optional alpha value for the histogram color. Default is 0.4.}
    }
  }
  \item{out.format}{a list. Specifies a file format and medium size for saving the plot. List elements are:
    \describe{
    \item{file}{character. Can be "pdf", "png", "jpeg" and "bmp" to create files in the corresponding format. Default is "pdf". When NULL, does not write a file.}
    \item{panels.per.page}{numeric value or "auto". The number of panels (i.e. regions) to plot per page. Default varies with the number of tracks. Default is "auto".}
    \item{global.scale}{numeric value or "auto". This is a parameter to scale the size of graphical elements with the paper size, and acts as a multiplier for the default sizes. It should scale approximately linear with the paper size and / or multiplied with the number of panels per page. Default is "auto".}
    \item{paper.width}{paper size in inches. Default is 8.3.}
    \item{paper.height}{paper size in inches. Default is 11.7.}
  }
  Can be NULL to not write a file. 
  }
  \item{max.logp}{numeric or FALSE. Threshold for the smallest p-value exponent (in -log, e.g. '12') to be shown on the y axis. When FALSE, y axis size corresponds to the minimum p-value found in gwas.resultfiles. Default is FALSE.}
  \item{cores}{integer. Specifies the allowed number of parallel processes. Requires the package 'multicore' to be installed when more than 1 process should be used for computation. Consumes a multiplicity of memory according to number parallel processes started. WORKS CURRENTLY ONLY ON LINUX SYSTEMS. Default is 1 (multithreading deactivated)}
  \item{ytracks}{a function. See 'examples' and \code{\link{ytracks.regionalplot}}.}
  \item{panel.add}{a function. See 'details'.}
}

\details{
This plotting tool has been specifically designed with use of lattice graphics to be able to produce a large number of plots at the same time (e.g. set a false discovery rate threshold and plot the best dozens or hundred SNPs of a study). When plotting a larger number of regions, you might want to select only SNPs at a certain distance between each other. The accompanying function 'removeNeighborSNPs' automatically selects those with the highest P value (if existent) within a specified window (does also work without supplying P value or position data). 
Via custom panel functions, it offers the opportunity to insert custom graphics or tracks, though it requires some effort to reserve the required space (see examples). It is further possible to manipulate the returned trellis object. What has additionally been proven quite useful is the ability to create pdf documents and search them for certain genes or SNPs.

Contents and construction of the data tracks in the plot: 
\describe{
  \item{LD track}{LD is calculated and visualized between SNPs contained in the region (i.e. occuring in any of the supplied gwas.resultfiles), given that genotype data exists. Genotypes are fetched by SNP rs ID, thus there is no need to consider position information (neither by HapMap nor ped/map retrieval). The LD is displayed using the commonly known red traingles, where the color intensity (using alpha transparency) scales linearily with the r^2 value of correlation between two SNPs. The color scale starts at the rsquare.min parameter an ends at 1. 
}
  \item{genes track}{Genes that have a HGNC symbol assigned (more precisely, when a value exists for the biomart.config$gene$attr$name field in biomart) are plotted as green arrows, where the direction corresponds to the strand (i.e. the arrow head corresponds to the 3' end'). When exon information is available, exons are drawn as black rectangles. This is determined by \code{any(is.na(unlist(biomart.config$exon)))}, thus when a single biomart.config$exon parameter is set to NA, exon information will neither be retrieved nor drawn or the existence of the buffer variable when use.buffer = TRUE. 
}
  \item{variant track}{When a vcf file with additional variants is available, its contents can be visualized in a separate track. Basically, a dotted line marking the position and an allele frequency histogram is plotted, and optionally the identifier and position number as text above (decided by the 'details' parameter). 
When two files are specified, comparative allele frequencies are drawn between the two datasets using arrows starting at the allele frequency of the second file and pointing towards the frequency value of the first file in a histogram-like style. Make sure that both files use the same reference allele (otherwise frequencies are not comparable). Variants that are only represented in one dataset will be set to frequency zero in the second dataset.}
  \item{pvalue track annotation}{For each SNP in the pvalue graph, a text annotation can be set. When draw.snpname = "auto", the identifier of the query SNP will be plotted. Custom layout specification can be done by supplying a data frame with columns 'snps' (determining a SNP ID as occuring in one of the pval files), 'text' (the text to be drawn), 'angle' (the connector line angle and label position in degree, where 0 degree is right to the SNPs p-value position, 90 is above and so on), 'length' (the length of the connector line, in multiples of the default), and 'cex' (the size of the text and line width. Currently, only a single value (the lowest) of the column is used for all labels.)}
  \item{custom tracks}{Definition of a custom track is shown in the examples section.}
}

Configuring biomart: 
SNP, gene and exon data is retrieved from biomart. A list that specifies the needed attribute and filter values in biomart has to be given as a parameter to the regionalplot function. For a large number of organisms, the configuration data is predefined. The available species can be listed by stating names(postgwas::biomartConfigs). All necessary fields in that list can be seen by stating \code{names(unlist(biomartConfigs$hsapiens))}. The values of these elements can be changed as long as they are valid in the corresponding biomart database. For example, setting \code{postgwas::biomartConfigs$gene$attr$name <- "uniprot_swissprot_accession"} will plot accession numbers instead of gene symbols. The \code{\link{listAttributes}} function from the biomart package can be used to display valid field names. For some organisms, exon information is not available by definition. When one of the 'exon' attribute elements is set to NA, exon information will be ignored. 

Using buffer data (custom data): 
Most functions in the whole package are able to store downloaded data in buffer variables in the global workspace. This is useful when the same plot (i.e. function call) is repeated several times with different graphical parameters, because the buffered data can be re-used then, which saves bandwith and computation time. In addition, it is possible to look up specific values in the buffer data or alter them, e.g. inserting a newly discovered or hypothetical gene. It is also possible to completely replace that data with custom or previously downloaded data when some online resources are not available. All buffer variables are independent from each other, thus it is possible to selectively delete or provide a single buffer variables.  
The following buffer variables exist: 'postgwas.buffer.snps.regionalplot' for all SNPs used (except for rare variants from the vcf file), 'postgwas.buffer.genes.regionalplot', 'postgwas.buffer.exons.regionalplot' and postgwas.buffer.ld.regionalplot. The contents and format of these variables can be explored by running \code{example(regionalplot)} which will created the buffer variables in the global environment. Remember to keep the buffer data current when plotting new/different regions and having use.buffer = TRUE.


When talking about (matching) base positions there is much to consider. These vary between genome assemblies, and the regionalplot function tries to handles them automatically whenever possible. In short, we have the following scenarios:  
\itemize{
  \item{Genes are not plotted: Base positions are directly used from the GWAS result files. Positions are checked to be the same within all GWAS result files. If this is not the case, the accompanying function bm.remapSnps can be used to create new files with synchronized base positions.}
  \item{Genes are plotted with automatic retrieval: Base positions are used from the GWAS result files to construct the regions, but are shifted afterwards to the base positions of the current ENSEMBL assembly (or the one given in the biomart.config argument). SNP positions that cannot be annotated this way are estimated by the offset of the closest annotated SNP. This allows the usage of SNPs in the GWAS that are not contained in dbSNP or biomart (except for the query SNPs itself).}
  \item{Genes are plotted with employment of custom buffer data: The approach is the same as above, but the position concordance is within the responsibility of the user by supplying correct buffer data.}
}

Drawing SNP names or custom text:
The data frame that has to be supplied to draw.snpname requires five columns named 'snps', 'text', 'angle', 'length' and 'cex'. The contents of these columns are self-explaining, 'snps' has to identify the snp at which to draw the text (by rs-id as in the gwas.resultfile), angle and length determine the connector line shape, and the 'text' argument is actually drawn. The character size factor 'cex' has to be identical for all SNPs (rows of the data frame). See the examples.

Finally, note that when a low SNP denistiy GWAS result file is used for plotting, the pvalue graph might be truncated at one or both sides of the graph. This is a natural effect because the graph line ends at the last SNP contained within the window border. 
}

\value{
A trellis object. As a side effect, a file is generated in the working directory named regionalplot1.* (numbered consecutively).
}


\seealso{
  \code{\link{biomartConfigs}}, \code{\link{bm.remapSnps}}, \code{\link{getGenotypes}}, \code{\link{removeNeighborSnps}}, \code{\link{ytracks.regionalplot}}
}


\examples{

# base parameter: 
# for each region, a centered SNP has to be be defined
# and  a GWAS resultfile for pvalue graph(s)

snps <- data.frame(SNP = c("rs4438106", "rs279565"))
gwas.resultfiles = system.file(
                     "extdata", 
                     c("example.gwas1", "example.gwas2"), 
                     package = "postgwas"
                   )

\dontrun{
# simple example showing genes and LD
# uses online database retrieval
regionalplot(
  snps = snps, 
  gwas.resultfiles = gwas.resultfiles
)

# this would save all downloaded data in buffer variables (given they do not yet exist)
# repeated execution (with use.buffer = TRUE) will then use the buffer data
regionalplot(
  snps = snps, 
  gwas.resultfiles = gwas.resultfiles,
  use.buffer = TRUE
)
}

\dontshow{
postgwas.buffer.snps <- read.table(
  system.file("extdata", "postgwas.buffer.snps", package = "postgwas"), header = TRUE
)
postgwas.buffer.genes.regionalplot <- read.table(
  system.file("extdata", "postgwas.buffer.genes.regionalplot", package = "postgwas"), header = TRUE
)
postgwas.buffer.exons.regionalplot <- read.table(
  system.file("extdata", "postgwas.buffer.exons.regionalplot", package = "postgwas"), header = TRUE
)
load(system.file("extdata", "postgwas.buffer.ld.regionalplot", package = "postgwas"))
}

# now use buffer data
# also use custom SNP annotation text
# LD is calculated from ped / map files 
# and stored in the buffer variable (postgwas.buffer.ld.regionalplot)
pedmap.files <- system.file(
                  "extdata", 
                  c("example.map", "example.ped"), 
                  package = "postgwas"
                )
invisible(regionalplot(
  snps = snps, 
  gwas.resultfiles = gwas.resultfiles, 
  use.buffer = TRUE, 
  draw.snpname = data.frame(
    snps = c("rs279565", "rs15750"), 
    text = c("rs279565", "this is some custom \n annotation text for SNP rs15750"),
    angle = c(20, 110),
    length = c(1, 2), 
    cex = c(0.8)
  ),
  ld.options = list(
    map = pedmap.files[1], 
    ped = pedmap.files[2]
  ),
  # do not write a pdf file - remove this line to do so
  out.format = NULL
))

# we can extract data or manipulate (customize) the buffer
# e.g. show LD rsquare for two SNPs from the first region
postgwas.buffer.ld.regionalplot[[1]]["rs13240", "rs8411"]
# also, remove the ARPC4-TTLL3 gene
postgwas.buffer.genes.regionalplot <- postgwas.buffer.genes.regionalplot[
                                        postgwas.buffer.genes.regionalplot$hgnc_symbol != "ARPC4-TTLL3", 
                                      ]


# it is possible to draw in the plot by defining a custom panel function
# e.g. draw a transcriptional cassette
invisible(regionalplot(
  snps = snps, 
  gwas.resultfiles = gwas.resultfiles, 
  use.buffer = TRUE, 
  panel.add = function(...) {
    panel.rect(
      xright = 46550000, 
      xleft = 46800000, 
      ytop = -1,
      ybottom = -4
    )
  }, 
  ld.options = list(
    map = pedmap.files[1], 
    ped = pedmap.files[2]
  ),
  # do not write a pdf file - remove this line to do so
  out.format = NULL
))

# same plots plus rare variants, filtering for de novo variants ('.')
reseq.file <- system.file("extdata", "example.vcf.gz", package = "postgwas")
regionalplot(
  snps = snps,
  gwas.resultfiles = gwas.resultfiles, 
  use.buffer = TRUE, 
  var.options = list(
    vcf = list(
            file = reseq.file, 
            remap.positions = FALSE
          ), 
    vcf.id.prune = c(":"), 
    vcf.info.color = list(EFF = "HIGH", EFF = "MODERATE", EFF = "LOW")
  ),
  # remove this argument or set file = "pdf" to plot to a file
  out.format = list(file = NULL, panels.per.page = 2, global.scale = 1.2)
)

# comparative frequencies with two different resequencing datasets
# also, exclude variants with MAF < 0.4 (is applied to both datasets)
reseq.file2 <- system.file("extdata", "example2.vcf.gz", package = "postgwas")
invisible(regionalplot(
  snps = snps,
  gwas.resultfiles = gwas.resultfiles, 
  use.buffer = TRUE, 
  var.options = list(
    vcf = list(
            list(
              file = reseq.file, 
              remap.positions = FALSE
            ), 
            list(
              file = reseq.file2, 
              remap.positions = FALSE,
              chrom.map = data.frame(CHR = 1:23, CHR.VCF = paste("chr", 1:23, sep = ""))
            )
          ),
    vcf.id.prune = c("rs622349", "rs1123078", "rs4965603", "rs878099", "rs7170610"), 
    vcf.info.color = list(AF = "."), 
    vcf.af.prune = c(AF=40)
  ),
  # do not write a pdf file - remove this line to do so
  out.format = NULL
))

# this adds a complete custom data track
# function to reserve space in the panel
ytracks <- function(...) {
  args <- list(...)
  tracks.original <- ytracks.regionalplot(...)
  rbind(
    tracks.original, 
    "cust" = list(
      ysize = args$ylim.upper * 0.5, 
      yspace = args$ylim.upper * 0.1, 
      ystart = -(sum(tracks.original$ysize, tracks.original$yspace) + args$ylim.upper * 0.1), 
      ystop = -(sum(tracks.original$ysize, tracks.original$yspace) + args$ylim.upper * 0.6), 
      name = "custom track"
    )
  )
}
# panel function that print some random data points in that track area
panel.add <- function(...) {
               # get the x and y bounds from the function arguments
               ystart <- list(...)$y.tracks["cust", "ystart"]
               ysize <- list(...)$y.tracks["cust", "ysize"]
               xstart <- list(...)$region.startbp
               xstop <- list(...)$region.endbp
               panel.text(
                 xstart + (xstop - xstart) / 2, 
                 ystart - ysize / 2, 
                 "custom track with random data", 
                 cex = list(...)$global.scale
               )
               panel.points(
                 x = runif(20, min = xstart, max = xstart + (xstop - xstart)), 
                 y = jitter(rep(ystart - ysize / 2, 20), factor = ysize / 2)
               )
             }
# we additionally apply custom out.format parameters 
# to scale the panels to a quadratic shape
regionalplot(
  snps = snps, 
  gwas.resultfiles = gwas.resultfiles, 
  use.buffer = TRUE, 
  ytracks  = ytracks, 
  panel.add = panel.add, 
  ld.options = list(
    map = pedmap.files[1], 
    ped = pedmap.files[2], 
    max.snps.per.window = 10, 
    show.rsquare.text = TRUE
  ),
  out.format = list(file = NULL, panels.per.page = 2, global.scale = 1.2)
)

\dontrun{
# an example using a different organism and smaller window size
# we have no LD database, but LD plot could be supplied with custom ped/mapfiles

snps <- data.frame(SNP = c("s15-470734", "s04-1462081"))
pval.file = system.file("extdata", "example.gwas3", package = "postgwas")

regionalplot(
  snps = snps,
  gwas.resultfiles = pval.file, 
  window.size = 15000, 
  ld.options = NULL, 
  biomart.config = biomartConfigs$scerevisiae
)
}

}
