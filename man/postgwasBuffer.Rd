% File man/postgwasBuffer.Rd
\name{postgwasBuffer}
\alias{postgwasBuffer}

\title{Postgwas Buffer Variable Usage}

\description{
 This help page contains explanatory text on buffer variables in postgwas.
}

\details{
Buffer data is ALWAYS stored in a environment named postgwasBuffer by the postgwas functions. 
The variables in that environment can be listed with \code{ls("postgwasBuffer")}.
Nevertheless, buffer variables are accessed (used) in order of the search path \code{\link{search}}. 
Thus, when a buffer variable exists in the global workspace (.GlobalEnv), it will be used by postgwas functions. 
When a buffer variable exists in the global workspace (.GlobalEnv) and at the same time in the postgwasBuffer environment, the value from the global workspace will be used because it is listed first in the search path. 
To remove buffer variables, it is possible to either explicitely specify the environment where it is located, e.g. \code{rm("postgwas.buffer.genes", pos = "postgwasBuffer")} or repeatedly call \code{rm("postgwas.buffer.genes", inherits = TRUE)} until all definitions of this variable are removed. 
Of course it is also possible to start a new R session to remove all buffer data.
}

\seealso{
  \code{\link{search}}
}
