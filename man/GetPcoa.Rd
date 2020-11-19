\name{GetPcoa}
\alias{GetPcoa}
\title{PCoA Implementation}
\usage{
GetPcoa(dist_matrix, k = 3)
}
\arguments{
\item{dist_matrix}{numeric indicating the distance matrix}

\item{k}{numeric variable indicating the number of dimension coordinate system (e.g. k = 2 or 3 space), default is 3}
}
\description{
A PCoA implementation optimized by multi-thread parallel computing to rapidly parse out the beta-diversity pattern for thousands of samples
}
\examples{
data(dist_matrix)
GetPcoa(dist_matrix)
GetPcoa(dist_matrix, k = 3)
}
