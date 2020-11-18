\name{CompFunc}
\alias{CompFunc}
\title{CompFunc}
\usage{
CompFunc(abds, rev = 0, dist_type = 0, is_sim = 0)
}
\arguments{
\item{abds}{string matrix indicating abundance matrix}

\item{rev}{logical variable indicating if the abundance matrix is reversed, default is FALSE}

\item{dist_type}{numeric indicating indicating which distance algorithm will be used, 0: HMS, 1: Cosine, 2: Euclidean, 3: Jensen-Shannon, 4: Bray-Curtis, default is 0}

\item{is_sim}{logical variable indicating whether the result matrix is similarity one or distance one, TRUE: similarity, FALSE: distance, default is FALSE}
}
\description{
Comprehensively calculates the distances among microbiome functional profiles by considering their multi-level metabolic pathway hierarchy among functional gene families
}
\examples{
data(abd_matrix)
CompFunc(abd_matrix)
CompFunc(abd_matrix, rev=0, dist_type=0, is_sim=0)
}
