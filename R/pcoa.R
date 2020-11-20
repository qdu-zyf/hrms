#' PCoA Implementation
#'
#' @param dist_matrix numeric matrix indicating the pairwise distance matrix
#' @param k numeric variable indicating the number of dimension coordinate system (e.g. k = 2 or 3 space), default is 3
#'
#' @return the PCoA coordinates matrix
#'
#' @export
#' @examples
#' data(dist_matrix)
#' pcoa(dist_matrix)
#' pcoa(dist_matrix, k = 3)
pcoa <- function(dist_matrix, k = 3){
	res = .Call('get_pcoa', dist_matrix, k, PACKAGE = "hrms")
	return (res)
}

