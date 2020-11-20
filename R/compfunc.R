#' Calculates the Distances Among Microbiome Functional Profiles
#'
#' @param abds string matrix indicating abundance matrix
#' @param rev logical variable indicating if the abundance matrix is reversed, default is FALSE
#' @param dist_type numeric indicating indicating which distance algorithm will be used, 0: HMS, 1: Cosine, 2: Euclidean, 3: Jensen-Shannon, 4: Bray-Curtis, default is 0
#' @param is_sim logical variable indicating whether the result matrix is similarity one or distance one, TRUE: similarity, FALSE: distance, default is FALSE
#'
#' @return the pairwise distance (similarity) matrix
#'
#' @export
#' @examples
#' data(abd_matrix)
#' compfunc(abd_matrix)
#' compfunc(abd_matrix, rev=0, dist_type=0, is_sim=0)
compfunc <- function(abds, rev = 0, dist_type = 0, is_sim = 0) { #rev: 0:norm, 1:rev 
	if (rev == 0) { #In R language, matrices are stored by columns, so the normal order matrices need to be transposed instead
		abds <- t(abds)
	}
	features <- rownames(abds)
	samples <- colnames(abds)
	data(gene)
	gene <- as.vector(gene)
	data(pw)
	pw <- as.vector(pw) 
        	result <- .Call('Multi_Comp_Table', gene, pw, features, samples, abds, dist_type, is_sim, PACKAGE = 'hrms')
	result <- matrix(unlist(result), nrow = length(samples))
	rownames(result) <- samples
	colnames(result) <- samples
	return (result)
}
