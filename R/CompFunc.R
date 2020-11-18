CompFunc <- function(abds, rev = 0, dist_type = 0, is_sim = 0) { #rev: 0:norm, 1:rev 
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

