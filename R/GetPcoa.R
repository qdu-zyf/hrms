GetPcoa <- function(dist_matrix, k = 3){
	res = .Call('get_pcoa', dist_matrix, k, PACKAGE = "hrms")
	return (res)
}

