#include "table_format.h"
#include "matrix.h"
#include <Rcpp.h>
#include <ctime>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
RcppExport SEXP get_pcoa(SEXP dis, SEXP k_) {
//	time_t now = time(0);
//	Rcpp::Rcout << "1: " << now << endl;
        Matrixx M(dis);
//	now = time(0);
//	Rcpp::Rcout << "2: " << now << endl;
        M.get_Deviation_Matrix();
//	now = time(0);
//      Rcpp::Rcout << "3: " << now << endl;
	int k = as<int> (k_);
        NumericMatrix pc = M.Get_PC_Matrix(k);
//	now = time(0);
//      Rcpp::Rcout << "4: " << now << endl;
	NumericVector per = M.get_percentage(k);
//	now = time(0);
//      Rcpp::Rcout << "5: " << now << endl;
        return (List::create(Named("pc") = pc, Named("per") = per));
}

