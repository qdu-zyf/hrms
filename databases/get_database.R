library(hrms)
library(Rcpp)
library(devtools)

sourceCpp("get_database.cpp")

gene <- Load_db_file("KO/ko_id.tab")
use_data(gene, hrms)
pw <- Load_db_file("KO/ko_pw.tab")
use_data(pw, hrms)
