#include <Rcpp.h>
#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <sstream>
#include "src/hash.h"
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
vector<string> Load_db_file(SEXP file_name_s) {
    string file_name = as<string> (file_name_s);
    vector<string> results;
    ifstream infile(file_name.c_str(), ifstream::in);
    if (!infile){
        Rcpp::Rcout << "Error: Open Pathway file error : " << file_name << endl;
        vector<string> tmp;
        return tmp;
    }
    string buffer;
    while(getline(infile, buffer)){
        results.push_back(buffer);
    }
    infile.close();
    infile.clear();
    
    return results;
}

