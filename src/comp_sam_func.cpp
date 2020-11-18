// Updated at Aug 13, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format
#include <iostream>

#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "utility.h"
#include "version.h"
#include "comp_sam_func.h"
//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

char Ref_db;

string Listfilename;
string Listprefix;

string Queryfile1;
string Queryfile2;

string Tablefilename;
string Outfilename;

//int Coren = 0;
int Coren = sysconf(_SC_NPROCESSORS_CONF);

int Dist_matrix = 0; //0: cos 1: eu;

bool Is_sim; //true: sim, false: dist;

int Mode = 0; //0: single, 1: multi, 2: multi_table
bool Reversed_table = false;
bool Is_heatmap;
int Cluster = 2;

/*
void test(){
	 
	 _Comp_Tree_Func comp_tree_func;
	 //comp_tree_func.Debug_output("debug.out");
	 }
*/


void check(int & Dist_matrix) {
	if ((Dist_matrix > 4) || (Dist_matrix < 0)){
                Rcpp::Rcout << "Warning: Distance matrix must be 0-4, change to default (0)" << endl;
                Dist_matrix = 0;
        }
}

/*
vector<float>  manage_matrix(int n, vector<float> & sim_matrix, int is_sim) {//exchange the 1-dim dist matrix to 2-dim one
	vector<float> result(n * n);
	for(int i = 0; i < n; i ++){

		if(is_sim) result[i * n + i] = 1.0;
		else result[i * n + i] = 0.0;

		for (int j = i+1; j < n; j ++){
                        unsigned long p = i * (long) n + j - (1 + i + 1) * (i + 1) / 2;
                        if (is_sim){
				result[i * n + j] = result[j * n + i] = sim_matrix[p];
                        }	
                        else {
				result[i * n + j] = result[j * n + i] = 1 - sim_matrix[p];
                        }
		}

	}
	return result;
}
*/

void Output_Matrixx(const char * outfilename, int n, vector<vector<float> > & sim_matrix, vector <string> & sam_name){
	 
	 FILE * outfile = fopen(outfilename, "w");

	 /*
	 if (outfile == NULL){
		Rcpp::Rcout << "Error: Cannot open output file : " << outfilename << endl;
		return; 
	 }
	 */

	 //Label
	 fprintf(outfile, "SampleID"); ////
	 for(int i = 0; i < n; i ++)
			 fprintf(outfile, "\t%s", sam_name[i].c_str());
	 fprintf(outfile, "\n");
	 
	 for(int i = 0; i < n; i ++){
		fprintf(outfile, "%s", sam_name[i].c_str());
		
		for (int j = 0; j < n; j ++){				
			fprintf(outfile, "\t%f", sim_matrix[i][j]);
		}
	
		fprintf(outfile, "\n");
	}
			 
	fclose(outfile);
}

/*
// [[Rcpp::export]]	 
RcppExport SEXP Single_Comp(SEXP Queryfile1, SEXP Queryfile2, SEXP is_sim_s, SEXP Dist_matrix_s){
	 
	//_Comp_Tree_Func comp_tree_func(Ref_db);
	_Comp_Tree_Func comp_tree_func;
	
	float * abd_1 = new float [comp_tree_func.Get_GeneN()];
	float * abd_2 = new float [comp_tree_func.Get_GeneN()];
	
	//Rcpp::Rcout << comp_tree_func.Load_Gene_Count(Queryfile1.c_str(), abd_1) << " loaded" << endl;
	//Rcpp::Rcout << comp_tree_func.Load_Gene_Count(Queryfile2.c_str(), abd_2) << " loaded" << endl;
	string Queryfile1_s = as<string> (Queryfile1);
	string Queryfile2_s = as<string> (Queryfile2);
	comp_tree_func.Load_Gene_Count(Queryfile1_s.c_str(), abd_1);
	comp_tree_func.Load_Gene_Count(Queryfile2_s.c_str(), abd_2);
	
	int Dist_matrix = as<int> (Dist_matrix_s);
	check(Dist_matrix);
	double sim = comp_tree_func.Calc_sim(abd_1, abd_2, Dist_matrix);

	int is_sim = as<int> (is_sim_s);
        if (is_sim) return (wrap(sim));
	else return (wrap(1.0-sim));
}

vector<vector<float> > multi_comp(const char * Listfilename, const char * Listprefix, vector<string> &samples, int Dist_matrix, int is_sim){
	
	//_Comp_Tree_Func comp_tree_func(Ref_db);
	_Comp_Tree_Func comp_tree_func;
		 
	 //load list
	//vector <string> sam_name;
	vector <string> file_list;
	
	//int file_count = Load_List(Listfilename.c_str(), file_list, sam_name, Listprefix);
	int file_count = Load_List(Listfilename, file_list, samples, Listprefix);
		
	//load abd
	float **Abd = new float * [file_count];
	for (int i = 0; i < file_count; i ++){
		Abd[i] = new float [comp_tree_func.Get_GeneN()];
		//Rcpp::Rcout << comp_tree_func.Load_Gene_Count(file_list[i].c_str(), Abd[i]) << " KOs in file " << i + 1 << endl;
		comp_tree_func.Load_Gene_Count(file_list[i].c_str(), Abd[i]);
		}
	
	//make order
	vector <int> order_m;
	vector <int> order_n;
	unsigned long iter = 0;
	for (int i = 0; i < file_count - 1; i ++)
		for (int j = i + 1; j < file_count; j ++){			
			order_m.push_back(i);
			order_n.push_back(j);
			iter ++;
			}
		
	vector <float>  sim_matrix;
	for (unsigned long i = 0; i < iter; i ++)
		sim_matrix.push_back(0);
		
	//openmp	
	omp_set_num_threads(Coren);
	//Rcpp::Rcout << "Coren = " << Coren << endl;
	//omp_set_num_threads(8);
	check(Dist_matrix);

	#pragma omp parallel for
	for (unsigned long i = 0; i < iter; i ++){
		unsigned long m = order_m[i];
		unsigned long n = order_n[i];
		unsigned long p = m * (long) file_count + n - (1 + m + 1) * (m + 1) / 2;

		sim_matrix[p] = comp_tree_func.Calc_sim(Abd[m], Abd[n], Dist_matrix);
	}

	vector<vector<float> > results = manage_matrix(file_count, sim_matrix, is_sim);

	return results;
}

//[[Rcpp::export]]
RcppExport SEXP Multi_Comp(SEXP Listfilename_s, SEXP Listprefix_s, SEXP Dist_matrix_s, SEXP is_sim_s){
	
	string Listfilename = as<string> (Listfilename_s);
	string Listprefix = as<string> (Listprefix_s);
	int is_sim = as<int> (is_sim_s);
	int Dist_matrix = as<int> (Dist_matrix_s);
	vector<string> samples;
	vector<vector<float> > results = multi_comp(Listfilename.c_str(), Listprefix.c_str(), samples, Dist_matrix, is_sim);
	return (wrap(results));	
}

//[[Rcpp::export]]
RcppExport SEXP Multi_Comp_Out(SEXP Listfilename_s, SEXP Listprefix_s, SEXP Outfilename_s, SEXP Dist_matrix_s, SEXP is_sim_s) {
	string Listfilename = as<string> (Listfilename_s);
	string Listprefix = as<string> (Listprefix_s);
	int is_sim = as<int> (is_sim_s);
	int Dist_matrix = as<int> (Dist_matrix_s);
	string outfilename = as<string> (Outfilename_s);
	vector<string> samples;
	vector<vector<float> > sim_matrix = multi_comp(Listfilename.c_str(), Listprefix.c_str(), samples, Dist_matrix, is_sim);
	Output_Matrixx(outfilename.c_str(), samples.size(), sim_matrix, samples);
	return (wrap(samples.size()));
}

*/


vector<float> Multi_Comp_Table(SEXP gene_s, SEXP pw_s, _Table_Format & abd_table, int Dist_matrix, int is_sim){
	StringVector gene(gene_s);
	StringVector pw(pw_s);
	//_Comp_Tree_Func comp_tree_func(Ref_db);
	
	_Comp_Tree_Func comp_tree_func(gene, pw);
	unsigned long file_count = (unsigned long) abd_table.Get_Sample_Size();
	//Rcpp::Rcout << "file_count = " << file_count << endl;	
	//load abd
	float **Abd = new float * [file_count];
	
	for (unsigned long i = 0; i < file_count; i ++){
		Abd[i] = new float [comp_tree_func.Get_GeneN()];
		//Rcpp::Rcout << comp_tree_func.Load_Gene_Count(&abd_table, Abd[i], i) << " KOs in file " << i + 1 << endl;
		comp_tree_func.Load_Gene_Count(&abd_table, Abd[i], i);
	}
	vector<float> results(file_count * file_count);
	
	omp_set_num_threads(Coren);        
        #pragma omp parallel for schedule(dynamic, 1)
	for(unsigned long i = 0; i < file_count; i++) {
		results[i*file_count+i] = 1.0;
		for(unsigned long j = i + 1; j < file_count; j++) {
			results[i*file_count+j] = results[j*file_count+i] = comp_tree_func.Calc_sim(Abd[i], Abd[j], Dist_matrix);			
		}
	}
	
	if(!is_sim) {
		#pragma omp parallel for schedule(dynamic, 1)
		for(unsigned long i = 0; i < file_count; i++) {
                	for(unsigned long j = 0; j <= i; j++) {
                        results[i*file_count+j] = results[j*file_count+i] = 1 - results[i*file_count+j];
        	        }
        	}
	}
	return results;
	/*
	//Rcpp::Rcout << file_count << " files loaded" << endl;
	
	//make order
	vector <int> order_m;
	vector <int> order_n;
	unsigned long iter = 0;
	for (int i = 0; i < file_count - 1; i ++)
		for (int j = i + 1; j < file_count; j ++){			
			order_m.push_back(i);
			order_n.push_back(j);
			iter ++;
		}
		
	vector <float>  sim_matrix(iter, 0);
  	
	unsigned long m, n, p; 
	//openmp	 
	omp_set_num_threads(Coren);
	
	#pragma omp parallel for private(m, n, p)
	for (unsigned long i = 0; i < iter; i ++){
		
		m = order_m[i];
		n = order_n[i];
		p = m * (unsigned long) file_count + n - (1 + m + 1) * (m + 1) / 2;
		
		sim_matrix[p] = comp_tree_func.Calc_sim(Abd[m], Abd[n], Dist_matrix);
	}

	vector<float> results = manage_matrix(file_count, sim_matrix, is_sim);
	return results;
	*/
}


//[[Rcpp::export]]	
RcppExport SEXP Multi_Comp_Table(SEXP gene_s, SEXP pw_s, SEXP features_s, SEXP samples_s, SEXP abds_s, SEXP Dist_matrix_s, SEXP is_sim_s) {
	StringVector features(features_s);
	StringVector samples(samples_s);
	NumericMatrix abds(abds_s);
	_Table_Format table(features, samples, abds);
	int Dist_matrix = as<int> (Dist_matrix_s);
	check(Dist_matrix);
	int is_sim = as<int> (is_sim_s);
	vector<float> results = Multi_Comp_Table(gene_s, pw_s, table, Dist_matrix, is_sim); 
	return (wrap(results));
}

/*
//[[Rcpp::export]]	
RcppExport SEXP Multi_Comp_Table_Out(SEXP Tablefilename_s, SEXP Outfilename_s, SEXP Reversed_table_s, SEXP Dist_matrix_s, SEXP is_sim_s) {
	string Tablefilename = as<string> (Tablefilename_s);
	int Reversed_table = as<int> (Reversed_table_s);
	_Table_Format table(Tablefilename.c_str(), Reversed_table);

	vector<string> samples = table.Get_Sample_Names();

	int Dist_matrix = as<int> (Dist_matrix_s);
	check(Dist_matrix);

	int is_sim = as<int> (is_sim_s);
	vector<vector<float>> sim_matrix = Multi_Comp_Table(table, Dist_matrix, is_sim);
 
	string outfilename = as<string> (Outfilename_s);
	Output_Matrixx(outfilename.c_str(), samples.size(), sim_matrix, samples);
	return (wrap(samples.size()));
}
*/

