// Updated at Aug 12, 2019
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
//version 3.1 or above with _Table_Format
// Updated by Xiaoquan Su
// Add_feature and Add abd
#pragma once
#ifndef TABLE_FORMAT_H
#define TABLE_FORMAT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

class _Table_Format{
    
    public:
    
    friend class _Table_Format_Seq;
    
    _Table_Format(){}    
    _Table_Format(StringVector & features, StringVector & samples, NumericMatrix & abds); //default is true
    
    void Filter_Max(float max); // Filter features by the Maximum Value
    void Filter_Min(float min); // Filter features by the Minimum Value
    void Filter_Ave(float ave); // Filter features by the average Value
    void Filter_Zero(float zero); // Filter features by the proportion of non-zero values
    void Filter_Empty(); //Filter empty features
        

    StringVector Get_Feature_Names();
    StringVector Get_Sample_Names();
    int Get_Sample_Size();//Get the size of Samples
    int Get_Feature_Size();//Get the size of Features
    
    float Get_Abd_By_Order(unsigned int s, unsigned int i);
    float Get_Abd_By_Feature(unsigned int s, string a_feature);
    NumericVector Get_Abd(unsigned int s);
    
    float Calc_Dist_Cos(int sam_m, int sam_n); //Calculate the Cosin distance between two samples
    float Calc_Dist_E(int sam_m, int sam_n); //Calculate the Euclidean distance between two samples
    float Calc_Dist_JSD(int sam_m, int sam_n); //Calculate the Jeason-Shannon distance between two samples
    float Calc_Dist_Bray_Curtis(int sam_m, int sam_n); //Calculate the Bray-Curtis distance between two samples
    
    float Calc_Corr_S(int sam_m, int sam_n);//Calculate Spearman correlation coefficient
    float Calc_Corr_P(int sam_m, int sam_n);//Calculate Pearson correlation coefficient
    
    float Calc_Corr_S(vector <float> & sam_m, vector <float> & sam_n);
    float Calc_Corr_P(vector <float> & sam_m, vector <float> & sam_n);
    
    void Calc_Dist_Matrixx(const char * outfilename, int metrics, int coren, bool is_sim); //0: cost 1: eu dist 2: jsd
    void Calc_Corr_Matrixx(const char * outfilename, int metrics, int coren); //0: s-corr 1: p-corr
    
	
    protected:
    StringVector Samples;
    StringVector Features;
    NumericMatrix Abd;
    
    map <string, int> Feature_hash;
    
    bool * Max_filtered;
    bool * Min_filtered;
    bool * Ave_filtered;
    bool * Zero_filtered;
    bool * Empty_filtered;
    
    void Init_Filter(); //init features
    bool Check_Filter(string a_feature);//check features to determine whether to filter
    void Build_Feature_Hash();// Build the hash numbers of features.
    void BubbleSort(float *array1, int *rank1, int len);//Bubble sort code, was used by the function named 'Calc_Corr_S'

};

#endif /* TABLE_FORMAT_H */
