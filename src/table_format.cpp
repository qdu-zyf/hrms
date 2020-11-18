// Updated at Aug 12, 2019
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
//version 3.1 or above with _Table_Format
// Updated by Xiaoquan Su
// Add_feature and Add abd

#include "table_format.h"

_Table_Format::_Table_Format(StringVector & features, StringVector & samples, NumericMatrix & abds){ 
	Features = features;
	Samples = samples;
	Abd = abds;
	Build_Feature_Hash();
	Init_Filter();	
}

StringVector _Table_Format::Get_Feature_Names(){
    return Features;
}
    
StringVector _Table_Format::Get_Sample_Names() {
    return Samples;
}

int _Table_Format::Get_Feature_Size(){
    return Features.size();
}

int _Table_Format::Get_Sample_Size(){
    return Samples.size();
}
 
void _Table_Format::Filter_Max(float max){

      float max_temp;
      unsigned long n = Features.size();

      for (int i = 0; i < (int) Features.size(); i ++){
          max_temp = 0.0;
          for (int j = 0; j < (int) Samples.size(); j++)
               max_temp = max_temp > Abd[j*n+i] ? max_temp : Abd[j*n+i];

          if (max_temp < max)
             Max_filtered[i] = true;
           } 
      }
//
       
void _Table_Format::Filter_Min(float min){

      float min_temp;
      unsigned long n = Features.size();

      for (int i = 0; i < (int)Features.size(); i ++){
          min_temp = 1.0;
          for (int j = 0; j < (int)Samples.size(); j++)
               min_temp = min_temp < Abd[j*n+i] ? min_temp : Abd[j*n+i];
              
          if (min_temp < min) 
             Min_filtered[i] = true;
             
          } 

      }
// 

void _Table_Format::Filter_Ave(float ave){

      float samplesum, sampleave;      
      unsigned long n = Features.size();

      for (int i = 0; i < (int) Features.size(); i ++){                                                                                                                               
          samplesum = 0;                                                                                                                                                      
          for (int j = 0; j < (int) Samples.size(); j++)                                                                                                                          
              samplesum = samplesum + Abd[j*n+i];                                                                                                               
                                
              sampleave = samplesum/Samples.size();

          if (sampleave < ave)
             Ave_filtered[i] = true;
                                                
            }                                                                                                                                                                
                                                                                                                                                                             
      }
//
 
void _Table_Format::Filter_Zero(float zero){
     
     unsigned long n = Features.size();

      for (int i = 0; i < (int) Features.size(); i ++){                                                                                                                               
          int samplesum = 0;                                                                                                                                         
          for (int j = 0; j < (int) Samples.size(); j++)    
              if (Abd[j*n+i] > 0)
                 samplesum ++;                                                                                                                             
           
          if ((float) samplesum / (float) Samples.size() < zero)
             Zero_filtered[i] = true;
            }                                                                                                                                                                            
      }

void _Table_Format::Filter_Empty(){
	
	unsigned long n = Features.size();

    for (int i = 0; i < (int) Features.size(); i ++){                                                                                                                               
          bool is_empty = true;                                                                                                                                         
          for (int j = 0; j < (int) Samples.size(); j++)    
              if (Abd[j*n+i] > 0)
                 is_empty = false;                                                                                                                            
           
             Empty_filtered[i] = is_empty;
            }                                                                                                                                                                              
      }

float _Table_Format::Calc_Dist_E(int sam_m, int sam_n){

     float abd_m_norm[Features.size()], abd_n_norm[Features.size()];
     float sum_m = 0;
     float sum_n = 0;
     float     f = 0;
     unsigned long n = Features.size(); 

     for (int i = 0; i < Features.size(); i++){
         abd_m_norm[i] = Abd[sam_m * n + i];
         abd_n_norm[i] = Abd[sam_n * n + i];
 
         sum_m += abd_m_norm[i];
         sum_n += abd_n_norm[i];
         }
     
     if (sum_m <= 0) return 1;
     if (sum_n <= 0) return 1;
     
     for (int i = 0; i < Features.size(); i++){
         if (sum_m > 0) abd_m_norm[i] /= sum_m;
         if (sum_n > 0) abd_n_norm[i] /= sum_n;
         }

     
     for (int i = 0; i < Features.size(); i++)
         f += pow(abd_m_norm[i] - abd_n_norm[i], 2);

     return sqrt(f);
     }

float _Table_Format::Calc_Dist_Cos(int sam_m, int sam_n){

     float f_m_sum = 0;
     float f_n_sum = 0;
     float f = 0;
     unsigned long n = Features.size();

     for (int i = 0; i < Features.size(); i++){
         
         float fm = Abd[sam_m * n + i];
         float fn = Abd[sam_n * n + i];

         f += fm * fn;

         f_m_sum += fm * fm;
         f_n_sum += fn * fn;
         }
     
     float ff = sqrt(f_m_sum) * sqrt(f_n_sum);
     if (ff == 0) return 0;
     
     f = f / ff;
     
     return (1 - f);
     }

float _Table_Format::Calc_Dist_JSD(int sam_m, int sam_n){

      float abd_m_norm[Features.size()];
      float abd_n_norm[Features.size()];
      
      float sum_m = 0;
      float sum_n = 0;
      
      unsigned long n = Features.size();

      //Norm
      for (int i = 0; i < Features.size(); i++){
         if (Abd[sam_m * n + i] > 0) abd_m_norm[i] = Abd[sam_m * n + i];
         else abd_m_norm[i] = 0;
         
         if (Abd[sam_n * n + i] > 0) abd_n_norm[i] = Abd[sam_n * n + i];
         else abd_n_norm[i] = 0;
 
         sum_m += abd_m_norm[i];
         sum_n += abd_n_norm[i];
         }
     
     if (sum_m <= 0) return 1;
     if (sum_n <= 0) return 1;
     
     for (int i = 0; i < Features.size(); i++){
         if (sum_m > 0) abd_m_norm[i] /= sum_m;
         if (sum_n > 0) abd_n_norm[i] /= sum_n;
         }
      
      //calc
      float dkl_m = 0;
      float dkl_n = 0;
      
      for (int i = 0; i < Features.size(); i ++){
          
          if ((abd_m_norm[i] == 0) && (abd_n_norm[i] == 0)) continue;
          
          float abd_q = (abd_m_norm[i] +  abd_n_norm[i]) / 2;
          
          if (abd_m_norm[i] > 0)
             dkl_m += abd_m_norm[i] * log(abd_m_norm[i] / abd_q);
          
          if (abd_n_norm[i] > 0)
             dkl_n += abd_n_norm[i] * log(abd_n_norm[i] / abd_q);

          }
          
      return (dkl_m + dkl_n)/2.0;
     }

float _Table_Format::Calc_Dist_Bray_Curtis(int sam_m, int sam_n){

      float abd_m_norm[Features.size()];
      float abd_n_norm[Features.size()];
      
      float sum_m = 0;
      float sum_n = 0;
      
      unsigned long n = Features.size();

      //Norm
      for (int i = 0; i < Features.size(); i++){
         if (Abd[sam_m * n + i] > 0) abd_m_norm[i] = Abd[sam_m * n + i];
         else abd_m_norm[i] = 0;
         
         if (Abd[sam_n * n + i] > 0) abd_n_norm[i] = Abd[sam_n * n + i];
         else abd_n_norm[i] = 0;
 
         sum_m += abd_m_norm[i];
         sum_n += abd_n_norm[i];
         }
     
     if (sum_m <= 0) return 1;
     if (sum_n <= 0) return 1;
     
     for (int i = 0; i < Features.size(); i++){
         if (sum_m > 0) abd_m_norm[i] /= sum_m;
         if (sum_n > 0) abd_n_norm[i] /= sum_n;
         }
      
      //calc
      float sum = 0;
      float diff = 0;
            
      for (int i = 0; i < Features.size(); i ++){
          sum += (abd_m_norm[i] + abd_n_norm[i]);
          float a_diff = abd_m_norm[i] - abd_n_norm[i];
          if (a_diff < 0) a_diff = a_diff * (-1.0);
          diff += a_diff;
          }
      if (sum <= 0) return 1;
      return diff / sum;
     }

float _Table_Format::Calc_Corr_S(vector <float> & sam_m, vector <float> & sam_n){

      float sumValue;
      float partValue, rho;
      int sample_num1, sample_num2, rank_tempm, rank_tempn;
      int mline = sam_m.size();
      unsigned long nline = sam_n.size();
      
      if (mline!=nline){
          Rcpp::Rcout << "Error: Different size of two vectors, check your input files! "<<endl;
          return 0;
                        }
      float abdm[mline], abdn[mline];
      int rankn[mline], rankm[mline];
      int rankn1[mline], rankm1[mline];
      float rankn2[mline], rankm2[mline];

      for (int i = 0; i < mline; i++){
          abdm[i] = sam_m[i];
          abdn[i] = sam_n[i];
          rankm[i] = i;
          rankn[i] = i;
          }
      
      BubbleSort(abdm, rankm, mline);
      BubbleSort(abdn, rankn, mline);
        
      for (int i = 0; i < mline; i++){
          for (int j = 0; j < mline; j++){
              if (sam_m[i] == abdm[j] && rankm[j] == i) rankm1[i]=j;
              if (sam_n[i] == abdn[j] && rankn[j] == i) rankn1[i]=j;
              }
          }
    
      for (int i = 0; i < mline; i++){

          rank_tempm = 0;
          rank_tempn = 0;
          sample_num1 = 0;
          sample_num2 = 0; 

          for (int j = 0; j < mline; j++){
              if (sam_m[i] == sam_m[j]) {
                 rank_tempm += rankm1[j];     
                 sample_num1 ++;
                 }
              if (sam_n[i] == sam_n[j]) {
                 rank_tempn += rankn1[j];
                 sample_num2 ++;
                 } 
              }
          rankm2[i] = (float) rank_tempm / sample_num1;
          rankn2[i] = (float) rank_tempn / sample_num2;
          
          }
      
      sumValue = 0.0; 

      for (int i = 0; i < mline; i++)
          sumValue += (rankm2[i] - rankn2[i])*(rankm2[i] - rankn2[i]); 
      
      partValue = mline*(mline*mline - 1);

      rho = 1 - ((6*sumValue)/partValue);
      
      return rho;

      }

float _Table_Format:: Calc_Corr_P(vector <float> & sam_m, vector <float> & sam_n){

      float ave_m = 0;
      float ave_n = 0;
      int mline = sam_m.size();
      unsigned long nline = sam_n.size();
      
      if (mline!=nline){
          Rcpp::Rcout << "Error: Different size of two vectors, check your input files! "<<endl;
          return 0;
                        }
      
      for (int i = 0; i < mline; i++){
          
          ave_m += sam_m[i];
          ave_n += sam_n[i];

          }
      
      ave_m = ave_m / (float) mline;
      ave_n = ave_n / (float) nline;


      float deno_m = 0;
      float deno_n = 0;
 
      for (int i = 0; i < mline; i++){
           
          deno_m += pow(sam_m[i] - ave_m, 2);
          deno_n += pow(sam_n[i] - ave_n, 2);

          }

      float deno = sqrt(deno_m) * sqrt(deno_n);
      
      if (deno == 0) return 0;

      float nume = 0;

      for (int i = 0; i < mline; i++)
          nume += (sam_m[i] - ave_m)*(sam_n[i] - ave_n);


      return nume/deno;
      }

float _Table_Format::Calc_Corr_P(int sam_m, int sam_n){

      return Calc_Corr_P(Abd[sam_m], Abd[sam_n]);

      }

float _Table_Format::Calc_Corr_S(int sam_m, int sam_n){

           
      return Calc_Corr_S(Abd[sam_m], Abd[sam_n]);

      }


void _Table_Format::Calc_Dist_Matrixx(const char * outfilename, int metrics, int coren, bool is_sim){
     
      //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < (int) Samples.size() - 1; i ++)
        for (int j = i + 1; j < (int) Samples.size(); j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
      
      ofstream outfile(outfilename, ofstream::out);                                 
      if (!outfile){                                                                
         Rcpp::Rcout << "Error: Cannot open output file: " << outfilename << endl;               
         return;                                                                   
         } 
    
      //calc dist 
      vector <float>  dist_matrix;
      for (long i = 0; i < iter; i ++)
        dist_matrix.push_back(0);
      
      long m, n, p;
      omp_set_num_threads(coren);
      
      #pragma omp parallel for private(m, n, p)
      for (long i = 0; i < iter; i ++){
        
         m = order_m[i];
         n = order_n[i];
         p = m * Samples.size() + n - (1 + m + 1) * (m + 1) / 2;    
        
         switch (metrics){
                case 0: dist_matrix[p] = Calc_Dist_Cos(m, n); break;
                case 1: dist_matrix[p] = Calc_Dist_E(m, n); break;
                case 2: dist_matrix[p] = Calc_Dist_JSD(m, n); break;
                case 3: dist_matrix[p] = Calc_Dist_Bray_Curtis(m, n); break;
                default: dist_matrix[p] = Calc_Dist_Cos(m, n); break;
                }
        }
      
       for (int i = 0; i < (int)Samples.size(); i++)
          outfile <<"\t" <<Samples[i];
       outfile <<endl;
      
       for (int i = 0; i < (int)Samples.size(); i++){
          outfile <<Samples[i];
          for (int j = 0; j < (int)Samples.size(); j++){
              long m = (i <= j) ? i : j;
              long n = (i > j) ? i : j;
              long p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2; 
              if (is_sim) {
                          if (i == j) outfile << "\t" << 1;
                          else outfile << "\t" << 1 - dist_matrix[p];
                          }
              else {
                   if (i == j) outfile << "\t" << 0;
                   else outfile << "\t" << dist_matrix[p];
                   }
              }
          outfile <<endl;
          }

      outfile.close();
      outfile.clear();
      }

void _Table_Format::Calc_Corr_Matrixx(const char * outfilename, int metrics, int coren){
     
      //make order
    vector <int> order_m;
    vector <int> order_n;
    long iter = 0;
    for (int i = 0; i < (int)Samples.size() - 1; i ++)
        for (int j = i + 1; j < (int)Samples.size(); j ++){            
            order_m.push_back(i);
            order_n.push_back(j);
            iter ++;
            }
      
      ofstream outfile(outfilename, ofstream::out);                                 
      if (!outfile){                                                                
         Rcpp::Rcout << "Error: Cannot open output file: " << outfilename << endl;               
         return;                                                                   
         } 
    
      //calc dist 
      vector <float>  corr_matrix;
      for (long i = 0; i < iter; i ++)
        corr_matrix.push_back(0);
      
	long m, n, p;
      omp_set_num_threads(coren);
      
      #pragma omp parallel for private(m, n, p)
      for (long i = 0; i < iter; i ++){
        
         m = order_m[i];
         n = order_n[i];
         p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2;    
        
         if (metrics == 0)
                    corr_matrix[p] = Calc_Corr_S(m, n);
         else corr_matrix[p] = Calc_Corr_P(m, n);

        }
      
       for (int i = 0; i < (int)Samples.size(); i++)
          outfile <<"\t" <<Samples[i];
       outfile <<endl;
      
       for (int i = 0; i < (int)Samples.size(); i++){
          outfile <<Samples[i];
          for (int j = 0; j < (int)Samples.size(); j++){
              long m = (i <= j) ? i : j;
              long n = (i > j) ? i : j;
              long p = m * (long) Samples.size() + n - (1 + m + 1) * (m + 1) / 2; 
              if (i == j) outfile << "\t" << 1;
                          else outfile << "\t" <<corr_matrix[p];
                          }             
          outfile <<endl;
          }

      outfile.close();
      outfile.clear();
      }

void _Table_Format::BubbleSort(float *array1, int *rank1, int len){

     float temp;
     int   temp_rank;

     for (int i = 0; i < len-1; i++){
         for (int j = 0; j < len-1-i; j++){
             if (array1[j] > array1[j+1]){
                temp        = array1[j];
                array1[j]   = array1[j+1];
                array1[j+1] = temp;
                temp_rank  = rank1[j];
                rank1[j]   = rank1[j+1];
                rank1[j+1] = temp_rank;
                }
             }
         }
     }
     
float _Table_Format::Get_Abd_By_Feature(unsigned int s, string a_feature){
            
      if (s >= Samples.size())
         return 0;
      if (Feature_hash.count(a_feature) == 0)
         return 0;
      if (Check_Filter(a_feature))
         return 0;
      
      unsigned long n = Features.size();
      return Abd[s * n + Feature_hash[a_feature]];
      }
      
float _Table_Format::Get_Abd_By_Order(unsigned int s, unsigned int i){
            
      if (s >= Samples.size() || (i >= Features.size()))
         return 0;

      unsigned long n = Features.size();
      return Abd[s * n + i];
      }      

NumericVector _Table_Format::Get_Abd(unsigned int s){
       
       if (s < Samples.size())
          return Abd.row(s);
       else {
	   NumericVector tmp;
	   return tmp;
           }
       }

void _Table_Format::Init_Filter(){
     Max_filtered = new bool[Features.size()];
     Min_filtered = new bool[Features.size()];
     Ave_filtered = new bool[Features.size()];
     Zero_filtered = new bool[Features.size()];
     Empty_filtered = new bool[Features.size()];
     
     for (int i = 0; i < (int)Features.size(); i ++){
         Max_filtered[i] = false;
         Min_filtered[i] = false;
         Ave_filtered[i] = false;
         Zero_filtered[i] = false;
         Empty_filtered[i] = false;
         }
     }

bool _Table_Format::Check_Filter(string a_feature){
     
     if (Feature_hash.count(a_feature) == 0) return true;
     int feature_index = Feature_hash[a_feature];    
     return (Max_filtered[feature_index] || Min_filtered[feature_index] || Ave_filtered[feature_index] || Zero_filtered[feature_index] || Empty_filtered[feature_index]);
     }
    
void _Table_Format::Build_Feature_Hash(){
         
     	for (int i = 0; i < (int)Features.size(); i ++) {
         	Feature_hash[as<string> (Features[i])] = i;
	}
}
     
// #############################################################################
//#endif /* TABLE_FORMAT_H */
