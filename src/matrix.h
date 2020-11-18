#pragma GCC optimize(3,"Ofast","inline")
#ifndef _MATRIX_H
#define _MATRIX_H

#include "utility.h"
#include <omp.h>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include <RcppArmadillo.h>
#include <ctime>

using namespace std;
using namespace Eigen;
using namespace Rcpp;

float dbEps = 1e-9;

class Matrixx {

public:
	Matrixx():n(0) {};
	Matrixx(const NumericMatrix &m);
	MatrixXf get_matrix();
	void get_Deviation_Matrix();
	NumericMatrix Get_PC_Matrix(int k);
	NumericVector get_percentage(int k);

private:
	int n;
	MatrixXf matrix;
	MatrixXf pdblVects; //特征向量
	VectorXf pdbEigenValues; //特征值矩阵（特征值在对角线上）
};

Matrixx::Matrixx(const NumericMatrix &m) {
	this->n = (int) sqrt(m.size());
	matrix.resize(n, n);
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i = 0; i < n; i++) {
		for(int j = 0; j <= i; j++) {
			matrix(i, j) = matrix(j, i) = m[i*n+j];
		}
	}
}

MatrixXf Matrixx::get_matrix() {
	return matrix;
}

void Matrixx::get_Deviation_Matrix() {
        /*
        MatrixXf one = MatrixXf::Constant(n, n, 1);
        MatrixXf mat(n, n);
        mat.setIdentity();
        mat = mat - one / n;
        float para = 1.0;

        #pragma omp parallel for schedule(dynamic, 1)
        for(int i = 0; i < n; i++) {
                for(int j = 0; j <= i; j++) {
                        matrix(i, j) = matrix(j, i) = matrix(i, j) * matrix(i, j) * para;
                }
        }

        matrix =  -mat * matrix * mat;
        */
        #pragma omp parallel for schedule(dynamic, 1)   
        for(int i = 0; i < n; i++) {
                for(int j = 0; j <= i; j++) {
                        matrix(i, j) = matrix(j, i) = matrix(i, j) * matrix(i, j);
                }
        }
        
        float * row_sums = (float *) malloc(sizeof(float) * n);
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
                row_sums[i] = 0;
        }
        //#pragma omp parallel for
        for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                        row_sums[i] += matrix(i, j);
                }
        }

        float sums = 0;
        #pragma omp parallel for reduction(+:sums)
        for(int i = 0; i < n; i++) {
                sums += row_sums[i];
        }
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                        matrix(i, j) = -0.5 * matrix(i, j) + 1.0 / (2 * n) * (row_sums[i] + row_sums[j]) - 1.0 / (2 * n * n) * sums;
                }
        }
}

NumericMatrix Matrixx::Get_PC_Matrix(int k) {
	if(k >= n) {
		Rcpp::Rcout << "'k' must be in {1, 2, ..  n - 1}" << endl;
		NumericMatrix tmp;
		return tmp;
	}
	
	SelfAdjointEigenSolver<MatrixXf> es(matrix);
        pdbEigenValues = es.eigenvalues();
        pdblVects = es.eigenvectors();
	
	NumericMatrix pc(n, k);
	
	#pragma omp parallel for
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < k; j++) {
			pc(i, j) = pdblVects(i, n-j-1) * sqrt(pdbEigenValues(n-j-1));
		}
	}

	return pc;
}

NumericVector Matrixx::get_percentage(int k) {
		NumericVector percentage(k);
		float sum_s = 0;
		#pragma omp parallel for reduction(+:sum_s) 
		for(int i = 0; i < n; i++) {
			sum_s += matrix(i, i);
		}
		#pragma omp parallel for 
		for(int i = 0; i < k; i++) {
			percentage(i) = pdbEigenValues(n-i-1)/sum_s;
		}
		return percentage;
}

#endif
