#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> 
#include <numeric>
#include <vector>
#include <algorithm>
#include <string>
#include <vector>
#include <random>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <RInside.h>

#include "stats.h"

using namespace std;

extern ofstream LOG;

///////////////////////////////////////////
// Type definitions 
///////////////////////////////////////////

// Define matrix_t as 2-D vector 
typedef vector<vector<double> > matrix_t;
// Type definitions from BOOST
typedef long long   int64;
typedef unsigned char BYTE;
typedef unsigned long long uint64;
#define FMT_INT64   "%lld"
#define FMT_UINT64   "%llu"
#define FMT_HEX64   "%llx"
#define pp_maxsnp 	6

///////////////////////////////////////////
// Function definitions 
///////////////////////////////////////////
void printLOG(string s);

string int2str(int n);

string dbl2str(double n, int prc);

string char2str(char *f);

double Abs(double a);

void sizeMatrix(vector<vector<double> > &m, int r, int c);

void GetSnpInfo(char *filename, vector<int> &snpchr, vector<string> &snpname);

void GetSetInfo(char *setname, vector<string> &snpname, vector<int> &sA, vector<int> &sB, bool &skip_symm, bool set_test, int p);

void GetDataSize(char *filename, int **DataSize, int &ndataset_out);

void GetData(char *filename, int *DataSize, int &n, int &p, int &ncase, int &nctrl, int ndataset, vector<bool> &pheno, BYTE ***geno, double ***geno_bar);

vector<double> CalcRegionInter(RInside &R, string fout, vector<bool> &pheno, BYTE **geno, double **geno_bar, vector<int> &snpchr, vector<string> &snpname, bool skip_symm, int p, int n, int ncase, int nctrl, vector<int> &sA, vector<int> &sB, double myth_pgates, double myth_trun, int reps, bool flagperm, int max_cov_cnt);

double CalcPgates(vector<double> &plist, Rcpp::NumericMatrix &corr_matrix, double myth_pgates);

void LDContrastTest(vector<bool> &pheno, vector<double> &zlist, vector<double> &plist, vector<int> &cov_index, BYTE **geno, double **geno_bar, vector<int> &sA, vector<int> &sB, bool skip_symm, int p, int n, int ncase, int nctrl);

vector<double> CalcPerm(int reps, vector<bool> &pheno, vector<double> &zlist, vector<double> &plist, BYTE **geno, double **geno_bar, vector<int> &sA, vector<int> &sB, bool skip_symm, int p, int n, int ncase, int nctrl, double myth_trun);

int CountSig(vector<double> &vec_p, double alpha);

#endif
