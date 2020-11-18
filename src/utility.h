// Updated at Aug 28, 2019
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <RcppArmadillo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/dir.h>
#include <sys/stat.h>

#define BUFFER_SIZE 5000

#include "hash.h"
using namespace std;

string Check_Env();


int dcmp(long double x);

int Check_Path(const char * path, int type);

bool Check_Path(const char * path);

bool Check_File(const char * file);

unsigned int Get_Count(const char * infilename);

int Check_Format(const char * infilename);

string Check_OTU(string otu);

string Check_SP(string sp);

int Load_ID(const char * idfilename, vector <string> & ID, int skip);

int Load_ID(const char * idfilename, vector <string> & ID);

int Load_List(const char * listfilename, vector <string> & list);

int Load_List(const char * listfilename, vector <string> & list, string prefix); //with prefix

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids); //with id

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids, string prefix); //with id

void Make_list(const char * listname, const char * outpathname, vector <string> ids, int mode); //0: classification.txt; 1: functions.txt; 2: id

void Add_list_prefix(const char * inlistname, const char * prefix, const char * outlistname); //for prefix


#endif
