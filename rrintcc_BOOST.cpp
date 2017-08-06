/* Copyright (C) (2017) (Sen ZHANG) <szhangat@connect.ust.hk

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

To collect contigency tables,
We first pre-calculate the number of 1's of 0~2^16 using the hamming weight to count the number of 1's in the string. (see http://en.wikipedia.org/wiki/Hamming_weight)
see the function bitCount. We store them in the vector "wordbits". Then use table looking method. see function popcount.

The contigency table collection part is modified from BOOSTx64.c (Can YANG, 2010)

*/

// Format: ./rrintcc_BOOST --config configNames.txt --silent --max-cov 10000

#include "utility.h"

using namespace std;

ofstream LOG;

int main(int argc, char* argv[])
{
	/* Declare variable */
	int    numSets 		 = 1;
	string foutpath, resname, logname, configname;
	string filename, mapname, setpath, setname;

	configname = "configNames.txt"

	// Used for set-set interaction tests
    // Warning: Should not be changed!!!
    bool show_message	 = true;
	bool skip_symm 		 = false;
	bool set_test 		 = true;

	// Used for combined p value calculation
	double myth_pgates 	 = 0.05;
	double myth_trun  	 = 0.05;
	int    reps 	     = 1000;   // How many permutations will be performed if flagperm = true
	bool   flagperm 	 = false;  // whether do permutations for ptts and ptprod or not 
    int    max_cov_cnt   = 20000;   // default: 20000 (means the cov matrix is at most 20000*20000)

    // Read data from argv[]
    for (int i = 1; i < argc; i++)
    {
    	if (strcmp(argv[i], "--config") == 0 || strcmp(argv[i], "-c") == 0)
    		configname = argv[i+1];

    	if (strcmp(argv[i], "--silent") == 0 || strcmp(argv[i], "-s") == 0) 
    	{
    		show_message = false;
    	}

    	if (strcmp(argv[i], "--max-cov") == 0) 
    	{
    		if (argv[i+1] == NULL || ! atoi(argv[i+1])) 
    		{
    			printf("err: illegal input for --max-cov\n");
    			exit(1);
    		}
    		max_cov_cnt = atoi(argv[i+1]);
    	}

    	//if (strcmp(argv[i], "--set") == 0 || strcmp(argv[i], "-s") == 0) 

    }

	int *DataSize;
	int ndataset;

	// Change pheno from pointer to vector (for shuffle reason to do permutation test). Keep geno and geno_bar as pointers because p is large while n is ok

	vector<bool> pheno; // pheno: length n (bool)
	BYTE **geno; // geno: length(n*p) (BYTE)
	double **geno_bar; //geno_bar: leng(p*2) (double) [0]: ctrl; [1]: case;

	vector<int> snpchr;
	vector<string> snpname;

	vector<int> sA, sB; // Use vector<int> rather than vector<bool> to speed up
	double pmin;

	int n, p, ncase, nctrl;;  // n: number of samples; p: number of varibles
	
	RInside R(argc, argv);
	R.parseEvalQ("library(mvtnorm); library(corpcor)");
	
	clock_t st, ed;

	/////////////////////////////////////////////////////
	// Calc file names
	if (show_message) 
	{
		printf("-----------------------------------------\n");
		printf("start getting the file names...\n");
	}
	st = clock();
	GetFileNames(configname, foutpath, resname, logname, filename, mapname, setpath, setname, show_message);
	ed = clock();
	if (show_message) 
	{
		printf("cputime for getting file names: %f seconds.\n", (double)(ed - st)/CLOCKS_PER_SEC);
	}

	///////////////////////////////////////////////////////
	// Initialize .log output
	// Only record results after GetFileNames
	LOG.open(logname.c_str(), ios::out);
	LOG.clear();

	/////////////////////////////////////////////////////
	// Calc data size
	if (show_message) 
	{
		printf("-----------------------------------------\n");
		printf("start getting the data size...\n");
	}
	st = clock();
	GetDataSize(filename, &DataSize, ndataset, show_message);
	ed = clock();
	if (show_message) 
	{
		printf("cputime for getting data size: %f seconds.\n", (double)(ed - st)/CLOCKS_PER_SEC);
	}

	// load .map data (snp chromosome and snp IDs)
	if (show_message) 
	{
		printf("-----------------------------------------\n");
		printf("start reading the map file...\n");
	}
	st = clock();
	GetSnpInfo(mapname, snpchr, snpname, show_message);
	ed = clock();
	if (show_message) 
	{
		printf("cputime for reading the map file: %f seconds.\n", (double)(ed - st)/ CLOCKS_PER_SEC);
	}

	// load BOOST.txt data to pheno (n), geno (n*p) and geno_bar (p*2)
	if (show_message) 
	{
		printf("-----------------------------------------\n");
		printf("start reading the BOOST file...\n");
	}
	st = clock();
	GetData(filename, DataSize, n, p, ncase, nctrl, ndataset, pheno, &geno, &geno_bar, show_message);
	ed = clock();
	if (show_message) 
	{
		printf("cputime for reading the BOOST file: %f seconds.\n", (double)(ed - st)/ CLOCKS_PER_SEC);
		printf("-----------------------------------------\n");
		printf("The number of snps: %d\n", p);
		printf("The number of samples: %d (ncase = %d; nctrl = %d)\n", n, ncase, nctrl);
	}

	//////////////////////////////////////////////////////////////////////////////
	if (show_message) 
	{
		printf("-----------------------------------------\n");
		printf("start calculating the region interactions...\n");
	}
//	time(&st);
	for(int i = 0; i < numSets; i++)
	{
		st = clock();
		if (show_message && i > 0 && i%1000 == 0)
		{
			printf("%d sets have been analyzed\n", i);
		}
		
		
		//string setname;
		//setname = setpath + "locipair" + to_string(i+1) + ".set";

		string fout = foutpath;
		fout.append("snp_pair_results");
		fout.append(to_string(i+1));
		fout.append(".txt");


		// load .set data: Write sA, sB, and skip_symm inside
		sA.clear();
		sB.clear();
		skip_symm = false; // Need to reset skip_symm, sA, sB!!!
		GetSetInfo(setname, snpname, sA, sB, skip_symm, set_test, p, show_message);

		pmin = CalcRegionInter(R, fout, pheno, geno, geno_bar, snpchr, snpname, skip_symm, p, n, ncase, nctrl, sA, sB, myth_pgates, myth_trun, reps, flagperm, max_cov_cnt, show_message);

		// ofstream::app: Appending to the last line of current file
		// ios::out: Rewrite current file (See utility.cpp) 
		ofstream EPI;
		EPI.open(resname.c_str(), ofstream::app);
		EPI.precision(4);
		EPI << setw(8)  << "Pair " << to_string(i+1) << " | " 
			<< setw(8)  << "pmin: " << " "
			<< setw(15) << pmin << "\n";
		EPI.flush();
		EPI.close();

		ed = clock();

		if (show_message) 
		{	printf("cputime for calculating the region interactions: %f seconds.\n", (double)(ed - st)/ CLOCKS_PER_SEC);
		}
	}	


	// time(&ed);
	// printf("cputime for calculating the region interactions: %d seconds.\n", (int)ed - st);

	//free Datasize, geno and geno_bar
	free(DataSize);

	for(int i = 0; i < n; i++)
		free(geno[i]);
	free(geno);

	for(int i = 0; i < p; i++)
		free(geno_bar[i]);
	free(geno_bar);

	LOG.close();
	return 1;

} // end main()






