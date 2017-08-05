#include "utility.h"

extern ofstream LOG;

void printLOG(string s)
{
	LOG << s;
	LOG.flush();
}

string int2str(int n)
{
	ostringstream s2(stringstream::out);
	s2 << n;
	return s2.str();
}

string char2str(char *f)
{
	ostringstream s2;
	s2 << f;
	return s2.str();
}

void GetSnpInfo(string filename, vector<int> &snpchr, vector<string> &snpname)
{
	string line;
	ifstream fp (filename);
	if(fp.is_open())
	{
		while(getline(fp, line))
		{
			istringstream iss(line);
			int chrtmp;
			string nametmp;
			iss >> chrtmp >> nametmp;
			snpchr.push_back(chrtmp);
			snpname.push_back(nametmp);
		}
	}
	else
	{
		fprintf(stderr,"can't open input file %s\n",filename.c_str());
		printLOG("can't open input file "+ filename +"\n");
		exit(1);
	}

}

// For vector as function input:
// (1) Need to use reference &, to avoid copying each element of the vector
// (2) If possible, use vector<int>::const_iteractor

// The indices are found based on the order in snpname. 
// snpname is read from .map file
// This means snp number and order in boost or ped should be consistent with .map file
// (.map file and BOOST file must have the same snps!!!)


void GetSetInfo(string setname, vector<string> &snpname, vector<int> &sA, vector<int> &sB, bool &skip_symm,  bool set_test, int p)
{

	vector<vector <int> > snpset;
	ifstream fset(setname);
	vector<int> tmpset; 
	string tmpsetline;
	bool flagset = false; 

	sA.clear();
	sB.clear();
	sA.resize(p, 0);
	sB.resize(p, 0);

	if (fset.is_open())
	{
		while (getline(fset, tmpsetline))
		{
            // Deal with diff end-of-line terminator in linux (\n) and windows(\r\n) files 
            if (tmpsetline.back() == '\r')
            	tmpsetline.erase(tmpsetline.length()-1);

			if (tmpsetline == "END")
			{	
				if (tmpset.size() > 0)
				{
					snpset.push_back(tmpset);
					tmpset.clear();
				}

				flagset = false;
			}

			if (flagset)
				tmpset.push_back(find(snpname.begin(), snpname.end(), tmpsetline) - snpname.begin());

			if (tmpsetline.substr(0, 3) == "SET")
				flagset = true;

		}

	}
	else
	{
		printf("cannot open the file: %s\n", setname.c_str());
		printLOG("cannot open the file: "+ setname +"\n");
		exit(1);
	}

	if(set_test)
	{
		if(snpset.size() > 2)
		{
			printf("Can only specify one or two SETs when testing for epistasis\n");
			printLOG("Can only specify one or two SETs when testing for epistasis\n");
			exit(1);
		}

		if(snpset.size() == 0)
		{
			printf("There are no valid sets specified\n");
			printLOG("There are no valid sets specified\n");
			exit(1);
		}

		for(int e = 0; e < snpset[0].size(); e++)
			sA[snpset[0][e]] = 1;

		// Has a second set been specified?

		if (snpset.size() == 2)
		{
			for (int e = 0; e < snpset[1].size(); e++)
				sB[snpset[1][e]] = 1;
		} 
		else if(snpset.size() == 1) //Otherwise, SET x SET (itself)
		{
		  	skip_symm = true;
		 	for (int e=0;e<snpset[0].size();e++)
		 	 	sB[snpset[0][e]] = 1;	
		}
	}
	else
	{
	  	skip_symm = true;
	  	for (int e=0;e<p;e++)
		{
	  		sA[e] = 1;
	  		sB[e] = 1;	  
		}
	}
}
void GetDataSize(string filename, int **DataSize, int &ndataset_out)
{
	FILE * fp, *fp_i;
	int c, ndataset;
	int n, p, i, flag, ii;
	char filename_i[100];


	fp = fopen(filename.c_str(), "r");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open input file %s\n", filename.c_str());
		printLOG("can't open input file "+filename+"\n");
		exit(1);
	}

	ndataset = 0;
	while (!feof(fp)) {
		ndataset++;
		fscanf(fp, "%s\n", &filename_i);
	}

	*DataSize = (int *)calloc(ndataset * 2, sizeof(int));

	ii = 0;
	rewind(fp);
	while (!feof(fp))
	{
		ii++;
		fscanf(fp, "%s\n", &filename_i);

		fp_i = fopen(filename_i, "r");
		if (fp_i == NULL)
		{
			fprintf(stderr, "can't open input file %s\n", filename_i);
			printLOG("can't open input file "+ char2str(filename_i) +"\n");
			exit(1);
		}
		//initialization
		if (ii == 1)
		{
			n = 0;//samples number

				  // find the number of samples: n
			while (1)
			{
				int c = fgetc(fp_i);//read a character from the data file
				switch (c)
				{
				case '\n'://the end of line
					n++;
					break;
					// fall through,
					// count the '-1' element
				case EOF://file end
					goto out;
				default:
					;
				}
			}

		}
	out:
		rewind(fp_i);//Repositions the file pointer to the beginning of a file

					 // find number of variables: p
		p = 0;
		i = 0;
		flag = 1;
		while (1)
		{
			c = getc(fp_i);
			if (c == '\n') goto out2;//end of line
			if (isspace(c))
			{
				flag = 1;
			}
			/*do {
			c = getc(fp);
			if(c=='\n') goto out2;//end of line
			} while(isspace(c));//space
			*/
			if (!isspace(c) && (flag == 1))
			{
				p++;//indicate the dimension of the vector
				flag = 0;
			}

		}
	out2:
		fclose(fp_i);

		//	DataSize[0] = n;
		(*DataSize)[ndataset * 0 + ii - 1] = n;
		(*DataSize)[ndataset * 1 + ii - 1] = p - 1;

	}
	ndataset_out = ndataset;
}

void GetData(string filename, int *DataSize, int &n, int &p, int &ncase, int &nctrl, int ndataset, vector<bool> &pheno, BYTE ***geno, double ***geno_bar)
{
	FILE * fp, *fp_i;
	char filename_i[100];
	int i, j, ii, k, tmp;
	

	n = DataSize[0];
	p = 0;
	printf("n = %d\n", n);
	printLOG("n = "+int2str(n)+"\n");

	for (i = 0; i<ndataset; i++)
	{
		p += DataSize[ndataset * 1 + i];
	}

	////////////////////////////////////////////////////////////////
	// get ncase and nctrl
	
	i = 0;
	j = 0;

	ncase = 0;
	nctrl = 0;

	fp = fopen(filename.c_str(), "r");
	if (fp == NULL)
	{
		fprintf(stderr, "can't open input file %s\n", filename.c_str());
		printLOG("Can't open input file "+filename+"\n");
		exit(1);
	}
	// only use the first file to get ncase and nctrl
	fscanf(fp, "%s\n", &filename_i);
	printf("%s\n", filename_i);
	printLOG(char2str(filename_i) + "\n");
	fp_i = fopen(filename_i, "r");

	while (!feof(fp_i)) {
		/* loop through and store the numbers into the array */
		if (j == 0)
		{
			//j = 0 means read ind class label y
			fscanf(fp_i, "%d", &tmp);

			if (tmp) // tmp=1 means case
				ncase++;
			else
				nctrl++;

			j++;
		}
		else
		{
			fscanf(fp_i, "%d", &tmp);
			j++; //column index
			if (j == (DataSize[ndataset] + 1)) // DataSize[ndataset] is the nsnp in the first dataset
			{
				j = 0;
				i++; // row index
			}

		}

		if (i >= n)
		{
			break;
		}
	}

	/////////////////////////////////////////////////////////////////
	//load data to pheno (n) and geno (n*p)
	//Allocate memories for geno
	// pheno now is rewritten as vector<bool>
	pheno.resize(n, false);

	*geno = NULL;
	*geno = (BYTE **)calloc(n, sizeof(BYTE*));
	if(*geno == NULL)
	{
		printf("Memory allocation error: no enough memory\n");
		printLOG("Memory allocation error: no enough memory\n");
		exit(1);
	}
	for(int i = 0; i < n; i++)
	{
		(*geno)[i]=(BYTE *)calloc(p, sizeof(BYTE));
		if ((*geno)[i] == NULL)
		{
			printf("Memory allocation error: no enough memory\n");
			printLOG("Memory allocation error: no enough memory\n");
			exit(1);
		}

	}

	*geno_bar = NULL;
	*geno_bar = (double **)calloc(p, sizeof(double *));
	if (*geno_bar == NULL)
	{
		printf("Memory allocation error: no enough memory\n");
		printLOG("Memory allocation error: no enough memory\n");
		exit(1);
	}
	for (int i = 0; i < p; i++)
	{
		(*geno_bar)[i] = (double *)calloc(2, sizeof(double));
		if ((*geno_bar)[i] == NULL)
		{
			printf("Memory allocation error: no enough memory\n");
			printLOG("Memory allocation error: no enough memory\n");
			exit(1);
		}
		(*geno_bar)[i][0] = 0;
		(*geno_bar)[i][1] = 0;
	}

	rewind(fp);

	j = 0; // column index
	ii = 0; // file index
	k = 0;

	while (!feof(fp)) {
		ii++;
		fscanf(fp, "%s\n", &filename_i);

		fp_i = fopen(filename_i, "r");
		if (fp_i == NULL)
		{
			fprintf(stderr, "can't open input file %s\n", filename_i);
			printLOG("can't open input file "+char2str(filename_i)+"\n");
			exit(1);
		}

		i = 0; //row index (individual index)
		while (!feof(fp_i)) {
			/* loop through and store the numbers into the array */

			if (j == 0)
			{
				//j = 0 means read class label y
				fscanf(fp_i, "%d", &tmp);

				if (tmp)// tmp=1 means case
					pheno[i] = true;

				j++;
			}
			else
			{
				fscanf(fp_i, "%d", &tmp);

				(*geno)[i][j-1] = tmp;
				if (pheno[i])
					(*geno_bar)[j-1][1] += tmp;
				else if(!pheno[i])
					(*geno_bar)[j-1][0] += tmp;
				else
				{
					printf("Incorrect phenotype\n");
					printLOG("Incorrect phenotype\n");
					exit(1);
				}

				j++; //column index
				if (j == (DataSize[ndataset + ii - 1] + 1))
				{
					j = 0;
					i++; // row index
				}

			}

			if (i >= n)
			{
				break;
			}
		}

		fclose(fp_i);
		k += DataSize[ndataset + ii - 1];
	}

	for (int i = 0; i < p; i++)
	{
		(*geno_bar)[i][0] = (*geno_bar)[i][0]/double(nctrl);
		(*geno_bar)[i][1] = (*geno_bar)[i][1]/double(ncase);
	}

	fclose(fp);

}

////////////////////////////////////////////////////////////////////// 
/*--------------------------------------------------------------------
From filename: "filenamelist.txt"; E.g. example_bt_tag_BOOST.txt
(1) pheno (bool* n): n samples; true: case, false: ctrl (This p is from BOOST.txt file)
(2) geno (BYTE** n*p): n samples, p snps; 0, 1, 2
(3) geno_bar (double** p*2): p snps; [0]: ctrl, [1]:case 
----------------------------------------------------------------------
From setname: E.g. "example_bt_tag.set"
(4) snpset (int** k*pi): k sets, pi snps in each set (i=0,..,k-1); 
	The index of the snp in snpname                   
----------------------------------------------------------------------
From mapname: E.g. "example_bt_tag.map"
(5) snpchr (int* p): p snps    (This p is from .map file)
(6) snpname (string* p): p snps
--------------------------------------------------------------------*/
//////////////////////////////////////////////////////////////////////
// vector<double> [0]: pgates; [1]: ptts; [2]: ptprod

double CalcRegionInter(RInside &R, string fout, vector<bool> &pheno, BYTE **geno, double **geno_bar, vector<int> &snpchr, vector<string> &snpname, bool skip_symm, int p, int n, int ncase, int nctrl, vector<int> &sA, vector<int> &sB, double myth_pgates, double myth_trun, int reps, bool flagperm, int max_cov_cnt)
{

	// Take a list of SNPs, or all SNPs 
	// Test these against either themselves, or all SNPs (vector<bool> sB)

	//  A     B
	//  ALL x ALL    skip e1>e2
	//  SET1 x ALL
	//  SET1 x SET1  skip e1>e2
	//  SET1 x SET2
	
	vector<double> zlist;
	vector<double> plist;
	vector<int> cov_index;	
	int cov_cnt = 0;

	//Note: Only need to store e1 and e2 (indices) <e1_1, e2_1, e1_2, e2_2, e1_3, e2_3, ...> in cov_index (length: 2*nInter)


	///////////////////////////////////////////////////////////
	// Begin iterating over pairs: SET x SET
	// Rewrite zlist, plist, and cov_index
	LDContrastTest(pheno, zlist, plist, cov_index, geno, geno_bar, sA, sB, skip_symm, p, n, ncase, nctrl);

	// Write results to output file
	ofstream EPI;
	EPI.open(fout.c_str(), ios::out);
	//printLOG("Writing epistasis pairwise results to [ "+fout+" ] \n");
	EPI.precision(4);
	EPI << setw(4) << "CHR1" << " "
	    << setw(pp_maxsnp) << "SNP1" << " "
	    << setw(4) << "CHR2" << " "
	    << setw(pp_maxsnp) << "SNP2" << " "
	    << setw(12) << "Z_LD" << " "
	    << setw(12) << "P_LD" << " "
	    << "\n";
	EPI.flush();
	for (int i = 0; i < zlist.size(); i++)
	{
		EPI << setw(4) << snpchr[cov_index[2*i]] << " " 
	        << setw(pp_maxsnp) << snpname[cov_index[2*i]] << " "
	        << setw(4) << snpchr[cov_index[2*i+1]] << " "
	        << setw(pp_maxsnp) << snpname[cov_index[2*i+1]]  << " "
	        << setw(12) << zlist[i]   << " "
	        << setw(12) << plist[i]  << " "
	        << "\n";
	  	EPI.flush();
	}

	EPI.close();
	fout += ".summary.cori";
	EPI.open(fout.c_str(), ios::out);
	EPI.clear();
	EPI.precision(12);
	/////////////////////////////////////////////////
	// Write covariance results
	if(cov_index.size() == 0)
	{
		printf("Error (cov mat): no interaction results\n");
		printLOG("Error (cov mat): no interaction results\n");
		exit(1);
	}

	cov_cnt = cov_index.size()/2;

	//////////////////////////////////////////////////////////////////
	/* // test for RInside

	// We can re-assign matrix to cori
	Rcpp::NumericMatrix corr_matrix1(6, 6);
	R["cori"] = corr_matrix1;
	R.parseEvalQ("print(cori)");

	Rcpp::NumericMatrix corr_matrix2(3, 3);
	R["cori"] = corr_matrix2;
	R.parseEvalQ("print(cori)");

	// cori in R is binded to corr_matrix2 (like a pointer), changing values in corr_matrix2 in C++ will also affect cori in R workspace!
	corr_matrix2(1,1) = 123;
	R.parseEvalQ("print(cori)");

	// end: test for RInside */
	//////////////////////////////////////////////////////////////////


	Rcpp::NumericMatrix corr_matrix(cov_cnt, cov_cnt);

	if (cov_cnt <= max_cov_cnt)
	{
		// Calculate covariance matrix
	    for(int i = 0;i < cov_cnt; i++)
	    { //snp-snp A1-A2
	    	for(int j = i; j < cov_cnt; j++)
	    	{ //snp-snp B1-B2
	    		int A1 = cov_index[2*i];
	    		int A2 = cov_index[2*i+1];
	    		int B1 = cov_index[2*j];
	    		int B2 = cov_index[2*j+1];

				// Calculate cov_case
				double LD_A1A2_case, LD_B1B2_case, LD_A1B1_case, LD_A2B2_case, LD_A1B2_case, LD_A2B1_case;
				double delta_4_case, delta_2_case, sigma_2_case, tau_2_case, cov_case;
				LD_A1A2_case = LD_B1B2_case = LD_A1B1_case = LD_A2B2_case = LD_A1B2_case = LD_A2B1_case = 0;
				delta_4_case = delta_2_case = sigma_2_case = tau_2_case = 0;

				// Calculate cov_control
				double LD_A1A2_control, LD_B1B2_control, LD_A1B1_control, LD_A2B2_control, LD_A1B2_control, LD_A2B1_control;
				double delta_4_control, delta_2_control, sigma_2_control, tau_2_control, cov_control;
				LD_A1A2_control = LD_B1B2_control = LD_A1B1_control = LD_A2B2_control = LD_A1B2_control = LD_A2B1_control = 0;
				delta_4_control = delta_2_control = sigma_2_control = tau_2_control = 0;

				for(int k = 0; k < n; k++)
				{
					if (pheno[k]) // case
					{
						LD_A1A2_case += (geno[k][A1] - geno_bar[A1][1])*(geno[k][A2] - geno_bar[A2][1]);
						LD_B1B2_case += (geno[k][B1] - geno_bar[B1][1])*(geno[k][B2] - geno_bar[B2][1]);
						LD_A1B1_case += (geno[k][A1] - geno_bar[A1][1])*(geno[k][B1] - geno_bar[B1][1]);
						LD_A2B2_case += (geno[k][A2] - geno_bar[A2][1])*(geno[k][B2] - geno_bar[B2][1]);
						LD_A1B2_case += (geno[k][A1] - geno_bar[A1][1])*(geno[k][B2] - geno_bar[B2][1]);
						LD_A2B1_case += (geno[k][A2] - geno_bar[A2][1])*(geno[k][B1] - geno_bar[B1][1]);

						delta_4_case += (geno[k][A1] - geno_bar[A1][1])*(geno[k][A2] - geno_bar[A2][1]) * (geno[k][B1] - geno_bar[B1][1])*(geno[k][B2] - geno_bar[B2][1]);
					}
					else if (!pheno[k]) // ctrl
					{
						LD_A1A2_control += (geno[k][A1] - geno_bar[A1][0])*(geno[k][A2] - geno_bar[A2][0]);
						LD_B1B2_control += (geno[k][B1] - geno_bar[B1][0])*(geno[k][B2] - geno_bar[B2][0]);
						LD_A1B1_control += (geno[k][A1] - geno_bar[A1][0])*(geno[k][B1] - geno_bar[B1][0]);
						LD_A2B2_control += (geno[k][A2] - geno_bar[A2][0])*(geno[k][B2] - geno_bar[B2][0]);
						LD_A1B2_control += (geno[k][A1] - geno_bar[A1][0])*(geno[k][B2] - geno_bar[B2][0]);
						LD_A2B1_control += (geno[k][A2] - geno_bar[A2][0])*(geno[k][B1] - geno_bar[B1][0]);

						delta_4_control += (geno[k][A1] - geno_bar[A1][0])*(geno[k][A2] - geno_bar[A2][0]) * (geno[k][B1] - geno_bar[B1][0])*(geno[k][B2] - geno_bar[B2][0]);
					}
					else
					{
						printf("Incorrect phenotype\n");
						printLOG("Incorrect phenotype\n");
						exit(1);
					}
					 
				}

				LD_A1A2_case = LD_A1A2_case/double(4*ncase);
				LD_B1B2_case = LD_B1B2_case/double(4*ncase);
				LD_A1B1_case = LD_A1B1_case/double(4*ncase);
				LD_A2B2_case = LD_A2B2_case/double(4*ncase);
				LD_A1B2_case = LD_A1B2_case/double(4*ncase);
				LD_A2B1_case = LD_A2B1_case/double(4*ncase);

				delta_4_case = delta_4_case/double(16*ncase);

				delta_2_case = LD_A1A2_case * LD_B1B2_case;
				sigma_2_case = LD_A1B1_case * LD_A2B2_case;
				tau_2_case = LD_A1B2_case * LD_A2B1_case;

				cov_case = delta_4_case - delta_2_case + (sigma_2_case + tau_2_case)/double(ncase-1);
				cov_case = cov_case/double(ncase);

				LD_A1A2_control = LD_A1A2_control/double(4*nctrl);
				LD_B1B2_control = LD_B1B2_control/double(4*nctrl);
				LD_A1B1_control = LD_A1B1_control/double(4*nctrl);
				LD_A2B2_control = LD_A2B2_control/double(4*nctrl);
				LD_A1B2_control = LD_A1B2_control/double(4*nctrl);
				LD_A2B1_control = LD_A2B1_control/double(4*nctrl);

				delta_4_control = delta_4_control/double(16*nctrl);

				delta_2_control = LD_A1A2_control * LD_B1B2_control;
				sigma_2_control = LD_A1B1_control * LD_A2B2_control;
				tau_2_control = LD_A1B2_control * LD_A2B1_control;

				cov_control = delta_4_control - delta_2_control + (sigma_2_control + tau_2_control)/double(nctrl-1);
				cov_control = cov_control/double(nctrl);

				if(i < j )
					corr_matrix(i, j) = corr_matrix(j, i) = cov_case + cov_control;
				else
					corr_matrix(i, i) = cov_case + cov_control;
		 
			} //end for(int j = i; j < cov_cnt; j++)
	   	} //end for(int i = 0;i < cov_cnt; i++)

	   	// Tranform covariance matrix into correlation matrix
	   	for(int i = 0; i < cov_cnt; i++)
	   	{
			for(int j = i+1; j < cov_cnt; j++)
                // Transform cov(\Delta LD_A, \Delta LD_B) to cov(T_A, T_B)
				corr_matrix(i, j) = corr_matrix(j, i)= corr_matrix(i, j)/sqrt(corr_matrix(i, i)*corr_matrix(j, j));
	   	}
	   	for(int i = 0; i < cov_cnt; i++)
	   		corr_matrix(i, i) = 1;

	    // Write the correlation matrix to result file
	    for(int i = 0; i < cov_cnt; i++)
	    {
			for(int j = 0; j < cov_cnt - 1; j++)
				EPI << setw(12) << corr_matrix(i, j)<< "\t";

			// Write the last one individually cuz need \n rather than \t at the end
			EPI << setw(12) << corr_matrix(i, cov_cnt-1) <<endl;
			EPI.flush();
	    }
	 
	} //end if (cov_int <= max_cov_cnt) 
	else
	{
		cout << "WARNING: covariance matrix is too large (max: " << max_cov_cnt  << "by " <<  max_cov_cnt << ")." << endl;
		cout << "You will not be able to run gene_based_interaction methods." << endl;
	}
	
	EPI.close();

	// Calc and return pmin
	clock_t st, ed;
	st = clock();
	R["numpv"] = plist.size();
	R["cori"] = corr_matrix;
	R["minpv"] = *min_element(plist.begin(), plist.end());

	double pmin = R.parseEval("1-pmvnorm(lower=qnorm(minpv/2),upper=-qnorm(minpv/2),mean=rep(0, numpv),corr=cori)");
	ed = clock();
	printf("cputime for calling R fucntions pmvnorm: %f seconds.\n", (double)(ed - st)/ CLOCKS_PER_SEC);
	return pmin;

}


// Clear and ReWrite zlist, plist, and cov_index
void LDContrastTest(vector<bool> &pheno, vector<double> &zlist, vector<double> &plist, vector<int> &cov_index, BYTE **geno, double **geno_bar, vector<int> &sA, vector<int> &sB, bool skip_symm, int p, int n, int ncase, int nctrl)
{
	int ii = 0;
	double z_ld, ld_aff, ld_unf, v_ld_aff, v_ld_unf;

	zlist.clear();
	plist.clear();
	cov_index.clear();
	
	// Begin iterating over pairs: SET x SET
	for (int e1 = 0; e1 < p; e1++)
	{
		if(sA[e1])
		{
			for(int e2 = 0; e2 < p; e2++)
			{
				///////////////////////////////////////////////
				// Skip this test under certain conditions

				// The SNP not in the set
				if (!sB[e2]) continue; 

		        // We've already performed this test
		        if (e1>=e2 && skip_symm) continue;

		        // Same SNP 
		        if (e1==e2) continue;

		        ///////////////////////////////////////////////
		        // Copy epi_Sen.cpp here
		        // The modified code based on pheno and geno can refer to epi_Sen.cpp

		        // coding scheme: 0, 1, 2 -> Need to transform to 0, 1/2, 1
		        z_ld = ld_aff = ld_unf = v_ld_aff = v_ld_unf = 0;

		        for (int i = 0; i < n; i++)
		        {
		        	if (pheno[i])
		        		ld_aff += (geno[i][e1] - geno_bar[e1][1])*(geno[i][e2] - geno_bar[e2][1]);
		        	else
		        		ld_unf += (geno[i][e1] - geno_bar[e1][0])*(geno[i][e2] - geno_bar[e2][0]);
		        	
		        }

		        ld_aff = ld_aff/double(4*ncase);
		        ld_unf = ld_unf/double(4*nctrl);

		        v_ld_aff = 1/double(64*ncase)*geno_bar[e1][1]*(2-geno_bar[e1][1])*geno_bar[e2][1]*(2-geno_bar[e2][1]) + 1/double(4*ncase)*(1-geno_bar[e1][1])*(1-geno_bar[e2][1])*ld_aff - 1/double(4*ncase)*pow(ld_aff, 2);
		        v_ld_unf = 1/double(64*nctrl)*geno_bar[e1][0]*(2-geno_bar[e1][0])*geno_bar[e2][0]*(2-geno_bar[e2][0]) + 1/double(4*nctrl)*(1-geno_bar[e1][0])*(1-geno_bar[e2][0])*ld_unf - 1/double(4*nctrl)*pow(ld_unf, 2);

		        z_ld = (ld_aff - ld_unf)/sqrt(v_ld_aff + v_ld_unf);

		        cov_index.push_back(e1);
		        cov_index.push_back(e2);
		        zlist.push_back(z_ld);
		        plist.push_back(normdist(-fabs(z_ld)) * 2);

				if ((ii+1)%10000==0)
					printf("SNP pair %d.\n", ii+1);

				ii ++;

			} // end e2
		} // end if(sA[e1])
	} // end e1
}
