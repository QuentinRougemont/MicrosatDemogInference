#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/regex.hpp>
#include <math.h>
#include "Genealogy.h"
#include "Neutral_Variation.h"
#include "Dataset.h"


using std::cout;
using std::endl;
using std::cin;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::map;


int main () {


	cout << endl << "MARKSIM **************" << endl;

	map<string, string> parameters;
	string line; 

	ifstream parameterfile;
	string infile("params");
	parameterfile.open(infile.c_str());

	vector<int> island_sizes;

	while (parameterfile >> line) {		
		string name = line;
		string toadd;
		if (name == "ISLANDS") {
			int isl;
			while (parameterfile >> isl) {
				island_sizes.push_back(isl/2);
			}
		} else {
			while (parameterfile >> line) {
				if (line != "//") 
					toadd += (line+' ');
				else 
					break;
			}	
		}
		parameters[name] = toadd;
	}

	parameterfile.close();

	cout << "Input Parameter Values" << endl;
	cout << "**********************" << endl;
	for (map<string, string>::iterator iter=parameters.begin(); iter!=parameters.end(); ++iter) 
		cout << iter->first << "\t" << iter->second << endl;
	cout << "# DIPLOIDS/ISLAND\t";
	for (vector<int>::iterator iter=island_sizes.begin(); iter!=island_sizes.end(); ++iter) 
		cout << *iter << "\t";
	cout << "\n\n" << endl;

	stringstream ss(stringstream::in | stringstream::out);
	int format, ds, loci, monomorphs, model, asc, verbose, ascsamp;
	double theta, p, aschet;
	string prefix;
 
	ss << parameters["VERBOSE"];
	ss >> verbose;
	ss << parameters["FORMAT"];
	ss >> format;
	ss << parameters["THETA"];
	ss >> theta;
	ss << parameters["DATASETS"];
	ss >> ds;
	ss << parameters["LOCI"];
	ss >> loci;
	ss << parameters["MONOMORPHS"];
	ss >> monomorphs;
	ss << parameters["MODEL"];
	ss >> model;
	ss << parameters["P"];
	ss >> p;
	ss << parameters["FILEPREFIX"];
	ss >> prefix;
	ss << parameters["ASCERTAINMENT"];
	ss >> asc;
	ss << parameters["ASCSAMPSIZE"];
	ss >> ascsamp;
	ss << parameters["ASCHET"];
	ss >> aschet;

	//populate stepsize probability table
	vector<double> gp;
	double cumulative = 0;
	for (int i=1; i<20000; i++) {
		double toadd = pow(1-p, i-1)*p;
		cumulative += toadd;
		gp.push_back(cumulative);
		if (cumulative > 0.99999) {
			gp[i] = 1.;
			break;
		}
	}


	double variance=0.;
	int varcount=0;
	for (int d=0; d<ds; d++) { 

		Dataset dataset(format, loci);

		if (verbose || verbose == 2) {cout << "dataset " << d << endl;}

		for (int loc=0; loc<loci; loc++) {
			
			if (verbose == 2) {cout << "\t...locus " << loc << endl; }

			// run MS
			string ms_command(parameters["MSCOMMAND"]);
			system(ms_command.c_str());
		
			ifstream infile;
			int seed = time(0);
			Ran ran(seed); // random number generator object
			int sample_size;
		
			string in("ms_outfile");
			infile.open(in.c_str());  
			string line;
			while (infile >> line) {
				if (line == "//") {
					infile >> line; // next line has the tree
					break;
				}
			}
			infile.close();
			infile.clear(); 	
			Genealogy genealogy(line);
			Genealogy &gen_ref = genealogy;
			sample_size = genealogy.get_sample_number();
			vector<double> &gpref = gp;
			Ran &rand = ran;

			switch (model) {
			
				case 0:
					{
					//*************** STEPWISE MUTATION MODEL
					SMM smm;
					Neutral_Variation<SMM> data(gen_ref, smm, theta, gpref);
					data.mutate(rand);
					if (asc==1) {
						double het=data.calchet(ascsamp); 
						if (het < aschet) {
							--loc; 
							break;
						}
					}
					int na = data.get_num_alleles();
					if (monomorphs==0 && na==1) { // no variation and "no monomorphs" specified
						--loc;
						cout << "monomorph" << endl;
						break;
					}
					Neutral_Variation<SMM> &dataref =data;
					dataset.addlocus(dataref, loc);
					break;
					}	
				case 1:
					{	
					//*************** GENERALIZED STEPWISE MODEL
					GSM gsm;
					Neutral_Variation<GSM> data(gen_ref, gsm, theta, gpref);
					data.mutate(rand);
					if (asc==1) {
						double het=data.calchet(ascsamp);
						if (het < aschet) {
							--loc; 
							break;
						}
					}
					int na = data.get_num_alleles();
					if (monomorphs==0 && na==1) { // no variation and "no monomorphs" specified
						--loc;
						cout << "monomorph" << endl;
						break;
					}
					Neutral_Variation<GSM> &dataref =data;
					dataset.addlocus(dataref, loc);
					break;
					}
				case 2:	
					{
					//*************** INFINITE ALLELES MODEL FOR MICROSATELLITES
					DSM dsm;
					Neutral_Variation<DSM> data(gen_ref, dsm, theta, gpref);
					data.mutate(rand);
					if (asc==1) {
						double het=data.calchet(ascsamp);
						if (het < aschet) {
							--loc; 
							break;
						}
					}
					int na = data.get_num_alleles();
					if (monomorphs==0 && na==1) { // no variation and "no monomorphs" specified
						--loc;
						cout << "monomorph" << endl;
						break;
					}
					Neutral_Variation<DSM> &dataref =data;
					dataset.addlocus(dataref, loc);
					break;
					}
				case 3:
					{
					//*************** SNPs under the infinite sites model
					ISM ism;
					theta=10; // high theta ensures there will be >0 SNPs to work with
					Neutral_Variation<ISM> data(gen_ref, ism, theta, gpref);
					data.mutate(rand);
					if (asc==1) {
						double het=data.calchet(ascsamp);
						if (het < aschet) {
							--loc; 
							break;
						}
					}
					Neutral_Variation<ISM> &dataref =data;
					dataset.addlocus(dataref, loc);
					break;
					}
				case 4:
					{
					//*************** SNPSTRs
					SNPSTR snpstr;
					Neutral_Variation<SNPSTR> data(gen_ref, snpstr, theta, gpref);
					data.mutate(rand);
					if (asc==1) {
						double het=data.calchet(ascsamp);
						if (het < aschet) {
							--loc; 
							break;
						}
					}
					data.create_compound_alleles();
					Neutral_Variation<SNPSTR> &dataref =data;
					dataset.addlocus(dataref, loc);
					break;			
					}
				case 5: 
					{
					//**************** SNP haplotypes
					SNPH snph;
					Neutral_Variation<SNPH> data(gen_ref, snph, theta, gpref);
					data.mutate(rand);
					Neutral_Variation<SNPH> &dataref = data;
					dataset.addlocus(dataref, loc);
					}
			}

		} // end locus loop

		dataset.printfile(prefix, d, format, loci, island_sizes);

	} // end data set loop
	return 0;

}
	
		
