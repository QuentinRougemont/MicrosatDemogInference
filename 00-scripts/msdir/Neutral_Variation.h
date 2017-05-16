#ifndef NEUTRALVARIATION_H
#define NEUTRALVARIATION_H

#include <string>
#include <iostream>
#include <vector>
#include <boost/regex.hpp>
#include <map>
#include <sstream>
#include <ctime>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include "ran.h"
#include "nr3.h"
#include "Genealogy.h"

using std::string;
using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::map;
using std::time;

template <class Mutational_Model>
class Neutral_Variation {

public:

	void mutate (Ran &ran) 
	{
		for (int i = num_groups; i>0; --i) {

			const double &brlen = genealogy.get_branch_lengths(i);
			int k = 0; // stores number of mutations along branch
			double rannum = ran.doub();
			double lambda = brlen*theta;
			double cumulative_P = 0.;
			while (true) {
				double factorial; 
				double lambda_to_theK;
				if (k==0) {
					factorial = 1;
					lambda_to_theK = 1;
				}else {
					factorial = 1;
					lambda_to_theK =1;
					for (int x =k; x>0; --x) {
						factorial *= x;
						lambda_to_theK *= lambda;
					}
				}
				double toadd = exp(-lambda);
				toadd *= lambda_to_theK;
				toadd /= factorial;
				
				cumulative_P += toadd; 
				if (rannum <= cumulative_P)
					break;
				else 
					++k;
			}
			Ran &rand = ran;
			int change = model(k, rand, genealogy, gp);
			if (change == -999) { // used in the case of an ISM
				const vector<int> &vecky = genealogy.get_group_makeup(genealogy.get_mutated_group()); // retrieve the samples descended from the mutated branch
				for (vector<int>::const_iterator iter = vecky.begin(); iter!=vecky.end(); ++iter)
					samples[*iter] = 1;			
				break; // using the ISM and don't need to loop anymore
			} else if (change == -899) { // used in the case of the first use of SNPSTR() call
				for (int q=1; q<=num_samples; ++q) // set all SNPs in the SNPSTR to ancestral state "0"
					snp_components[q] = 0; 
				const vector<int> &vecky = genealogy.get_group_makeup(genealogy.get_mutated_group());
				for (vector<int>::const_iterator iter = vecky.begin(); iter!=vecky.end(); ++iter)
					snp_components[*iter] = 1;   
				int special_change = model(k, rand, genealogy, gp); // use the SNPSTR() operator again to get the STRP mutation started			
				const vector<int> &vecky2 = genealogy.get_group_makeup(i);
				for (vector<int>::const_iterator iter=vecky2.begin(); iter != vecky2.end(); ++iter) 	
					samples[*iter] += special_change;
			} else if (change == -799) { //used in the case of SNPH()
				string toadd("");
				string toaddNON("");
				for (int snp =0; snp<k; ++snp) {
					toadd.append("1");
					toaddNON.append("0");
				}
				const vector<int> &vecky = genealogy.get_group_makeup(i);
				for (int g=1; g<=num_samples; ++g) {
					int check=0;
					for (vector<int>::const_iterator iter=vecky.begin(); iter!=vecky.end(); ++iter) {
						if (g==*iter) {
							check =1;							
							break;
						}
					}
					if (check ==0) 
						snphs[g].append(toaddNON);
					else 
						snphs[g].append(toadd);
				}
			} else {
				const vector<int> &vecky = genealogy.get_group_makeup(i);
				for (vector<int>::const_iterator iter=vecky.begin(); iter != vecky.end(); ++iter){
					samples[*iter] += change;
				}
			}		
		}
	}

	double calchet(int samp) 
	{
		map<int, double> frequencies;
		for (int q=0; q<samp; ++q) 
			frequencies[samples[q]]++;
		double exp_hom = 0.;
		for (map<int, double>::iterator iter=frequencies.begin(); iter!=frequencies.end(); ++iter) 
			exp_hom += (iter->second/static_cast<double>(samp)) * (iter->second/static_cast<double>(samp));
		return 1-exp_hom;
	}

	void create_compound_alleles()  // if using snpstr mutation model, must call this after mutate to modify samples by snp phenotype
	{
		for (int q = 1; q<=num_samples; ++q) 
			if (snp_components[q] == 1) 
				samples[q] *= 100;
	}	


	const vector<int> get_dataset() 
	{
		vector<int> vecky;
		vecky.reserve(num_samples);
		for (int x=1; x<=num_samples; ++x) 
			vecky.push_back(samples[x]);
		return vecky;
	}

	const vector<string> get_snphs()
	{
		vector<string> vecky;
		for (int x=1; x<=num_samples; ++x) 
			vecky.push_back(snphs[x]);
		return vecky;
	}

	const vector<int> get_snp_components()  // used with the SNPSTR mutational model to return the SNP components of the SNPSTR
						// the STRP component is obtained using the usual get_dataset() method above 
	{
		vector<int> vecky;
		vecky.reserve(num_samples);
		for (int x=1; x<=num_samples; ++x) 
			vecky.push_back(snp_components[x]);
		return vecky;
	}


	const int get_num_alleles()		 // returns number of alleles at the locus
	{
		vector<int> unique_alleles;
		for (map<int, int>::iterator iter=samples.begin(); iter!=samples.end(); ++iter) {
			int test = 1;
			for (vector<int>::iterator iter2=unique_alleles.begin(); iter2!=unique_alleles.end(); ++iter2) {
				if (iter->second == *iter2) {
					test = 0;
					break;
				}				
			}
			if (test==1) 
				unique_alleles.push_back(iter->second);
		}
		vector<int>::size_type na = unique_alleles.size();
		int num_alleles = na;
		return num_alleles;
	} 

	inline const int get_range()			{return range;} // returns allele size range

	Neutral_Variation (Genealogy &ggenealogy, Mutational_Model &mmodel, const double ttheta, vector<double> &gptable): genealogy(ggenealogy), model(mmodel), theta(ttheta), gp(gptable) 
	{
		num_groups = genealogy.get_group_number();
		num_samples = genealogy.get_sample_number();
		for (int q=1; q<=num_samples; ++q)
			samples[q] = 0;
		range=-1;
	}

private:
	Genealogy &genealogy;
	Mutational_Model &model;
	map<int, int> samples;
	map<int, int> snp_components; // [SNPSTR model] a separte array to store the linked SNP variation
	vector<double> &gp;
	map<int, string> snphs; //snp haplotype samples
	int num_groups;
	int num_samples;
	int range;
	int num_alleles;
	double theta;  
};

struct SMM { // stepwise mutation model
	int operator () (const int k, Ran &random, Genealogy &genealogy, vector<double> &gp) {
		int change=0;
		double rannum2;
		if (k > 0) {
			for (int x=0; x<k; ++x) {
				rannum2 = random.doub();
				if (rannum2 >= 0.5) 
					change += 1;
				else 
					change -= 1;
			}
		}
		return change;
	}
	SMM() {	}
};

struct SNPH {
	int operator () (const int k, Ran &random, Genealogy & genealogy, vector<double> &gp) {
		int change=-799;
		return change;
	}
	SNPH() { }
};

struct GSM {  // generalized stepwise model
	int operator () (const int k, Ran &random, Genealogy &genealogy, vector<double> &gp) {
		int change=0;
		int change_modifier = 0;
		double rannum2;
		if (k > 0) {
			for (int x=0; x<k; ++x) {
				rannum2 = random.doub();
				if (rannum2 <= 0.5) 
					change += 1;
				else 
					change -= 1;
			}
			rannum2 = random.doub();
			for (int x=0; x<20000; ++x) {
				if (rannum2 <= gp[x]) {
					change_modifier += x+1;	
					break;
				}
			}
			change *= change_modifier;
		}
		return change;
	}
	GSM() {	}
};

struct DSM { // diffuse stepwise model (an infinite alleles model for STRPs)
	int operator () (const int k, Ran &random, Genealogy &genealogy, vector<double> &gp) {
		int change=0;
		double rannum2 = random.doub();
		if (k > 0) {
			for (int x=0; x<k; ++x) {
				rannum2 = random.doub();
				if (rannum2 <= 0.5) 
					change += 1;
				else 
					change -= 1;
			}
			rannum2 = random.doub();
			double modifier = 1000*rannum2;
			change *= static_cast<int>(modifier);
		}
		return change;
	}
	DSM() {	}
};

struct ISM { 
	int operator () (const int k, Ran &random, Genealogy &genealogy, vector<double> &gp) {
		const int change = -999;
		int group_to_change;
		double rannum2 = random.doub();
		double brlen_total = 0.;
		const int groups = genealogy.get_group_number(); 
		for (int i=groups; i>0; --i) {
			const double &brlen = genealogy.get_branch_lengths(i);
			brlen_total += brlen;
		}
		double local_total = 0.;
		for (int i=groups; i>0; --i) {
			const double &brlen = genealogy.get_branch_lengths(i);
			if ( (local_total += brlen/brlen_total) >=rannum2 ){
				group_to_change = i;
				break;
			}
		}	
		genealogy.mark_mutated_group(group_to_change); // set the mutated_group variable of the Genealogy &genealogy object		
		return change;
 	}	
	ISM() {	}
};

struct SNPSTR { // produce a SNP and SMM-modeled STRP on the same genealogy
	int operator () (const int k, Ran &random, Genealogy &genealogy, vector<double> &gp) {
		bool snpstr_status = genealogy.get_snpstr_status();
		int group_to_change;
		if (snpstr_status == false) { // establish SNP variation
			int change = -899;
			while (true) {
				int change = -899;
				genealogy.trip_snpstr_switch();
				double rannum2 = random.doub();
				double brlen_total = 0.;
				const int groups = genealogy.get_group_number(); 
				for (int i=groups; i>0; --i) {
					const double &brlen = genealogy.get_branch_lengths(i);
					brlen_total += brlen;
				}
				double local_total = 0.;
				for (int i=groups; i>0; --i) {
					const double &brlen = genealogy.get_branch_lengths(i);
					if ( (local_total += brlen/brlen_total) >=rannum2 ){
						group_to_change = i;
						break;
					}
				}
			
				vector<int>::size_type na;
				const vector<int> &vecky = genealogy.get_group_makeup(group_to_change);
				na = vecky.size();
				double group_size = na;
				int ns = genealogy.get_sample_number();
				if (na/static_cast<double>(ns) > 0.) // set test to ">0.0", if you don't want SNP ascertainment bias 
					break;
				else 
					continue; // try again; this variant found in <10% of samples
				
			}
			genealogy.mark_mutated_group(group_to_change);		
			return change;
		} else { // not the first time the SNPSTR() operator was used, so continue work on the microsatellite locus
			int change=0;
			double rannum2;
			if (k > 0) {
				for (int x=0; x<k; ++x) {
					rannum2 = random.doub();
					if (rannum2 >= 0.5) 
						change += 1;
					else 
						change -= 1;
				}
			}
			return change;
		}
	}
	SNPSTR() {}
};

#endif
