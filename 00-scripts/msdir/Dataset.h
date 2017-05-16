#ifndef DATASET_H
#define DATASET_H


class Dataset {


public:
	void printfile(string prefix, int rep, int format, int loci, vector<int> islands) 
	{

		stringstream ss(stringstream::in | stringstream::out);
		string filenum, finalname;
		ss << rep;
		ss >> filenum;
		finalname=prefix+filenum;
		ofstream outstar;
		string out(finalname);
		outstar.open(out.c_str());

		this->makesamples();
		
		switch (format) {

			case 0: // STRUCTURE format
				{
				int count=0;
				int samplenum=0;
				int pop = 1;
				int check_value = islands[0] + 1;
				for (map<int, vector<int> >::iterator iter=samples.begin(); iter!=samples.end(); ++iter) {
					if (count%2==0)
						++samplenum;
					if (samplenum == check_value) {
						++pop;
						check_value += islands[pop-1];
					}	
					outstar << samplenum << " " << pop << " 0 "; // !!!! change this line and recompile to tweak structure formatting
					for (vector<int>::iterator iter2= (iter->second).begin(); iter2!= (iter->second).end(); ++iter2) 
						outstar << *iter2 << " ";
					outstar << endl;
					++count;					
				}
				break;
				}
		
			case 1: // SMARTPCA format   SNP data ONLY!!!!
				{
				//.geno file
				int num_individuals = 0;
				int testflag =0;
				for (map<int, vector<int> >::iterator iter=chromosomes.begin(); iter!=chromosomes.end(); ++iter) {
					for (vector<int>::iterator iter2=(iter->second).begin(); iter2!= (iter->second).end(); ++iter2) {
						int genotype_count = 0;
						if (testflag == 0) 
							num_individuals+=1;
						genotype_count += *iter2;
						++iter2;
						genotype_count += *iter2;
						outstar << genotype_count;
					}
					++testflag;
					outstar << endl;
				}
				
				string system_commandpt1 = "mv ";
				string system_commandpt2 = " ";
				string system_commandpt3 = ".geno";
				string system_command = system_commandpt1 + finalname + system_commandpt2 + finalname + system_commandpt3;
				system(system_command.c_str());
				//.snp file
				int physical_position = 0;
				ofstream snpout;
				string snpfile = finalname + ".snp";
				snpout.open(snpfile.c_str());
				for (int i=0; i<loci; ++i) {
					physical_position += 50;
					snpout << "snp" << i << "\t" << "11 " << "\t" << "0.0" << "\t" << physical_position << "\t" <<  endl;
				}
				snpout.close();
				//.ind file
				ofstream indout;
				string indfile = finalname + ".ind";
				indout.open(indfile.c_str());
				int island_counter = 1;
				int check_value = islands[0];
				for (int i=0; i<num_individuals; ++i)  {
					if (i == check_value) {
						++island_counter;
						check_value += islands[island_counter-1];
					}	
					indout << "ind" << i << " U pop" << island_counter << endl;
				} 
				indout.close();
				break;
				}
			

			case 2: // whitespace-delimited file

				{
				for (map<int, vector<int> >::iterator iter=samples.begin(); iter!=samples.end(); ++iter) {
					for (vector<int>::iterator iter2= (iter->second).begin(); iter2!= (iter->second).end(); ++iter2) 
						outstar << *iter2 << " ";
					outstar << endl;					
				}
				break;	
				}

			case 3: // Arlequin format (for microsatellites !!!)
				{
				outstar << "#file written using simulation program MARKSIM" << endl << endl;
				outstar << "[Profile]" << endl;
				outstar << "\tTitle=\"A series of simulated samples\"" << endl;
				outstar << "\tNbSamples=2" << endl << endl;
				outstar << "\tGenotypicData=1" << endl;
				outstar << "\tGameticPhase=0" << endl;
				outstar << "\tRecessiveData=0" << endl;
				outstar << "\tDataType=MICROSAT" << endl;
				outstar << "\tLocusSeparator=TAB" << endl;
				outstar << "\tMissingData='?'" << endl << endl;
				outstar << "[Data]" << endl;
				outstar << "\t[[Samples]]" << endl;

				int count=0;
				int samplenum=0;
				int pop = 1;
				int check_value = islands[0] + 1;
				int numpops = (int) islands.size();

				outstar << "\t\tSampleName=\"Sample " << pop << "\"" << endl;
				outstar << "\t\tSampleSize=" << islands[pop-1] << endl;
				outstar << "\t\tSampleData= {" << endl;

				for (map<int, vector<int> >::iterator iter=samples.begin(); iter!=samples.end(); ++iter) {
					if (count%2==0)
						++samplenum;
					if (samplenum == check_value) {
						++pop;
						check_value += islands[pop-1];
						outstar << endl << "}" << endl;
						outstar << "\t\tSampleName=\"Sample " << pop << "\"" << endl;
						outstar << "\t\tSampleSize=" << islands[pop-1] << endl;
						outstar << "\t\tSampleData= {" << endl;
					}	
					if (count%2 == 0) outstar << pop << "_" << samplenum << "\t" << "1" << "\t"; 
					else outstar << "\t\t";
					for (vector<int>::iterator iter2= (iter->second).begin(); iter2!= (iter->second).end(); ++iter2) 
						outstar << *iter2 + 30 << "\t";
					outstar << endl;
					++count;	
				}

				outstar << endl << "}" << endl;
				outstar << "\n[[Structure]]" << endl << endl;
				outstar << "\tStructureName=\"Simulated data\"" << endl;
				outstar << "\tNbGroups=1" << endl;
				outstar << "\tGroup={" << endl;
				for (int q=1; q<=numpops; q++) outstar << "\t\"Sample " << q << "\""<< endl;
				outstar << "\t}" << endl;

				break;
				}


				case 4: // Arlequin format (for SNPs !!!)
				{
				outstar << "#file written using simulation program MARKSIM" << endl << endl;
				outstar << "[Profile]" << endl;
				outstar << "\tTitle=\"A series of simulated samples\"" << endl;
				outstar << "\tNbSamples=2" << endl << endl;
				outstar << "\tGenotypicData=1" << endl;
				outstar << "\tGameticPhase=0" << endl;
				outstar << "\tRecessiveData=0" << endl;
				outstar << "\tDataType=DNA" << endl;
				outstar << "\tLocusSeparator=TAB" << endl;
				outstar << "\tMissingData='?'" << endl << endl;
				outstar << "[Data]" << endl;
				outstar << "\t[[Samples]]" << endl;

				int count=0;
				int samplenum=0;
				int pop = 1;
				int check_value = islands[0] + 1;
				int numpops = (int) islands.size();

				outstar << "\t\tSampleName=\"Sample " << pop << "\"" << endl;
				outstar << "\t\tSampleSize=" << islands[pop-1] << endl;
				outstar << "\t\tSampleData= {" << endl;

				for (map<int, vector<int> >::iterator iter=samples.begin(); iter!=samples.end(); ++iter) {
					if (count%2==0)
						++samplenum;
					if (samplenum == check_value) {
						++pop;
						check_value += islands[pop-1];
						outstar << endl << "}" << endl;
						outstar << "\t\tSampleName=\"Sample " << pop << "\"" << endl;
						outstar << "\t\tSampleSize=" << islands[pop-1] << endl;
						outstar << "\t\tSampleData= {" << endl;
					}	
					if (count%2 == 0) outstar << pop << "_" << samplenum << "\t" << "1" << "\t"; 
					else outstar << "\t\t";
					for (vector<int>::iterator iter2= (iter->second).begin(); iter2!= (iter->second).end(); ++iter2) 
						if (*iter2 == 0) outstar << "A\t";
						else outstar << "G\t";
					outstar << endl;
					++count;	
				}

				outstar << endl << "}" << endl;
				outstar << "\n[[Structure]]" << endl << endl;
				outstar << "\tStructureName=\"Simulated data\"" << endl;
				outstar << "\tNbGroups=1" << endl;
				outstar << "\tGroup={" << endl;
				for (int q=1; q<=numpops; q++) outstar << "\t\"Sample " << q << "\""<< endl;
				outstar << "\t}" << endl;

				break;
				}

		}

		outstar.close();
	}

	void makesamples() 
	{
		for (map<int, vector<int> >::iterator iter=chromosomes.begin(); iter != chromosomes.end(); ++iter) {
			int count = -1;
			for(vector<int>::iterator iter2=(iter->second).begin(); iter2 != (iter->second).end(); ++iter2) {
				++count;
				samples[count].push_back(*iter2);
			}
		}

	}

	inline void addlocus (Neutral_Variation<SMM> &data, int locus) {chromosomes[locus] = data.get_dataset();}
	inline void addlocus (Neutral_Variation<GSM> &data, int locus) {chromosomes[locus] = data.get_dataset();}
	inline void addlocus (Neutral_Variation<ISM> &data, int locus) {chromosomes[locus] = data.get_dataset();}
	inline void addlocus (Neutral_Variation<DSM> &data, int locus) {chromosomes[locus] = data.get_dataset();}
	inline void addlocus (Neutral_Variation<SNPSTR> &data, int locus) {chromosomes[locus] = data.get_dataset();}
	void addlocus (Neutral_Variation<SNPH> &data, int locus) 
	{
		vector<string> haplotypes = data.get_snphs();
		int label = 1;
		map<string, int> transform;
		for (vector<string>::iterator iter=haplotypes.begin(); iter!= haplotypes.end(); ++iter) {
			int check = 1;
			for (map<string, int>::iterator iter2=transform.begin(); iter2!=transform.end(); ++iter2) {
				if (*iter == iter2->first) {
					chromosomes[locus].push_back(transform[*iter]);
					check = 0;
					break;
				}
			}
			if (check == 1) {
				transform[*iter] = label;
				++label;
				chromosomes[locus].push_back(transform[*iter]);
			}
		}
			
	}

	Dataset (int fformat, int lloci): format(fformat), loci(lloci) 
	{
		;
	}

private:
	int format;
	int loci;	
	map<int, vector<int> > chromosomes; // locus --> all sampled alleles at the locus
	map<int, vector<int> > samples; // individual --> alleles of the individual at each locus

};

#endif
