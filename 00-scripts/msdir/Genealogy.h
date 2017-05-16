#ifndef GENEALOGY_H
#define GENEALOGY_H

#include <string>
#include <iostream>
#include <vector>
#include <boost/regex.hpp>
#include <map>
#include <sstream>

using std::string;
using std::vector;
using std::cout;
using std::cin;
using std::endl;
using std::map;

class Genealogy {

public:
	void initialize_samples ()  // initializes samples by storing external branch length and creating the first "groups"
	{
		boost::regex re;
		re.assign("\\D(\\d+):(\\d*.\\d\\d\\d)");
		boost::smatch matches;	
		string::const_iterator iterBeg = newick.begin();
		string::const_iterator iterEnd = newick.end();
		while (boost::regex_search(iterBeg, iterEnd, matches, re)) {
	
			string s1(matches[1].first, matches[1].second);
			string s2(matches[2].first, matches[2].second);
			string bind_string = s1+" " + s2;
			std::istringstream converter(bind_string); // convert strings to numerical types
			int sample; double blen;
			converter >> sample >> blen;
			branch_lengths[sample] = blen; // populate branch_lengths map
			vector<int> to_pass(1, sample); 
			groups[sample]=to_pass; // populate groups map
			samples[sample] = 0; // initialize each sample allele size to 0
			++num_samples;
			++num_groups;
			
			iterBeg = matches[0].second;

		}

	}
	void group_samples () // adds internal branch lengths to map "branch_lengths" and samples to map "groups"
	{
		boost::regex re("(\\d+):\\d*.\\d\\d\\d");
		boost::smatch matches;
		string::const_iterator iterBeg = newick.begin();
		string::const_iterator iterEnd = newick.end();
		vector<string> svec;
		// remove external branch lengths by "hand" building string instead of using boost::regex_replace: SIGNIFICANT increase in speed
		while (boost::regex_search(iterBeg, iterEnd, matches, re)) { 
			string replacement(matches[1].first, matches[1].second);
			string pre_match(iterBeg, matches[0].first);			
			iterBeg=matches[0].second;   // advance temp_newick iterator to facilitate searching of a smaller string for matches
			string piece= pre_match+replacement;
			svec.push_back(piece);
		}
		string final_piece(iterBeg, iterEnd); // add the trailing substring
		svec.push_back(final_piece);

		newick = "";
		for (vector<string>::iterator iter=svec.begin(); iter != svec.end(); ++iter) 
			newick += *iter;
		// create more and more inclusive groups
		re.assign("\\((\\d+),(\\d+)\\):(\\d*.\\d\\d\\d)");
		while (boost::regex_search(newick, matches, re)) {  // changed newick to modified_string

			++num_groups;			
			string s1(matches[1].first, matches[1].second); // first sample/group
			string s2(matches[2].first, matches[2].second); // second sample/group
			string s3(matches[3].first, matches[3].second); // branch length
			std::stringstream out;
			out << num_groups;
			string replacement = out.str();	 // new group # as a string
			string bind_string = s1 + " " + s2 + " " + s3; 
			std::istringstream converter(bind_string); // convert strings to numerical types
			int samp1, samp2; double blen;
			converter >> samp1 >>  samp2 >> blen;

			vector<int> to_pass = groups[samp1];
			vector<int> moar = groups[samp2];
			vector<int>::iterator iter = moar.begin();
			for (; iter!=moar.end(); ++iter)
				to_pass.push_back(*iter);
			groups[num_groups]=to_pass;
			branch_lengths[num_groups]=blen;

			// Build new string in lieu of boost::regex_replace algorithm: significant speed increase
			string::const_iterator beg = newick.begin();	
			string::const_iterator end = newick.end();
			string nnewick(beg, matches[0].first);
			nnewick += replacement;
			string tmp(matches[0].second, end);
			nnewick += tmp;
			newick = nnewick;			


		}
	}
	void get_group_info (int group) 
	{
		cout << "Group #" << group << ":" << endl;
		cout << "branch length = " << branch_lengths[group] << endl;
		vector<int>& members = groups[group];
		cout << "group members --> ";
		vector<int>::const_iterator iter=members.begin();
		for (; iter != members.end(); ++iter) 
			cout << *iter << " ";
		cout << endl;
	}		
	void mark_mutated_group (int group)
	{
		mutated_group = group;
	}

	inline const string get_topology () 				{return topology;}
	inline void get_newick () 					{cout << newick << "\n is the specified newick\n" << endl;}
	inline const double& get_branch_lengths (int key) 		{return branch_lengths[key];}	
	inline const vector<int>& get_group_makeup (int key) 		{return groups[key];}
	inline const int get_group_number () 				{return num_groups;}
	inline const int get_sample_number ()				{return num_samples;}
	inline const int get_mutated_group()				{return mutated_group;}
	inline const bool get_snpstr_status()				{return snpstr_switch;}
	inline void trip_snpstr_switch()				{snpstr_switch = true;}

	~ Genealogy () {} 

	Genealogy(const string newick_output): newick(newick_output)   // the only Constructor so far
	{
		num_samples = 0;	
		num_groups = 0;
	
		// save topology by removing branch lengths
		boost::regex re, re_semi;
		re.assign("\\:\\d*.\\d\\d\\d"); // regexp for branch length format
		re_semi.assign("\\;"); 
		topology = boost::regex_replace (newick, re, ""); // remove branch lengths
		topology = boost::regex_replace (topology, re_semi, ""); // remove trailing semi-colon

		// create groups
		initialize_samples();
		group_samples();

		snpstr_switch = false; // will only be tripped true if a SNPSTR object is later used

	}

	
private:
	string newick; // original, newick-formatted tree output by MS
	string topology; // topology of the geneaology; empty string unless get_topology is called explicitly
	map<int, double> branch_lengths; // an associative array that holds the branch lengths of each external and group branch KEY: group# 
	map<int, vector<int> > groups; // an associative array that holds the members of each group KEY: group#
	map<unsigned int, signed int> samples; // KEY: sample#, VALUE: allele size
	unsigned int num_samples;	
	unsigned int num_groups;
	unsigned int mutated_group; // this variable is assigned to when the ISM model from Neutral_Variation.h invokes it through mark_mutated_group(int i);
	bool snpstr_switch;
};

#endif
