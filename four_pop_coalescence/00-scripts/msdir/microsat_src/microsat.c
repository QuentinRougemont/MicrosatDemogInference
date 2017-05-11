/*
  name: microsat.c
  
  author:        murray cox <mpcox@u.arizona.edu>
  date:          august 2007
  
  original code: richard hudson (university of chicago)
  
  this code takes ms output and converts it to microsatellite length variation
  data.  the output has on each line the set of lengths of the nsam individuals
  (relative to the ancestral length).  mpc:- the code now simulates fully linked
  STRs.  the proportion of the total coalescent theta assignable to each STR must
  be specified.  the ancestral character state can be defined (defaults to zero).
  only a single step mutation model (SMM) is implemented.
  
  usage:         microsat [-a anc_state] [-l n_linked_sites theta_prop1 ... theta_prop_n]
  example:       ms 10 5 -t 4.0 | microsat -a 30 -l 2 0.4 0.6 > msat.dat
  compilation:   gcc -o microsat microsat.c rand1.c -lm
  
  Copyright (C) 2007 Murray P. Cox

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// standard libraries
#include <stdio.h>
#include <stdlib.h>

// define statements
#define MAX_STRS 1000
#define MAX_LINE_LENGTH 1000
#define TOLERANCE 0.00000000000001

// globals and function prototypes
int MAX_SITES = 1000;
double ran1();
char **cmatrix(int, int);
void biggerlist(int nsam, unsigned MAX_SITES, char ** list);

int main(int argc, char **argv) {
	
	// local variables
	int nsam, i, ms_datasets, segsites, count, step, ind, flag_c;
	int ancstate = 0, linked_n = 1, im_flag = 0;
	double *posit, prob, thetas[MAX_STRS];
	double theta_sum = 0;
	char **list, line[MAX_LINE_LENGTH+1], dummy[20], astr[100];
	FILE *pfin;
	
	/* read in first two lines (parameters and random number seed) */
	pfin = stdin;
	fgets(line, MAX_LINE_LENGTH, pfin);
	sscanf(line, " %s  %d %d", dummy,  &nsam, &ms_datasets);
	fgets(line, MAX_LINE_LENGTH, pfin);
	
	/* read command line arguments, if any */
	if( argc > 1 ) {
		for(flag_c = 1; flag_c < argc; ++flag_c ) {
			if(argv[flag_c][0] == '-') {

				switch( argv[flag_c][1] ) {
					
					// set ancestral state
					case 'a':
						++flag_c;
						if(argv[flag_c] == NULL){
							fprintf(stderr, "error: no values following -a flag\n");
							exit(1);									
						}
						ancstate = atoi(argv[flag_c]);
						break;
					
					case 'i':
						im_flag = 1;
						break;

					// number of linked STRs and theta proportions
					case 'l':
						++flag_c;
						if(argv[flag_c] == NULL){
							fprintf(stderr, "error: no values following -l flag\n");
							exit(1);									
						}
						// collect number of STRs to simulate
						linked_n = atoi (argv[flag_c]);
						
						// collect proportion of total theta for each STR	
						if(linked_n != 1){
							int j;
							for(j = 0; j < linked_n; ++j){
								++flag_c;
								
								if(argv[flag_c] == NULL){
									fprintf(stderr, "error: only %i thetas for %i STRs\n", j, linked_n);
									exit(1);									
								}
								
								thetas[j] = atof (argv[flag_c]);
								theta_sum += thetas[j];
								//printf("theta_sum is %f\n", theta_sum);
							}
							if(1 - theta_sum > TOLERANCE || theta_sum - 1 > TOLERANCE){
								fprintf(stderr, "error: sum of thetas (= %f) exceeds 1\n", theta_sum);
								exit(1);
							}
						}
						break;
						
					default:
						fprintf(stderr, "error: option %s unknown\n", argv[flag_c]);
						exit(1);
				}
				
			}else ++flag_c;
		}
	}
	
	// set default theta[0]
	if(linked_n == 1) {
		thetas[0] = 1.0;
	}
	
	// allocate data arrays 
	list = cmatrix(nsam, MAX_SITES+1);
	posit = (double *)malloc( MAX_SITES*sizeof(double) );
	if (posit == NULL) {
		fprintf(stderr, "error: malloc failure 1 for posit\n");
		exit(1); 
	}
	
	// steps through all ms datasets
	count = 0;
	while( ms_datasets - count++ ) {

		// read in a sample
		do {
			if( fgets( line, MAX_LINE_LENGTH, pfin) == NULL ){
				perror("error: line retrieval fault\n");
				exit(1);
			}
		}while ( line[0] != '/' );
	
		// read in number of segsites
		fscanf(pfin,"  segsites: %d", &segsites);
	
		// reallocate memory if segsites greater than MAX_SITES
		if( segsites >= MAX_SITES){
			MAX_SITES = segsites + 10 ;
			posit = (double *)realloc( posit, MAX_SITES*sizeof( double) ) ;
			if (posit == NULL) {
				fprintf(stderr, "error: realloc failure 1 for posit\n");
				exit(1);
			}
			void biggerlist(int nsam, unsigned MAX_SITES, char ** list);
		}
		
		// if dataset contains segregating sites
		if( segsites > 0) {
			
			// positioning read statement
			fscanf(pfin," %s", astr);
			
			// not sure what this routine does
			if( astr[1] == 'r' ){
				fscanf(pfin," %lf", &prob );
				fscanf(pfin," %*s");
			} 
		
		// posit array captures site positions, list array captures segsite strings
		for( i = 0; i < segsites; i++) fscanf(pfin, " %lf", posit+i);
		for( i = 0; i < nsam; i++)     fscanf(pfin, " %s", list[i]);
		}
		
		// analyse sample ( do stuff with segsites and list)
		// initialize nrepeats array with ancestral state
		int nrep[nsam][linked_n];
		for(ind = 0; ind < nsam; ind++) {
			int j;
			for(j = 0; j < linked_n; j++){
				nrep[ind][j] = ancstate;
			}
		}
		
		// randomly choose and assign +1/-1 single step mutation
		double random, value;
		int this_str;
		for(i = 0; i < segsites; i++){
			
			//printf("segsite is %d\t", i);
			
			if(ran1() < 0.5) step = -1;
			else step = 1;
			//printf("step is %d\t", step);
			
			random = ran1();
			//printf("ran is %f\t", random);
			
			this_str = 0;
			value = 0.0;
			int j = 0;
			for(j = 0; j < linked_n; ++j){
				value += thetas[j];
				if(random > value){
					continue;	
				}else{
					this_str = j;
					break;
				}
			}
			//printf("str is %d\n", this_str);
			
			for(ind = 0; ind < nsam; ind++){
				
				if( list[ind][i] == '1' ) {
					nrep[ind][this_str] += step;
					//printf("%d %d %d\n", i, ind, nrep[ind][this_str]);
				}
			}
		}

		// print datasets to stdout
		int a = 0, b = 0;
		if(im_flag == 0) {		
			for(a = 0; a < linked_n; ++a){
				for(b = 0; b < nsam; ++b){
					if(b == nsam - 1 && a == linked_n -1){
						printf("%d\n", nrep[b][a]);
					}else{
						printf("%d\t", nrep[b][a]);
					}
				}
			}
		}else{
			for(a = 0; a < nsam; ++a){
				for(b = 0; b < linked_n - 1; ++b){
					printf("%d\t", nrep[a][b]);
				}
				printf("%d\n", nrep[a][linked_n - 1]);
			}
			printf("//\n");
		}
		
	} // end of while dataset loop
 
	return(0);
} // end of main function


/* allocates space for gametes (character strings) */
char ** cmatrix(int nsam, int len) {
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) ){
		perror("error: malloc failure 1 in cmatrix\n");
		exit(1);
	}
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) ))){
			perror("error: malloc failure 2 in cmatrix\n");
			exit(1);
		}
	}
	return(m);
}

/* creates more memory for list */
void biggerlist(int nsam, unsigned nmax, char ** list ) {

	int i;
	MAX_SITES = nmax;
	for(i = 0; i < nsam; i++){
		list[i] = (char *)realloc( list[i], MAX_SITES*sizeof(char) ) ;
		if( list[i] == NULL ){
			perror("error: realloc failure 1 in biggerlist\n");
			exit(1);
		}
	}
}
