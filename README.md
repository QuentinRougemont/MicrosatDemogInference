# DemographicInference
Contains  Script To perform Demographic inference using microsatelittes data: 
based on :
* coalescent simulation usgin ```ms``` 
* ABC model selection and parameter estimation. 
* RandomForest model selection 

R scripts for simulating data and convertion into R were modified from previous developpement by Illera et al. (2014). Original script made by M. Nasvascues can be found here:
https://www.researchgate.net/publication/273445598_ABCSylvia

## new features includes:
The comparaison of 4 models of population divergence namely 
* SI (Strict Isolation) , 
* Isolation with Migration (IM), 
* Ancient Migration (AM) and 
* Secondary Contact (SC)

Computation of other summary statistics:
* Allelic Richness (Ar)
* Expected Heterozygosity
with the previous stats being include in function to reduce the complexity of the scripts

* inclusion of R script for:
1. model comparaison
2. parameter estimation
3. computation of robutness
4. goodness of fit
5. random forest (also I should update this script one day)

* inclusion in bash scripts for parallelisation and runs easily using any microstallite data.
microstallite data that have to be stored in the `01-data` folder together with a file containing the `repeat_motif` for each microsatellites markers.

I also extended the pipeline to perform :  
`3 pops` colaescent derivation of our two population models for:
* Isolation w. Migration (IM) 
* Strict Isolation (SI)

 `4 pops` colaescent derivation for:
* Isolation w. Migration (IM) (with several possible configurations for directions and symetries of introgression) 
* Strict Isolation (SI)

## References:

* Illera, J. C., A. M. Palmero, P. Laiolo, F. Rodríguez, Á. C. Moreno, and M. Navascués. 2014. Genetic, Morphological and Acoustic Evidence Reveals Lack of Diversification in the Colonization Process in an Island Bird. Evolution Volume 68, Issue 8, pages 2259-2274.DOI: 10.1111/evo.12429
available here: http://onlinelibrary.wiley.com/doi/10.1111/evo.12429/full

* Rougemont, Q., C. Roux, S. Neuenschwander, J. Goudet, S. Launey, G. Evanno. 2016. Reconstructing the demographic history of divergence between European river and brook lampreys using approximate Bayesian computations. PeerJ 4:e1910 https://doi.org/10.7717/peerj.1910

## Dependencies

`ms`coalescent simulator from Hudson 2002 available  [here](https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13)

## R dependencies :
* `pegas` package more info [here](https://cran.r-project.org/web/packages/pegas/index.html)
* `abc` package  more info [here](https://cran.r-project.org/web/packages/abc/index.html)
* `RandomForest SRC` package more info [here](https://cran.r-project.org/web/packages/randomForestSRC/index.html)

## Major steps to run the pipeline

### Prearing your data

You'll need to provide a separate input file for each population. Each file contains one markers in raw, with one row for each gene copies (2 rows per individuals) and loci in columns  
A file containing the length of the repeat_motif for each microsatellite markers is also needed.  
One row per markers.   
These input file are stored in the `01-data` folder

then go to `00-scripts/models/` and edit the script `model.\*.sh` to provide the name of the input file   
the script `model.1` sh is for SI  
`model.2.sh` for  IM   
`model.3.sh` for SC  
`model.4.sh` for AM  

the corresponding Rscripts are found in the `00-scripts/rscript/` with name `Simul_*_parallel.R `

## Chooose the prior

* go in ``00-scripts/rscript/` and edit the Simul\*.R scripts to set prior according to what you think will fit the data.  
I recommand to use large and uninformative priors first.  
Wou'll have to choose a fixed thetaRef, according to thetaRef=4\*Nref\*µ meaining that you need to have an idea of the mutation rate (ideally)
then you have to choose priors for :  
* effective population size (N1, N2, Nancestral, etc)
* migration rates (M1, M2, etc)
* Split time (T=Tsplit/4Nref)

beware that in `ms` from Hudson, all is scaled by 4Nref so that you need to thinks mostly in terms of ratios for the effective population size and Split times.
I strongly advise reading of `msdoc.pdf` included in the `00-scripts/msdir ` 

## Runing the pipeline

once the prior are set coalescents simulations can be run and needs ideally to be run in parallel if you want to save time.
to do so I provided a series of 4 scripts (one for each model) in `00-scripts/models/model_job_array\*.sh` These scripts are rather 'cluster-specific' and depending on your clustering machine, you'll have to change that but this provides a clue.

