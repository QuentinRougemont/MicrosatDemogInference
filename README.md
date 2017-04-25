# DemographicInference
Contains  Script To performed Demographic inference using microsatelittes data with coalescent theory, ABC model selection and parameter estimation. I also explore the use of Random Forest for model selection.
All R scripts for simulating data and convertion into R were modified from previous developpement by Illera et al. (2014).
Original script made by M. Nasvascues can be found here:
https://www.researchgate.net/publication/273445598_ABCSylvia

Ref: 
Illera, J. C., A. M. Palmero, P. Laiolo, F. Rodríguez, Á. C. Moreno, and M. Navascués. 2014. Genetic, Morphological and Acoustic Evidence Reveals Lack of Diversification in the Colonization Process in an Island Bird. Evolution Volume 68, Issue 8, pages 2259-2274.DOI: 10.1111/evo.12429
available here: http://onlinelibrary.wiley.com/doi/10.1111/evo.12429/full

Rougemont, Q., C. Roux, S. Neuenschwander, J. Goudet, S. Launey, G. Evanno. 2016. Reconstructing the demographic history of divergence between European river and brook lampreys using approximate Bayesian computations. PeerJ 4:e1910 https://doi.org/10.7717/peerj.1910

## Dependencies

```ms``` coalescent simulator from hudson available at [download ms](https://uchicago.app.box.com/s/l3e5uf13tikfjm7e1il1eujitlsjdx13)

## R dependencies :
```pegas``` packages

## Recent update : 
 implemented a 3 pops colaescent derivation of our two population models for Isolation w. Migration (IM) or Strict Isolation (SI)

##Original implementation
Allows the comparaison of 4 models of population divergence namely 
* SI (Strict Isolation) , 
* Isolation with Migration (IM), 
* Ancient Migration (AM) and 
* Secondary Contact (SC)
