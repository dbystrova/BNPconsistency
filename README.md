R code to replicate figures of the article: Bayesian mixture models (in)consistency for the number of clusters. 
https://arxiv.org/abs/2210.14201

Joint work by Louise Alamichel, Daria Bystrova, Julyan Arbel, Guillaume Kon Kam King.

## Abstract

Bayesian nonparametric mixture models are common for modeling complex data. 
While these models are well-suited for density estimation, their application for clustering has some limitations. Recent results proved posterior inconsistency of the number of clusters when the true number of clusters is finite for the Dirichlet process and Pitman--Yor process mixture models. 
We extend these results to additional Bayesian nonparametric priors such as Gibbs-type processes and finite-dimensional representations thereof. The latter include the Dirichlet multinomial process, the recently proposed Pitman--Yor, and normalized generalized gamma multinomial processes. 
We show that mixture models based on these processes are also inconsistent in the number of clusters and discuss possible solutions. Notably, we show that a post-processing algorithm introduced for the Dirichlet process can be extended to more general models and provides a consistent method to estimate the number of components.

To reproduce the results of the paper run `make figures` in the BNPconsistency folder.
