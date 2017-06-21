# Lplot
Local polynomial fitting (LOESS) of pairwise relatedness to pairwise geographic distance

In order to test if the average observed relatedness (r) predicted by LOESS at a given distance differs from the null model, 
row and column indices for the relatedness matrix are permuted n times; and at each permutation a LOESS model is re-fitted using 
the permuted relatedness and geographic distance matrix. 
95 % percentiles of the permutation-derived LOESS predictions are used to generate confidence envelopes around the null expectation of r = 0. 
Short vertical lines at the bottom of the figure are observed pairwise distances.

For details see: Bruno C, Macchiavelli R, Balzarini M (2008) Non-parametric smoothing of multivariate genetic distances in the analysis of spatial population structure at fine scale. Theor Appl Genet 117:435–447

Function created by Nathaniel S. Pope: https://github.com/nspope
