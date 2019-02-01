############################################################################################################################################################################
############################################# Local polynomial fitting (LOESS) of pairwise relatedness to pairwise geographic distance #####################################
############################################################################################################################################################################

# In order to test if the average observed relatedness (r) predicted by LOESS at a given distance differs from the null model, 
# row and column indices for the relatedness matrix are permuted n times; and at each permutation a LOESS model is re-fitted using 
# the permuted relatedness and geographic distance matrix. 
# 95 % percentiles of the permutation-derived LOESS predictions are used to generate confidence envelopes around the null expectation of r = 0. 
# Short vertical lines at the bottom of the figure are observed pairwise distances.
# For details see: Bruno C, Macchiavelli R, Balzarini M (2008) Non-parametric smoothing of multivariate genetic distances in the analysis of spatial population structure at fine scale. Theor Appl Genet 117:435–447
# Function created by Nathaniel S. Pope: https://github.com/nspope

##### Loess function
mDistoLoess <- function(mat1,mat2,nperm=999){
ltd <- mat1[lower.tri(mat1)]
ltk <- mat2[lower.tri(mat2)]
disto <- matrix(nrow = nperm, ncol = length(ltk))
for(n in 1:nperm){
	perm <- sample(1:nrow(mat2))
	mat2_p <- mat2[perm,perm]
	lt <- mat2_p[lower.tri(mat2_p)]
	disto[n,] <- predict(loess(lt~ltd))
	if(n %% 100 == 0) print(n)
	}
obs_fit <- loess(ltk~ltd)
obs <- predict(obs_fit)
CI <- apply(disto, 2, quantile, probs = c(0.025,0.975))
attr(disto, "obs") <- obs
attr(disto, "CI") <- CI
return(list(disto=disto))
}

###### Loess plot function
Lplot <- function(distance, relatedness, permutations){
  library(ggplot2)
  ot <- mDistoLoess(distance, relatedness, permutations) ## Takes a long time!! Leave overnight
  dd <- data.frame(cov = distance[lower.tri(distance)], obs = attr(ot$disto, "obs"), lci = attr(ot$disto, "CI")[1,], uci = attr(ot$disto, "CI")[2,])
  P <- ggplot(dd, aes(x=cov,y=obs)) + geom_line(size=1.2) + geom_ribbon(aes(ymax = uci, ymin = lci), alpha = 0.3) + theme_minimal() + 
    geom_hline(yintercept = mean(ot$disto), lty = 3) + ylab("Mean relatedness") + xlab("Distance (Km)") + geom_rug(sides = "b", alpha = 0.02) +
    theme(axis.title.y = element_text(size=20, color = "black", face = "bold"),
          axis.title.x = element_text(size=20, color = "black", face = "bold"))
return(P)
  }

################ Usage
D <- matrix(ncol=100, nrow=100, rnorm(100)) # D is a geographic distance matrix (a full matrix NOT a dist object)
R <- matrix(ncol=100, nrow=100, rnorm(100)) # R is a relatedness matrix (a full matrix NOT a dist object)
P1 <- Lplot(distance=D, relatedness=R, permutations=999) ##
P1
