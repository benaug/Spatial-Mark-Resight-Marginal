# Spatial-Mark-Resight-Marginal
SMR samplers marginalizing out latent individual IDs. Poisson observation model only. 
To speed up computation, I use the approach of Herliansyah et al. (2024, section 4.3) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

All samplers allow all latent ID observation types: marked without ID, unmarked, unknown marked status (except One Stage SMR).

All samplers allow density covariates and a habitat mask, but density covariates can be excluded.

There are 6 types of models: 

1) known number of marked individuals (Chandler and Royle 2013, Sollmann et al. 2013)

https://www.jstor.org/stable/23566419
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-1256.1

2) known number of marked individuals with interspersed marking and sighting (Whittington et al. 2018 with no marking process)

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954

3) unknown number of marked individuals (for natural marks or premarked scenario when you no longer know number of marked inds in population)
(Rich et al. 2014)

https://academic.oup.com/jmammal/article/95/2/382/866592

4) generalized SMR (gSMR) with known number of marked individuals. This includes a marking process to account for different spatial distributions of marked and unmarked individuals (Whittington et al. 2018):

5) generalized SMR with known number of marked individuals and interspersed marking and sighting (Whittington et al. 2018)

6) One Stage SMR where the marked individual data is used twice (Whittington et al. 2025): 

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.70246

This approach provides a means to account for different spatial distributions of marked and unmarked individuals without having to model the marking process (gSMR).
It is a modification of the Two Stage approach of Margenau et al. (2022) that does both stages at once:

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/eap.2553

My version here differs from Whittington et al. (2025) in that I allow for marked but no ID detections to be included. Margenau et al. (2022) say 
"In the instance a marked individual cannot be reliably recognized in camera photographs, perhaps due to body positioning, vegetation obstruction, or blurriness, the data record should be discarded from stage one of the model but included as a record in stage two."
, but this suggestion will introduce bias, whether they are talking about marked with no ID or unknown marked status individuals (I can't tell 100%). 
One and Two Stage approaches require a solution for unknown marked status samples if they occur. 
In the specific case where marked and unmarked individuals' detections are equally likely to be unknown mark status events, you can use
One and Two Stage approaches without introducing bias by removing these samples from *both stages* of the analysis. More concretely, I mean when
theta.marked[3]=theta.unmarked[3] (these are explicitly defined at top of single session test scripts). Otherwise, bias will be introduced by not considering these samples, but we have no way to include them with this approach.

By using the data twice, there is a concern we will get biased estimates and underestimate posterior standard deviations.
Simulations in Margenau et al. (2022) and Whittington et al. (2025) show minimal bias and roughly nominal coverage when considering only marked with ID and unmarked sample types and no density variation.
I ran 1 scenario each of single and multisession. There was -6.5% bias for the density covariate beta in the single session (but nominal coverage) and minimal bias in the multisession
scenario. I probably needed to use more than 100 data sets for the single session scenario given how high the CV for the density covariate beta
was. In both, the coverage of lam0 was around 0.90, so perhaps a bit low from using the data twice. Overall, looks good.

There are single and multisession versions of each model/sampler.
 
SMR models that allow observation models other than Poisson can be found here:

https://github.com/benaug/Spatial-Mark-Resight

SMR models with categorical partial IDs can be found here:

https://github.com/benaug/Spatial-Mark-Resight-IDCov

These are both more limited (e.g., no habitat mask, density covariates), but can be modified.

Analogous repositories in the "Marginal Unmarked Trilogy" can be found here:

Unmarked SCR: https://github.com/benaug/Unmarked-SCR-Marginal

SCR with random thinning: https://github.com/benaug/Random-Thin-Marginal