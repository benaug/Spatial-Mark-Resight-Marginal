# Spatial-Mark-Resight-Marginal
SMR samplers marginalizing out latent individual IDs. Poisson obervation model only. 
To speed up computation, I use the approach of Herliansyah et al. (2024, section 4.3) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

All samplers allow all latent ID observation types: marked without ID, unmarked, unknown marked status.

All samplers allow density covariates and a habitat mask, but density covariates can be excluded.

There are 6 types of models: 

1) known number of marked individuals

2) known number of marked individuals with interspersed marking and sighting

3) unknown number of marked individuals (for natural marks or premarked scenario when you no longer know number of marked inds in population)

4) generalized SMR (gSMR) with known number of marked individuals. This includes a marking process to account for different spatial distributions of marked and unmarked individuals (Whittington et al. 2018):

5) generalized SMR with known number of marked individuals and interpsersed marking and sighting (Whittington et al. 2018)
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954


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
theta.marked[3]=theta.unmarked[3]. Otherwise, bias will be introduced by not considering these samples, but we have no way to include them with this approach.

By using the data twice, there is a concern we will get biased estimates and underestimate posterior standard deviations.
Simulations in Margenau et al. (2022) and Whittington et al. (2025) show minimal bias and roughly nominal coverage when considering only marked with ID and unmarked sample types and no density variation.
I ran 1 simulation scenario with 25% of marked individual samples being unidentifiable and a density covariate and saw -7% bias in the density covariate beta, but effectively no bias in expected and realized N and other parameters. 
The CV for the density beta was 53% so this could just be due to limited information in the trap-level counts about spatial density variation compared to data with individual IDs where activity centers can be more precisely estimated.
It seems obvious that switching to One/Two Stage SMR will reduce statistical power with respect to the density covariate betas because activity center localization is very uncertain. But if you don't have the data
to model the marking process in gSMR, you will get biased abundance estimates (and I assume density covariate betas) with typical SMR approaches as well.


There are single and multisession versions of each model/sampler.
Need to add multisession gSMR with interspersed marking and sighting and One Stage. 

SMR models that allow observation models other than Poisson can be found here:

https://github.com/benaug/Spatial-Mark-Resight

SMR models with categorical partial IDs can be found here:

https://github.com/benaug/Spatial-Mark-Resight-IDCov

These are both more limited (e.g., no habitat mask, density covariates), but can be modified.