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


#disclaimer, number 6 below needs more thought. My data simulator didn't lead to different activity center distributions for marked and unmarked individuals. Will see if it still works sharing marked guy z's and s's across marked and combined data and report back.

6) One Stage SMR where the marked individual data is used twice (Whittington et al. 2025). This approach provides a means to account for different spatial distributions of marked and unmarked individuals. Simulations in Whittington et al. show minimal bias and roughly nominal coverage as does the single simulation scenario I ran (only issue was slightly low coverage for lam0, 0.90 instead of 0.95).
My version here differs from Whittington et al. in that 1) I allow for marked but no ID detections to be included (ignoring these introduces bias), 2) I estimate single activity centers for marked individuals instead of one each for the data with and without individual IDs, 
and 3) I regard the marked individuals to be included in the population with certainty (Whittington et al. estimates z indicator for marked individuals).
Note, unknown marked status samples that came from marked individuals cannot be accounted for and if you have these, to my knowledge, there is no way to correct the bias introduced using One or Two Stage approaches.
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.70246

There are single and multisession versions of each model/sampler.
Need to add multisession gSMR with interspersed marking and sighting and One Stage. 

SMR models that allow observation models other than Poisson can be found here:

https://github.com/benaug/Spatial-Mark-Resight

SMR models with categorical partial IDs can be found here:

https://github.com/benaug/Spatial-Mark-Resight-IDCov

These are both more limited (e.g., no habitat mask, density covariates), but can be modified.