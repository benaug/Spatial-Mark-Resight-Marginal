# Spatial-Mark-Resight-Marginal
SMR samplers marginalizing out latent individual IDs. Poisson obervation model only. 
To speed up computation, I use the approach of Herliansyah et al. (2024, section 4.3) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

All samplers allow all latent ID observation types: marked without ID, unmarked, unknown marked status.

All samplers allow density covariates and a habitat mask, but density covariates can be excluded.

There are 4 types of models: 
1) known number of marked individuals
2) known number of marked individuals with interspersed marking and sighting
3) unknown number of marked individuals (for natural marks or premarked scenario when you no longer know number of marked inds in population)
4) generalized SMR (gSMR) with known number of marked individuals. This includes a marking process to account for different spatial distributions of marked and unmarked individuals (Whittington et al. 2018):
5) generalized SMR with known number of marked individuals and interpsersed marking and sighting (Whittington et al. 2018)
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954

There are single and multisession versions of each model/sampler.
Could add multisession gSMR with interspersed marking and sighting.

SMR models that allow observation models other than Poisson can be found here:

https://github.com/benaug/Spatial-Mark-Resight

SMR models with categorical partial IDs can be found here:

https://github.com/benaug/Spatial-Mark-Resight-IDCov

These are both more limited (e.g., no habitat mask, density covariates), but can be modified.