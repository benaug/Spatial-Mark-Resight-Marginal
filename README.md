# Spatial-Mark-Resight-Marginal
SMR samplers marginalizing out latent individual IDs. Poisson obervation model only. 
To speed up computation, I use the approach of Herliansyah et al. (2024) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

These models use count prior data augmentation: https://github.com/benaug/SCR-Count-Prior-Data-Augmentation

All samplers allow all latent ID observation types: marked without ID, unmarked, unknown marked status.

All samplers allow density covariates and a habitat mask, but density covariates can be excluded.

There are 3 types of models: premarked (known number of marked individuals), premarked with interspersed marking and sighting, and natural marks (unknown number of marked individuals).

For generalized SMR, you can add the marking process to the current code.

There are single and multisession versions of each model/sampler.

SMR models that allow observation models other than Poisson can be found here:

https://github.com/benaug/Spatial-Mark-Resight

SMR models with categorical partial IDs can be found here:

https://github.com/benaug/Spatial-Mark-Resight-IDCov

These are both more limited, but can be modified.