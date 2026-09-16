# Spatial-Mark-Resight-Marginal

Closed-population spatial mark-resight (SMR) models using N-prior data augmentation and marginalization over latent 
individual identities. Poisson sighting observation model only.

## Overview

This repository contains single and multisession closed-population spatial mark-resight (SMR) models with Poisson
sighting observation models and marginalization over latent individual identities.

For the Poisson sighting model, individual identity can be marginalized analytically for unidentified observations using
Poisson thinning and superposition. This avoids sampling latent individual identities during MCMC and generally improves
mixing relative to conditional formulations that explicitly update those identities.

These models use [N-prior data augmentation](https://github.com/benaug/SCR-N-Prior-Data-Augmentation). To improve
computational efficiency for latent identity observation types, the approach described in Section 4.3 of
[Herliansyah et al. (2024)](https://link.springer.com/article/10.1007/s13253-023-00598-3) is used in the
custom z and activity center samplers.

All samplers allow the following observation types given detection:

- marked with ID,
- marked no ID,
- unmarked, and
- unknown mark status,

except for the one-stage SMR models, which have additional restrictions described below.

The repository includes several SMR formulations. Standard models assume the number of marked individuals is known,
whereas the `Natural` models estimate marked and unmarked abundance separately when the current number of marked
individuals is unknown. Generalized SMR models explicitly model the marking process, and `Interspersed` versions allow
mark status to vary among sighting occasions within a closed session. Additional models accommodate multiple marked
states, such as working versus failed GPS collars, or implement the one-stage SMR approach.

All samplers allow density covariates and a habitat mask, although density covariates can be excluded. The multisession
models are collections of closed-population sessions rather than an open-population model: abundance, density, and
activity centers can differ among sessions, but there is no demographic linkage among sessions. Integrated telemetry
locations can also be included to inform activity centers and the spatial scale parameter.

## Related repositories

SMR models that allow observation models other than Poisson are available in the
[Spatial-Mark-Resight-Conditional](https://github.com/benaug/Spatial-Mark-Resight-Conditional) repository. SMR
models incorporating categorical partial identities are available in the
[Spatial-Mark-Resight-IDCov](https://github.com/benaug/Spatial-Mark-Resight-IDCov) repository.

Analogous repositories using the marginal observation model for other combinations of marked and unmarked data types
include:

- [Unmarked SCR](https://github.com/benaug/Unmarked-SCR-Marginal)
- [SCR with random thinning](https://github.com/benaug/Random-Thin-Marginal)
- [SCR with integrated occupancy data](https://github.com/benaug/SCR_Dcov_IntegratedOccupancy)

Open-population versions of the marginalized generalized SMR models are available here:

- [Spatial-Mark-Resight-Open-Marginal](https://github.com/benaug/Spatial-Mark-Resight-Open-Marginal)


## Spatial mark-resight observation model

The observation model is a Poisson sighting process and, for the generalized SMR models, a separate marking process.
Individual identity is known for marked-ID detections and latent for the unidentified sighting types. The equations
below first describe models in which the number of marked individuals is known. The same Poisson thinning and individual ID
marginalization are used in the `Natural` models, but marked abundance is estimated rather than known. For clarity, 
the equations below describe the standard single-session model first. Multisession and interspersed versions
are subsequently described.

### Sighting process

For individual $i$ and sighting detector $j$, let

$$
\lambda_{i,j}=z_i\lambda_0\exp\left(-\frac{\|\mathbf{s}_i-\mathbf{x}_j^{S}\|^2}{2\sigma^2}\right),
$$

where $z_i$ is the population inclusion indicator, $\mathbf{s}_i$ is the activity center of individual $i$,
$`\mathbf{x}_j^{S}`$ is the location of sighting detector $j$, $\lambda_0$ is the baseline sighting rate,
and $\sigma$ is the spatial scale parameter. The latent true number of sightings is

$$
Y_{i,j}^{\mathrm{true}}\sim\mathrm{Poisson}\left(K_j^{S}\lambda_{i,j}\right),
$$

where $K_j^{S}$ is sighting effort. Let $m_i=1$ if individual $i$ is marked during the sighting period and $m_i=0$
otherwise. For marked individuals, define

$$
\boldsymbol{\theta}^{M}=\left(\theta_{\mathrm{ID}}^{M},\theta_{\mathrm{noID}}^{M},\theta_{\mathrm{unk}}^{M}\right),
$$

where the three elements are the probabilities that a sighting is recorded as marked with ID, marked no ID,
or unknown mark status. For unmarked individuals, define

$$
\boldsymbol{\theta}^{U}=\left(0,\theta_{\mathrm{um}}^{U},\theta_{\mathrm{unk}}^{U}\right),
$$

where

$$
\theta_{\mathrm{unk}}^{U}=1-\theta_{\mathrm{um}}^{U}.
$$

The first element is zero because an unmarked individual cannot produce an individually identifiable observation. 
Conditional on the latent true count, sightings are thinned into observed sample types using a multinomial distribution. If individual $i$ is marked,

$$
\left(Y_{i,j}^{\mathrm{ID}},Y_{i,j}^{\mathrm{noID}},Y_{i,j}^{\mathrm{unk},M}\right)\mid
Y_{i,j}^{\mathrm{true}},m_i=1\sim\mathrm{Multinomial}
\left(Y_{i,j}^{\mathrm{true}},
\left(\theta_{\mathrm{ID}}^{M},\theta_{\mathrm{noID}}^{M},\theta_{\mathrm{unk}}^{M}\right)
\right).
$$

If individual $i$ is unmarked,

$$
\left(Y_{i,j}^{\mathrm{ID}},Y_{i,j}^{\mathrm{um}},Y_{i,j}^{\mathrm{unk},U}\right)\mid
Y_{i,j}^{\mathrm{true}},m_i=0\sim\mathrm{Multinomial}
\left(Y_{i,j}^{\mathrm{true}},
\left(0,\theta_{\mathrm{um}}^{U},\theta_{\mathrm{unk}}^{U}\right)
\right).
$$

This is the conditional data generating model, which we use for estimation in the conditional SMR repository. In
the models in this repository, both the latent true sighting counts and the latent individual identities of
unidentified samples are marginalized analytically.

### Marking process in generalized SMR

The generalized SMR models include the process by which individuals are captured and marked. This accounts for
different spatial distributions of marked and unmarked individuals that arise when marking effort is spatially
nonrandom.

For individual $i$ and marking detector $j$, detection probability is modeled using a halfnormal detection function,

$$
p_{i,j}=p_0\exp\left(-\frac{\|\mathbf{s}_i-\mathbf{x}_j^{M}\|^2}{2\sigma^2}\right),
$$

where $`\mathbf{x}_j^{M}`$ is the location of marking detector $j$ and $p_0$ is baseline marking detection probability.

The marking observations follow

$$
Y_{i,j}^{M}\sim\mathrm{Binomial}\left(K_j^{M},p_{i,j}\right),
$$

where $K_j^{M}$ is marking effort.

The supplied model files use the same $\sigma$ for the marking and sighting processes, although this can be modified.

The single-session `Mb` version allows baseline detection during the marking process to differ between first capture
and subsequent captures, using separate parameters `p0.p` and `p0.c`.

### Multisession models

For multisession models, the same closed-population observation model is repeated for each session $g$. Quantities 
therefore become session specific, for example,

$$
\lambda_{i,g,j},\qquad\mathbf{s}_{i,g},\qquad K_{g,j}^{S},\qquad N_g.
$$

The supplied multisession examples estimate density parameters separately by session while sharing $\lambda_0$,
$\sigma$, and the sample type probabilities across sessions. These parameters can instead be estimated separately
by session or modeled hierarchically.

### Interspersed marking and sighting

In the standard models, mark state is constant over the sighting period. Sightings can therefore be summed across
occasions before fitting the model. The `Interspersed` versions retain the occasion dimension because an individual's mark state can change
during the sampling period, either because they are marked after the first sighting occasion, or because marks are lost. Let

$$
m_{i,k}=\begin{cases}
1, & \text{if individual } i \text{ is marked at sighting occasion } k,\\
0, & \text{otherwise}.\end{cases}
$$

An individual can therefore contribute unmarked sightings before capture and marked sightings after capture within 
the same closed-population session. Further, an individual may lose a mark and then only contribute unmarked sightings
again, be remarked, etc. Sighting effort and observed counts are retained by occasion, $K_{j,k}^{S}$, and the Poisson
thinning and individual ID marginalization described below are applied separately for each $k$ using the
appropriate mark state $m_{i,k}$. For multisession interspersed models, both indices are retained, giving quantities
such as $m_{i,g,k}$ and $K_{g,j,k}^{S}$.

### Unknown-mark-status samples

To my knowledge, no other SMR software considers the possibility of unknown-mark-status samples aside from my SMR
GitHub repositories. These detection types are equivalent to latent identity detections in unmarked SCR (Chandler 
and Royle, 2013) and carry
the least information of all observation types. Marked-no-ID and unmarked samples each carry more information
than unknown-mark-status samples because a subset of individuals can be excluded as the source. The only partial
identity information these samples contain is the spatial location of detection. 

Unknown-mark-status samples can be safely discarded, often with very little effect on parameter estimate precision, if marked and unmarked
individuals are equally likely to have their detections classified as unknown mark status, or, more specifically, when
$\theta_{\mathrm{unk}}^{M}=\theta_{\mathrm{unk}}^{U}$. If this condition does not hold, discarding these samples can
introduce bias, though the magnitude may be small in many scenarios. Relevant factors are how different these probabilities
are and what proportion of total observations are of this type.

A secondary benefit of explicitly including these samples in the model is to better differentiate them from unmarked
detections and improve how SMR detection types are classified. While anecdotal, I have noticed some researchers
classify any observation where a mark cannot be seen as unmarked, even when the individual could have been marked and 
the mark was simply not seen. For example, with GPS-collared deer, a camera image may not show the neck area where the collar is located. Classification
can be even more difficult with natural marks because a mark may be localized to one part of the animal, such as a scar,
and photographs may not show the relevant body region for every marked individual. These misclassifications will bias parameter estimates and
I suspect they are common. Therefore, I want to stress that **"mark not seen" does not necessarily mean "unmarked"**.

## Marginalization over sample type and individual identity

This section describes the standard single-session model. The same logic applies separately within each session of
the multisession models and within each sighting occasion of the interspersed models. Because the latent true count
is Poisson and the observed sample types are generated by multinomial thinning, the 
individual-level counts in each sample type are independent Poisson random variables after marginalizing over
$Y_{i,j}^{\mathrm{true}}$.

For a marked individual,

$$
Y_{i,j}^{\mathrm{ID}}\sim\mathrm{Poisson}\left(K_j^{S}\lambda_{i,j}\theta_{\mathrm{ID}}^{M}\right),
$$

$$
Y_{i,j}^{\mathrm{noID}}\sim\mathrm{Poisson}\left(K_j^{S}\lambda_{i,j}\theta_{\mathrm{noID}}^{M}\right),
$$

and

$$
Y_{i,j}^{\mathrm{unk},M}\sim\mathrm{Poisson}\left(K_j^{S}\lambda_{i,j}\theta_{\mathrm{unk}}^{M}\right).
$$

For an unmarked individual,

$$
Y_{i,j}^{\mathrm{um}}\sim\mathrm{Poisson}\left(K_j^{S}\lambda_{i,j}\theta_{\mathrm{um}}^{U}\right),
$$

and

$$
Y_{i,j}^{\mathrm{unk},U}\sim\mathrm{Poisson}\left(K_j^{S}\lambda_{i,j}\theta_{\mathrm{unk}}^{U}\right).
$$

Marked-ID sightings retain individual identity and remain in the likelihood at the individual level.
For the unidentified observation types, individual identity is not observed. Define the total underlying sighting
intensity from marked individuals as

$$
\Lambda_j^{M}=\sum_{i=1}^{M}m_i\lambda_{i,j},
$$

and from unmarked individuals as

$$
\Lambda_j^{U}=\sum_{i=1}^{M}(1-m_i)\lambda_{i,j}.
$$

The observed marked-no-ID count is

$$
Y_j^{\mathrm{noID}}=\sum_{i=1}^{M}Y_{i,j}^{\mathrm{noID}}.
$$

By the superposition property of independent Poisson random variables,

$$
Y_j^{\mathrm{noID}}\sim\mathrm{Poisson}\left(K_j^{S}\theta_{\mathrm{noID}}^{M}\Lambda_j^{M}\right).
$$

Similarly,

$$
Y_j^{\mathrm{um}}=\sum_{i=1}^{M}Y_{i,j}^{\mathrm{um}},
$$

with

$$
Y_j^{\mathrm{um}}\sim\mathrm{Poisson}\left(K_j^{S}\theta_{\mathrm{um}}^{U}\Lambda_j^{U}\right).
$$

Sightings with unknown mark status can arise from either marked or unmarked individuals,

$$
Y_j^{\mathrm{unk}}=\sum_{i=1}^{M}
\left(Y_{i,j}^{\mathrm{unk},M}+Y_{i,j}^{\mathrm{unk},U}\right),
$$

giving

$$
Y_j^{\mathrm{unk}}\sim\mathrm{Poisson}
\left(K_j^{S}
\left(\theta_{\mathrm{unk}}^{M}\Lambda_j^{M}+\theta_{\mathrm{unk}}^{U}\Lambda_j^{U}\right)
\right).
$$

The marginalized likelihood therefore follows from two properties of the Poisson distribution:

1. **Poisson thinning:** multinomial classification of a Poisson count produces independent Poisson counts for the 
resulting sample types.
2. **Poisson superposition:** summing independent Poisson counts across latent individual identities produces another
Poisson count with rate equal to the sum of the individual rates.

No other count distribution has both of these properties.

The custom N/z and activity center samplers use the approach of 
[Herliansyah et al. (2024), Section 4.3](https://link.springer.com/article/10.1007/s13253-023-00598-3) so that the
total unidentified detection intensity does not have to be completely resummed across all individuals for every
proposed update.

## Mark states

For the models in which the number of marked individuals is known, the standard models assume that mark state
is constant over the sighting period. The `Interspersed` versions allow each individual's mark status to vary
by occasion, but these states must be known. This requirement limits the types of marks that can be used. One example where this information can be
known is when marks are GPS collars and the individual identity of detections is determined by comparing GPS locations
to sighting events. Researchers generally know when collars die, fall off, or stop collecting locations via the VHF beacon,
mortality signals, remote communication with the collar, or telemetry data from remote communication or collar recovery,
however, this is often not perfectly known. These factors need to be considered carefully before using SMR.

The interspersed approach is also required to account for known deaths during the resighting period because these individuals
should be removed from the sighting process after death. While more flexible, the 
interspersed version requires modeling the occasion dimension and therefore can be much slower, depending on the number of occasions.
One might choose to aggregate occasions to speed up run time at the cost of less accuracy of the mark states relative to sightings.


### Multiple marked states and GPS collar failure

One feature of using GPS collars as marks is that a dead GPS collar may not provide a marked-ID
observation type. Then, if one uses the same thinning rates as for marked individuals with functioning collars, the 
thinning rates will not be the same and bias will be introduced (though it may be negligible depending on the scenario).
An approach to accommodate this feature of the data is to treat these individuals as unmarked when their
collar fails and classify these detection types as "unmarked". This discards data, but avoids bias. 

A better approach is to include a second marked state that represents "dead GPS collar", and estimate the
thinning rates separately by marked class. The `GPSfail` model extends the binary marked/unmarked representation
to distinguish three states:

1. unmarked,
2. marked with a working GPS collar, and
3. marked with a failed GPS collar.

A working GPS collar can produce all three marked observation types,

$$
\boldsymbol{\theta}^{M_1}=\left(\theta_{\mathrm{ID}}^{M_1},\theta_{\mathrm{noID}}^{M_1},\theta_{\mathrm{unk}}^{M_1}\right).
$$

A failed GPS collar can no longer produce an individually identifiable observation, so

$$
\boldsymbol{\theta}^{M_2}=\left(0,\theta_{\mathrm{noID}}^{M_2},\theta_{\mathrm{unk}}^{M_2}\right).
$$

Truly unmarked individuals similarly cannot produce an individually identifiable observation,

$$
\boldsymbol{\theta}^{U}=\left(0,\theta_{\mathrm{um}}^{U},\theta_{\mathrm{unk}}^{U}\right).
$$

This avoids treating an individual whose GPS unit has failed as unmarked and losing information. The approach is more
general than the GPS failure example and can be used for two or more mark types with different sample type observation
probabilities, such as GPS collars and ear tags.

## Integrated telemetry data

Telemetry locations can be included directly in the spatial likelihood to inform activity centers and $\sigma$.

For telemetry location $l$ from individual $i$, the two spatial coordinates are modeled as

$$
L_{i,l,d}\sim\mathrm{Normal}\left(s_{i,d},\sigma\right),\qquad d=1,2.
$$

Thus, telemetry locations are centered on the same activity center used by the SMR observation model and use the same
spatial scale parameter $\sigma$. In multisession models, telemetry locations are associated with the activity center and spatial scale parameter
for the corresponding session,

$$
L_{i,l,d}\sim\mathrm{Normal}\left(s_{i,g,d},\sigma_g\right).
$$

The supplied code links the telemetry observations to the relevant individual and, for multisession models,
the relevant session using  `tel.ID` and `tel.session`.

## Model versions

### Known number of marked individuals

These models assume that the number of marked individuals in the population is known, following
[Chandler and Royle (2013)](https://www.jstor.org/stable/23566419) and 
[Sollmann et al. (2013)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-1256.1).

1. **Single session — `SMR Poisson Dcov DA2 Marginal`**

   Standard closed-population marginalized SMR model with a known number of marked individuals. Mark state is constant
   during the sighting period and observations are summed over occasions.

2. **Multisession — `SMR Multisession Poisson Dcov Marginal`**

   Multisession version of the same model. Each session has its own closed-population, abundance,
   and activity centers. The supplied example shares detection and sample type parameters across sessions.

### Known number of marked individuals with interspersed marking and sighting

These models allow mark status to change within a closed-population session but do not include the marking process
in the likelihood. They are analogous to the interspersed design considered by
[Whittington et al. (2018)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954) without
explicitly modeling capture.

1. **Single session — `SMR Poisson Dcov DA2 Marginal Interspersed`**

   Retains sighting observations and mark states by occasion so individuals can contribute unmarked 
   sightings before marking and marked sightings afterward.

2. **Multisession — `SMR Multisession Poisson Dcov Marginal Interspersed`**

   Multisession version with both session and occasion dimensions retained.

### Unknown number of marked individuals

These models are intended for natural marks [(e.g., Rich et al. 2014)](https://academic.oup.com/jmammal/article/95/2/382/866592)
or researcher-deployed marks when the current number of marked individuals
is unknown, [(e.g., Rutledge et al. 2015)](https://connectsci.au/wr/article/41/5/447/40589/Using-novel-spatial-mark-resight-techniques-to).
Marked and unmarked abundance are both estimated. The supplied models use separate baseline density parameters
for the marked and unmarked components while sharing the spatial density covariate effect.

1. **Single session — `SMR Poisson Dcov DA2 Marginal Natural`**

   Estimates marked and unmarked abundance simultaneously within a single closed session.

2. **Multisession — `SMR Multisession Poisson Dcov Marginal Natural`**

   Multisession version estimating the marked and unmarked components separately within each session.

### Generalized SMR with a known number of marked individuals

Generalized SMR (gSMR) explicitly models the marking process to account for different spatial distributions
of marked and unmarked individuals caused by spatially nonrandom marking effort, following
[Whittington et al. (2018)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954).

1. **Single session — `SMR Poisson Dcov DA2 Marginal Generalized`**

   Standard gSMR model with a binomial marking process followed by the marginalized Poisson sighting model.

2. **Single session Mb — `SMR Poisson Dcov DA2 Marginal Generalized Mb`**

   Modification of the single-session gSMR model that allows baseline marking detection probability to differ
   between first capture and subsequent captures.

3. **Multisession — `SMR Multisession Poisson Dcov Marginal Generalized`**

   Multisession gSMR model with a separate marking and sighting process in each closed session.

### Generalized SMR with interspersed marking and sighting

These models combine the explicit gSMR marking process with mark states that vary over occasions, allowing marking and
sighting to be interspersed within a closed session.

1. **Single session — `SMR Poisson Dcov DA2 Marginal Generalized Interspersed`**

   Explicitly models the marking process and retains sighting observations and mark states by occasion.

2. **Multisession — `SMR Multisession Poisson Dcov Marginal Generalized Interspersed`**

   Multisession version of the generalized interspersed model.

3. **Single session with two marked states — `SMR Poisson Dcov DA2 Marginal Generalized Interspersed GPSfail`**

   Extends the generalized interspersed model to allow two marked states with different sample type observation
   probabilities. The supplied simulator is set up for working versus failed GPS collars, but the model can represent
   more general scenarios with two mark types.

### One-stage SMR

The One-stage SMR models use the marked individual data twice, following
[Whittington et al. (2025)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.70246).
This approach provides a means of accounting for different spatial distributions of marked and unmarked individuals 
without explicitly modeling the marking process as in gSMR. It is a modification of the two-stage approach of
[Margenau et al. (2022)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/eap.2553) that models both 
stages simultaneously.

The implementation here differs from Whittington et al. (2025) by allowing marked-no-ID detections
to be included. Marked individual detections are first used to estimate the marked individual spatial distribution, 
with `theta.thin` describing the probability that a marked sighting retains individual identity. All detections are
then used again in the population-level component, treating individuals as unmarked. 

Margenau et al. (2022) state:

> "In the instance a marked individual cannot be reliably recognized in camera photographs, perhaps due to body
positioning, vegetation obstruction, or blurriness, the data record should be discarded from stage one of the model
but included as a record in stage two."

However, discarding these observations from stage one will introduce bias, whether this recommendation refers to
marked-no-ID detections or unknown-mark-status detections. Using the random thinning model
in stage one solves the problem for marked-no-ID samples, but there is currently no solution for unknown-mark-status
samples generally. In the specific case where $\theta_{\mathrm{unk}}^{M}=\theta_{\mathrm{unk}}^{U}$, these samples
can safely be discarded from both stages.

Because the marked-individual data are used twice, there is also a potential concern that this approach could
produce biased parameter estimates or underestimate posterior uncertainty. Simulations in Margenau et al. (2022)
and Whittington et al. (2025) showed minimal bias and approximately nominal coverage when considering only
marked-ID and unmarked sample types and no spatial variation in density. This is reassuring, but I would not
conclude that using the data twice in this manner does not introduce bias or less than nominal coverage in all 
scenarios.

1. **Single session — `SMR Poisson Dcov DA2 Marginal OneStage`**

   Single session implementation of the one-stage approach.

2. **Multisession — `SMR Multisession Poisson Dcov Marginal OneStage`**

   Multisession implementation in which the one-stage approach is applied separately within each closed session.

