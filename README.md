# SimSR

This is an R package written to facilitate simulating Ricker spawner-recruit dynamics, particularly in the case of age-structured, mixed-stock fisheries. It was designed to generate input information for used by the functions in the `FitSR` [package](<https://github.com/bstaton1/FitSR/>). It may be installed using:

```R
devtools::install_github("bstaton1/SimSR")
```

Specific versions can be installed using:

```R
devtools::install_github("bstaton1/SimSR@v0.1.3")
```

### Example Usage

The functions are all reasonably-well documented, see the help pages for more details on what each one does, or the source code for exact details.

The general workflow is laid out as follows:

1. Generate the parameters to drive the populations
2. Simulate the true dynamics of the populations, and store the states at each time step
3. Simulate the sampling of the populations

#### Step 1: Generate the true parameters

Use the `init_sim()` function to create a list of dimensional variables and driving parameters for use in many of the other functions in this package. Only the arguments `U_msy` and `S_msy` do not have defaults:

```R
params = init_sim(
  U_msy = c(0.64, 0.46, 0.55, 0.37, 0.68, 0.5, 0.43, 0.5, 0.71, 0.61, 0.43, 0.53),
  S_msy = c(8000, 3300, 10800, 1100, 2600, 6300, 7300, 600, 2500, 1700, 300, 1100)
)
names(params)
```

The `params` object is a list storing the information needed to simulate the true and observed states.

#### Step 2: Simulate population dynamics

Use the `ricker_sim()` function to drive the populations forward in time. A simple harvest control rule is used: constant exploitation rate (with implementation error) equal among populations and set at the level that will maximize yield without overfishing more than p*100% of the populations, where p is either 0.1, 0.3, or 0.5.

```R
true = ricker_sim(params, Ustrategy = "Ustar_0.5")
```

The `true` object is a list storing the true states, including recruitment by brood year and stock (`R_ys`), the realized explotiation rate (`U_real`), total run, escapement, and age composition by calendar year and stock (`N_ts`, `S_ts`, and `q_tas`, respectively). 

#### Step 3: Simulate sampling the populations

Sampling proceeds in two steps: first, introduce measurement errors (the state that would be observed if that quantity was monitored for that stock and year), and second, impose a frequency of sampling. To introduce measurement errors in (_i_) population-specific escapement, (_ii_) population-specific age composition, and (_iii_) aggregate harvest, use the `obs_sim()` function:

```
obs = obs_sim(params, true)
```

The `obs` object contains the states that are possibly observable, with inserted measurement errors.

To impose a sampling frequency, use the `obs_filter()` function:

```R
obs = obs_filter(params, obs, mimic = NULL, minSRobs = 3, p_age = 1)
```

This will impose a sampling frequency where each population has escapement monitored for 60% of the years, 100% of the populations are sampled for age composition, and a minimum of 3 observed pairs of spawners and recruits are available for each population. The sampling frequency can be structured to mimic a specified schedule, see `?obs_filter` for more details.

### Notes

The functions are not set up to handle age structures with fewer than possible ages of maturity, but can handle scenarios with any number of ages greater than this. 

The functions used for observing populations are not set up for fewer than two populations. The simulation of true dynamics can handle only one population.

Version `v0.1.3` was used in the analyses for B. Staton's dissertation and the subsequent manuscript.
