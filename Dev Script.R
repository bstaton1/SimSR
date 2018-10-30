
params = init_sim(U_msy = rbeta(12, 30, 70),
                  S_msy = rlnorm(12, log(30000), 0.2))

true = ricker_sim(params)
obs = obs_sim(params, true)
obs_filter(params, obs)
