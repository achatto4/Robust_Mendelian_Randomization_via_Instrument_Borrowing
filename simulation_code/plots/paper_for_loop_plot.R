setwd(".")  # set to repo simulation_code/plots directory
source("Paper_plots_multi.R")


for(overlap_vec in c(0.5, 0.75, 1)){
  for(scenario in c("BI","BN","DI","DN")){
    for(est_theta in c(-0.2,-0.1,0,0.1,0.2)){
      for(thetaUvec in c(0.3)){


analyze_all_scenarios_multi_error(est_theta,
                                scenario,          # Specify the scenario as input
                                thetaUvec,  # Vector
                                Nvec = c("50000", "80000", "100000", "150000", "200000"),
                                prop_invalid_vec = c(0.5, 0.3, 0.7),# Single value
                                overlap_vec,    # Single value
                                sctfc = FALSE)

analyze_all_scenarios_multi_mse(est_theta,
                                scenario,          # Specify the scenario as input
                                thetaUvec,  # Vector
                                Nvec = c("50000", "80000", "100000", "150000", "200000"),
                                prop_invalid_vec = c(0.5, 0.3, 0.7),# Single value
                                overlap_vec,    # Single value
                                sctfc = FALSE)

analyze_all_scenarios_multi_coverage(est_theta,
                                scenario,          # Specify the scenario as input
                                thetaUvec,  # Vector
                                Nvec = c("50000", "80000", "100000", "150000", "200000"),
                                prop_invalid_vec = c(0.5, 0.3, 0.7),# Single value
                                overlap_vec,    # Single value
                                sctfc = FALSE)
      }
    }
  }
}

# setwd(".")  # set to repo simulation_code/plots directory
# source("view_func.R")

# analyse_new_scenarios_eff(power_dir = -1, thetaUvec = c(0.3, 0.5),
#                           Nvec = c("50000", "80000", "100000", "150000", "200000"),
#                           prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7),
#                           overlap_vec = c(0, 0.1, 0.3, 0.5, 0.75, 1),
#                           phi_range = exp(seq(log(0.1), log(10), length.out = 10)),
#                           scenario_vec = c("BI", "BN", "DI", "DN"), color = "overlap",plot_type = 0)
#
# analyse_mrmix_ibmix_power(thetaUvec = c(0.3,0.5),
#                                       Nvec = c("50000", "80000", "100000", "150000", "200000"),
#                                       prop_invalid_vec = c(0.3, 0.5, 0.7),
#                                       overlap_vec = c(0.5, 0.75, 1),
#                                       scenario_vec = c("BI", "BN", "DI", "DN"),
#                                       color_var = "overlap", plot_type = 1)

# analyse_new_scenarios_phi_hist(thetaUvec = c(0.3, 0.5),
#                                            Nvec = c(5e4, 8e4, 1e5, 2e5),
#                                            prop_invalid_vec = c(0.1, 0.3, 0.5, 0.7),
#                                            overlap_vec = c(0, 0.1, 0.3, 0.5, 0.75, 1),
#                                            phi_range = exp(seq(log(0.1), log(50), length.out =30)),
#                                            scenario_vec = c("BI", "BN", "DI", "DN"))