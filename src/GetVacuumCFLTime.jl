#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 22/03/2026

function get_vacuum_cfl_time(config)
    eps_0 = 8.8541878188e-12
    mu_0  = 4π * 1e-7

    numerator = sqrt(eps_0 * mu_0)

    Edge_U = config.edges_u
    Edge_V = config.edges_v
    Edge_W = config.edges_w

    denominator = sqrt(
        maximum(x -> inv(x*x), Edge_U) +
        maximum(x -> inv(x*x), Edge_V) +
        maximum(x -> inv(x*x), Edge_W)
    )
    return numerator / denominator
end