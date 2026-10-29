#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 22/03/2026

"""
    get_vacuum_cfl_time(config) -> Float64

The time step, in s, at the stability limit of the explicit leapfrog scheme in vacuum
on the grid of `config`, the Courant–Friedrichs–Lewy condition of the Yee scheme:
`Δt = 1 / (c₀ √(1/Δu² + 1/Δv² + 1/Δw²))`, with the smallest cell width in each
direction and `c₀ = 1/√(ε₀ μ₀)`.

On an equidistant grid this is exactly the stability limit, where the scheme is only
marginally stable: take a slightly smaller step, e.g. `0.99 * get_vacuum_cfl_time(config)`.
On a graded grid it lies below the exact limit, on the safe side, since the smallest
widths of the three directions rarely meet in one cell. Materials with `ε_r μ_r > 1`
slow the waves down and would allow a larger step; they are not taken into account.
"""
function get_vacuum_cfl_time(config)
    numerator = sqrt(ε₀ * μ₀)

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