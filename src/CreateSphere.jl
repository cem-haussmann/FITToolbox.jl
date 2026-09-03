#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 14/01/2026

function create_sphere!(domain::FITDomain, u_o, v_o, w_o, radius; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r = Float64(μ_r)
    unitToMeter = check_units(units)
    start_pos_u = u_o * unitToMeter
    start_pos_v = v_o * unitToMeter
    start_pos_w = w_o * unitToMeter
    r = radius * unitToMeter

    r > 0 || throw(ArgumentError("radius must be positive, got $radius"))

    # The centre must lie inside the grid; otherwise the sphere is missed entirely
    # or only partly captured, and without this check that goes unnoticed.
    (domain.nodes_u[1] <= start_pos_u <= domain.nodes_u[end] &&
     domain.nodes_v[1] <= start_pos_v <= domain.nodes_v[end] &&
     domain.nodes_w[1] <= start_pos_w <= domain.nodes_w[end]) ||
        throw(ArgumentError("sphere centre ($u_o, $v_o, $w_o) $units lies outside the domain"))

    # A sphere reaching past the boundary is clipped. That can be intentional —
    # a half sphere at the edge — so warn rather than throw.
    (start_pos_u - r < domain.nodes_u[1] || start_pos_u + r > domain.nodes_u[end] ||
     start_pos_v - r < domain.nodes_v[1] || start_pos_v + r > domain.nodes_v[end] ||
     start_pos_w - r < domain.nodes_w[1] || start_pos_w + r > domain.nodes_w[end]) &&
        @warn "sphere extends beyond the domain and will be clipped" radius=r centre=(start_pos_u, start_pos_v, start_pos_w)

    Nu, Nv, Nw = domain.Nu, domain.Nv, domain.Nw
    md = domain.material

    # Node positions relative to the sphere centre
    du = domain.nodes_u .- start_pos_u
    dv = domain.nodes_v .- start_pos_v
    dw = domain.nodes_w .- start_pos_w

    # k outermost for threading (stride Nu*Nv, so threads write far apart),
    # i innermost because md is contiguous in i.
    @threads for k in 1:Nw-1
        ww = (dw[k+1] - dw[k]) / 6.0
        w_samples = (dw[k] + ww, dw[k] + 3.0*ww, dw[k] + 5.0*ww)

        @inbounds for j in 1:Nv-1
            vw = (dv[j+1] - dv[j]) / 6.0
            v_samples = (dv[j] + vw, dv[j] + 3.0*vw, dv[j] + 5.0*vw)

            for i in 1:Nu-1
                uw = (du[i+1] - du[i]) / 6.0
                u_samples = (du[i] + uw, du[i] + 3.0*uw, du[i] + 5.0*uw)

                # 27 sub-voxel centres at 1/6, 3/6, 5/6 of the cell in each
                # direction; hypot is robust against overflow and underflow.
                count_inside = 0
                for us in u_samples, vs in v_samples, ws in w_samples
                    hypot(us, vs, ws) <= r && (count_inside += 1)
                end

                # 50% fill rule: assign the material if 14 or more are inside
                if count_inside >= 14
                    md[i, j, k, 1] = σ
                    md[i, j, k, 2] = ε_r
                    md[i, j, k, 3] = μ_r
                end
            end
        end
    end

    (any(md[:, :, :, 1] .!= domain.σ)  ||
     any(md[:, :, :, 2] .!= domain.ε_r) ||
     any(md[:, :, :, 3] .!= domain.μ_r)) ||
        @warn "no cells were filled — is the radius smaller than one cell?"
    return domain
end