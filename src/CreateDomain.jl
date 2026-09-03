#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

function _cell_count(extent, res; rtol=1e-9)
    res > 0 || throw(ArgumentError("resolution must be positive, got $res"))
    n = round(Int, extent / res)
    n >= 1 || throw(ArgumentError("domain size $extent is smaller than resolution $res"))
    isapprox(n * res, extent; rtol=rtol) ||
        @warn "Resolution does not divide domain size" extent res leftover=(extent - n*res)
    return n
end

struct FITDomain{T<:AbstractFloat, A<:AbstractVector{T}, M<:AbstractMatrix{T}}
    # Grid nodes
    nodes_u::A        # length Nu
    nodes_v::A        # length Nv
    nodes_w::A        # length Nw
    # Primary edge lengths (one fewer than nodes)
    edges_u::A        # length Nu-1
    edges_v::A        # length Nv-1
    edges_w::A        # length Nw-1
    # Dual edge lengths
    dual_edges_u::A        # length Nu
    dual_edges_v::A        # length Nv
    dual_edges_w::A        # length Nw
    # Primary edge centres
    edges_u_center::A # length Nu-1
    edges_v_center::A # length Nv-1
    edges_w_center::A # length Nw-1
    # Primal facet size
    primal_facets_u::M      
    primal_facets_v::M 
    primal_facets_w::M 
    # Dual facet size
    dual_facets_u::M
    dual_facets_v::M
    dual_facets_w::M
    # Grid dimensions
    Nu::Int           # number of nodes in u
    Nv::Int           # number of nodes in v
    Nw::Int           # number of nodes in w
    Np::Int           # = Nu*Nv*Nw (total number of nodes)
    # Homogeneous background material
    σ::T
    ε_r::T
    μ_r::T
    # Per-cell material distribution ((Nu-1) × (Nv-1) × (Nw-1) × 3)
    # index 4: 1 = σ, 2 = ε_r, 3 = μ_r
    material::Array{T, 4}

        function FITDomain(
            nodes_u::A, nodes_v::A, nodes_w::A,
            edges_u::A, edges_v::A, edges_w::A,
            dual_edges_u::A, dual_edges_v::A, dual_edges_w::A,
            edges_u_center::A, edges_v_center::A, edges_w_center::A,
            primal_facets_u::M, primal_facets_v::M, primal_facets_w::M,
            dual_facets_u::M, dual_facets_v::M, dual_facets_w::M,
            σ::T, ε_r::T, μ_r::T,
            material::Array{T, 4}
        ) where {T<:AbstractFloat, A<:AbstractVector{T}, M<:AbstractMatrix{T}}

        Nu = length(nodes_u)
        Nv = length(nodes_v)
        Nw = length(nodes_w)
        Np = Nu * Nv * Nw

        length(edges_u)        == Nu-1 || throw(ArgumentError("edges_u must have length Nu-1 = $(Nu-1)"))
        length(edges_v)        == Nv-1 || throw(ArgumentError("edges_v must have length Nv-1 = $(Nv-1)"))
        length(edges_w)        == Nw-1 || throw(ArgumentError("edges_w must have length Nw-1 = $(Nw-1)"))
        length(dual_edges_u) == Nu || throw(ArgumentError("dual_edges_u must have length Nu = $Nu"))
        length(dual_edges_v) == Nv || throw(ArgumentError("dual_edges_v must have length Nv = $Nv"))
        length(dual_edges_w) == Nw || throw(ArgumentError("dual_edges_w must have length Nw = $Nw"))
        length(edges_u_center) == Nu-1 || throw(ArgumentError("edges_u_center must have length Nu-1 = $(Nu-1)"))
        length(edges_v_center) == Nv-1 || throw(ArgumentError("edges_v_center must have length Nv-1 = $(Nv-1)"))
        length(edges_w_center) == Nw-1 || throw(ArgumentError("edges_w_center must have length Nw-1 = $(Nw-1)"))
        size(primal_facets_u) == (Nv-1, Nw-1) || throw(ArgumentError("primal_facets_u must be (Nv-1)×(Nw-1)"))
        size(primal_facets_v) == (Nu-1, Nw-1) || throw(ArgumentError("primal_facets_v must be (Nu-1)×(Nw-1)"))
        size(primal_facets_w) == (Nu-1, Nv-1) || throw(ArgumentError("primal_facets_w must be (Nu-1)×(Nv-1)"))
        size(dual_facets_u)   == (Nv,   Nw)   || throw(ArgumentError("dual_facets_u must be Nv×Nw"))
        size(dual_facets_v)   == (Nu,   Nw)   || throw(ArgumentError("dual_facets_v must be Nu×Nw"))
        size(dual_facets_w)   == (Nu,   Nv)   || throw(ArgumentError("dual_facets_w must be Nu×Nv"))
                
        size(material) == (Nu-1, Nv-1, Nw-1, 3) ||
            throw(ArgumentError("material must be ((Nu-1)×(Nv-1)×(Nw-1)×3), got $(size(material))"))

        new{T, A, M}(
            nodes_u, nodes_v, nodes_w,
            edges_u, edges_v, edges_w,
            dual_edges_u, dual_edges_v, dual_edges_w,
            edges_u_center, edges_v_center, edges_w_center,   # ← move before facets
            primal_facets_u, primal_facets_v, primal_facets_w,
            dual_facets_u, dual_facets_v, dual_facets_w,
            Nu, Nv, Nw, Np,
            σ, ε_r, μ_r,
            material
        )
    end
end

function create_domain(domain_size, resolution; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r  = Float64(μ_r)

    if !(length(domain_size) in (3)) || length(resolution) != 3
        @warn "The length of the domain vector (3) or resolution (3) is wrong"
        return NaN
    end

    o = zeros(6)
    if length(domain_size) == 3
        o[2] = domain_size[1]
        o[4] = domain_size[2]
        o[6] = domain_size[3]
    else
        o=copy(domain_size)
    end

    _create_domain(o[1],o[2],o[3],o[4],o[5],o[6],resolution[1],resolution[2],resolution[3]; units=units, σ=σ, ε_r=ε_r, μ_r=μ_r)
end

#this is a basic implementation for equi-distant grids
function _create_domain(x_min, x_max, y_min, y_max, z_min, z_max, res_x, res_y, res_z; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r  = Float64(μ_r)

    unitToMeter = check_units(units)

    n_x = _cell_count(x_max - x_min, res_x)
    n_y = _cell_count(y_max - y_min, res_y)
    n_z = _cell_count(z_max - z_min, res_z)

    #virtually shift the edges to zero,zero,zero origin and switch to SI units
    Edges_U = fill(res_x * unitToMeter, n_x)
    Edges_V = fill(res_y * unitToMeter, n_y)
    Edges_W = fill(res_z * unitToMeter, n_z)

    return create_domain(Edges_U,Edges_V,Edges_W;units="m",σ,ε_r,μ_r)
end            

#here you can specify non-equi-distant grids by yourself. Important this takes Edges as input! Input needs to in meters!
function create_domain(Edges_U, Edges_V, Edges_W; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r  = Float64(μ_r)
    unitToMeter = check_units(units)

    Edges_U = vec(Float64.(Edges_U))
    Edges_V = vec(Float64.(Edges_V))
    Edges_W = vec(Float64.(Edges_W))

    Edges_U.*=unitToMeter
    Edges_V.*=unitToMeter
    Edges_W.*=unitToMeter

    Dual_Edges_U = vcat(Edges_U[1]/2.0, (Edges_U[1:end-1] .+ Edges_U[2:end])/2.0, Edges_U[end]/2.0)
    Dual_Edges_V = vcat(Edges_V[1]/2.0, (Edges_V[1:end-1] .+ Edges_V[2:end])/2.0, Edges_V[end]/2.0)
    Dual_Edges_W = vcat(Edges_W[1]/2.0, (Edges_W[1:end-1] .+ Edges_W[2:end])/2.0, Edges_W[end]/2.0)

    Primal_Facets_U = Edges_V  .* Edges_W'   # (Nv-1) × (Nw-1), normal to u
    Primal_Facets_V = Edges_U  .* Edges_W'   # (Nu-1) × (Nw-1), normal to v
    Primal_Facets_W = Edges_U  .* Edges_V'   # (Nu-1) × (Nv-1), normal to w

    Dual_Facets_U = Dual_Edges_V .* Dual_Edges_W'   # Nv × Nw
    Dual_Facets_V = Dual_Edges_U .* Dual_Edges_W'   # Nu × Nw
    Dual_Facets_W = Dual_Edges_U .* Dual_Edges_V'   # Nu × Nv

    Nodes_U = cumsum(vcat(0.0, Edges_U))
    Nodes_V = cumsum(vcat(0.0, Edges_V))
    Nodes_W = cumsum(vcat(0.0, Edges_W))

    Edges_U_Center = Nodes_U[1:end-1] + (Nodes_U[2:end]-Nodes_U[1:end-1])/2.0
    Edges_V_Center = Nodes_V[1:end-1] + (Nodes_V[2:end]-Nodes_V[1:end-1])/2.0
    Edges_W_Center = Nodes_W[1:end-1] + (Nodes_W[2:end]-Nodes_W[1:end-1])/2.0


    Nu = length(Nodes_U)
    Nv = length(Nodes_V)
    Nw = length(Nodes_W)
    material = Array{Float64}(undef, Nu-1, Nv-1, Nw-1, 3)
    #the ordering is just randomly selected - each primary cell gets three material properties
    material[:,:,:,1] .= σ
    material[:,:,:,2] .= ε_r
    material[:,:,:,3] .= μ_r

    domain = FITDomain(
        Nodes_U, Nodes_V, Nodes_W,
        Edges_U, Edges_V, Edges_W,
        Dual_Edges_U, Dual_Edges_V, Dual_Edges_W,
        Edges_U_Center, Edges_V_Center, Edges_W_Center,      # ← centres first
        Primal_Facets_U, Primal_Facets_V, Primal_Facets_W,
        Dual_Facets_U, Dual_Facets_V, Dual_Facets_W,
        σ, ε_r, μ_r,
        material
    )
    return domain
end            