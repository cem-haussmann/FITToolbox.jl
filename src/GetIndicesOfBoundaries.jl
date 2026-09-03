#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 07/08/2026

using SparseArrays


abstract type DomainFace end
struct Positive <: DomainFace end
struct Negative <: DomainFace end


abstract type BoundaryComponent end
struct NormalComponent     <: BoundaryComponent end
struct TangentialComponent <: BoundaryComponent end
struct NodalComponent <: BoundaryComponent end

const DomainBoundary = Tuple{Normals, DomainFace}

# --- Nodal Components ---
get_component_offsets(::Normals, ::NodalComponent, Np) = (0,)

# --- Normal Components ---
get_component_offsets(::X, ::NormalComponent, Np) = (0,)
get_component_offsets(::Y, ::NormalComponent, Np) = (Np,)
get_component_offsets(::Z, ::NormalComponent, Np) = (2*Np,)

# --- Tangential Components ---
get_component_offsets(::X, ::TangentialComponent, Np) = (Np, 2*Np)
get_component_offsets(::Y, ::TangentialComponent, Np) = (0, 2*Np)
get_component_offsets(::Z, ::TangentialComponent, Np) = (0, Np)


function _create_projection_matrix(indices::Vector{Int}, total_size::Int)
    diag_vec = ones(Float64, total_size)
    diag_vec[indices] .= 0.0
    
    return spdiagm(0 => diag_vec)
end

function get_boundary_matrix(config::FITDomain, faces::Vector{<:DomainBoundary}, component::BoundaryComponent)
    indices = get_custom_boundary_indices(config, faces, component)
    total_size = component isa NodalComponent ? config.Np : 3 * config.Np
    return _create_projection_matrix(indices, total_size)
end

# --- Ghost / Dead Component Matrix ---
function get_ghost_matrix(config::FITDomain)
    indices = get_ghost_indices(config)
    return _create_projection_matrix(indices, 3*config.Np)
end

function get_ghost_indices(config::FITDomain)
    #indentify dead edges/facets
    idx_X = get_boundary_indices(config, X(), Positive(), NormalComponent())[1][1]
    idx_Y = get_boundary_indices(config, Y(), Positive(), NormalComponent())[1][1]
    idx_Z = get_boundary_indices(config, Z(), Positive(), NormalComponent())[1][1]
    
    indices = sort!(unique!(vcat(idx_X, idx_Y, idx_Z)))
        
    return indices
end

function get_all_tangential_boundary_indices(config::FITDomain)
    tx_pos = get_boundary_indices(config, X(), Positive(), TangentialComponent())[1]
    tx_neg = get_boundary_indices(config, X(), Negative(), TangentialComponent())[1]
    
    ty_pos = get_boundary_indices(config, Y(), Positive(), TangentialComponent())[1]
    ty_neg = get_boundary_indices(config, Y(), Negative(), TangentialComponent())[1]
    
    tz_pos = get_boundary_indices(config, Z(), Positive(), TangentialComponent())[1]
    tz_neg = get_boundary_indices(config, Z(), Negative(), TangentialComponent())[1]
    
    indices = vcat(
        tx_pos..., tx_neg...,
        ty_pos..., ty_neg...,
        tz_pos..., tz_neg...
    )
    
    return sort!(unique!(indices))
end

function get_all_normal_boundary_indices(config::FITDomain)
    nx_pos = get_boundary_indices(config, X(), Positive(), NormalComponent())[1]
    nx_neg = get_boundary_indices(config, X(), Negative(), NormalComponent())[1]
    
    ny_pos = get_boundary_indices(config, Y(), Positive(), NormalComponent())[1]
    ny_neg = get_boundary_indices(config, Y(), Negative(), NormalComponent())[1]
    
    nz_pos = get_boundary_indices(config, Z(), Positive(), NormalComponent())[1]
    nz_neg = get_boundary_indices(config, Z(), Negative(), NormalComponent())[1]
    
    indices = vcat(
        nx_pos..., nx_neg...,
        ny_pos..., ny_neg...,
        nz_pos..., nz_neg...
    )
    
    return sort!(unique!(indices))
end

function get_boundary_indices(config::FITDomain, normals::Normals, side::DomainFace, component::BoundaryComponent)
    p_base, i, j, k = get_boundary_indices(config, normals, side)
    
    offsets = get_component_offsets(normals, component, config.Np)
    
    p_components = Tuple(p_base .+ offset for offset in offsets)
    
    return p_components, i, j, k
end

function get_custom_boundary_indices(config::FITDomain, faces::Vector{<:DomainBoundary}, component::BoundaryComponent)
    indices = Int[]
    
    for (normal, side) in faces
        comp_tuple = get_boundary_indices(config, normal, side, component)[1]
        
        for arr in comp_tuple
            append!(indices, arr)
        end
    end
    
    return sort!(unique!(indices))
end


function _generate_boundary_indices(i_range, j_range, k_range, Mu, Mv, Mw)
    return vec([1 + (i - 1)*Mu + (j - 1)*Mv + (k - 1)*Mw for i in i_range, j in j_range, k in k_range])
end

# 1. X Boundaries (+u, -u)
function get_boundary_indices(config::FITDomain, ::X, ::Positive)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = config.Nu:config.Nu, 1:config.Nv, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, config.Nu, j, k
end

function get_boundary_indices(config::FITDomain, ::X, ::Negative)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:1, 1:config.Nv, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, 1, j, k
end

# 2. Y Boundaries (+v, -v)
function get_boundary_indices(config::FITDomain, ::Y, ::Positive)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, config.Nv:config.Nv, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, config.Nv, k
end

function get_boundary_indices(config::FITDomain, ::Y, ::Negative)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, 1:1, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, 1, k
end

# 3. Z Boundaries (+w, -w)
function get_boundary_indices(config::FITDomain, ::Z, ::Positive)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, 1:config.Nv, config.Nw:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, j, config.Nw
end

function get_boundary_indices(config::FITDomain, ::Z, ::Negative)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, 1:config.Nv, 1:1
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, j, 1
end