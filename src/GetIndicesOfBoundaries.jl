#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 07/08/2026

using SparseArrays


"""
    DomainFace

A side of the domain along a `Direction`: `Positive()`, the face on the last node plane,
or `Negative()`, the face on the first; see `DomainBoundary`.
"""
abstract type DomainFace end

"""
    Positive()

The face of the domain on the last node plane of a direction, e.g. `(DirX(), Positive())`
at u = `nodes_u[end]`; see `DomainBoundary`.
"""
struct Positive <: DomainFace end

"""
    Negative()

The face of the domain on the first node plane of a direction, e.g. `(DirZ(), Negative())`
at w = `nodes_w[1]`; see `DomainBoundary`.
"""
struct Negative <: DomainFace end


"""
    BoundaryComponent

Which entries on a face the boundary functions select: `NodalComponent()` for nodal
vectors of length `Np`, `NormalComponent()` or `TangentialComponent()` for edge and facet
vectors of length `3Np`.
"""
abstract type BoundaryComponent end

"""
    NormalComponent()

The component along the normal of a face: on a face normal to u, the u-entries at its
nodes. For edges these leave the face, and on the positive faces they leave the domain
and are dead (see `get_ghost_indices`); for facets they lie in the face, e.g. the flux
`B·n` through it.
"""
struct NormalComponent     <: BoundaryComponent end

"""
    TangentialComponent()

The two components in the plane of a face: on a face normal to u, the v- and w-entries
at its nodes. For edges these lie in the face, e.g. the edges of a PEC wall, where
`n × E = 0`; for facets they stand on it, and on the positive faces they leave the
domain and are dead.
"""
struct TangentialComponent <: BoundaryComponent end

"""
    NodalComponent()

The nodes of a face, for nodal vectors of length `Np` such as a potential, e.g. for a
grounded wall.
"""
struct NodalComponent <: BoundaryComponent end

"""
    DomainBoundary

A face of the domain as a tuple `(normal, side)`, e.g. `(DirX(), Positive())`. The
boundary functions take a vector of them:

```julia
faces = [(n, s) for n in (DirX(), DirY(), DirZ()) for s in (Positive(), Negative())]   # all six
```
"""
const DomainBoundary = Tuple{Direction, DomainFace}

# --- Nodal Components ---
get_component_offsets(::Direction, ::NodalComponent, Np) = (0,)

# --- Normal Components ---
get_component_offsets(::DirX, ::NormalComponent, Np) = (0,)
get_component_offsets(::DirY, ::NormalComponent, Np) = (Np,)
get_component_offsets(::DirZ, ::NormalComponent, Np) = (2*Np,)

# --- Tangential Components ---
get_component_offsets(::DirX, ::TangentialComponent, Np) = (Np, 2*Np)
get_component_offsets(::DirY, ::TangentialComponent, Np) = (0, 2*Np)
get_component_offsets(::DirZ, ::TangentialComponent, Np) = (0, Np)


function _create_projection_matrix(indices::Vector{Int}, total_size::Int)
    diag_vec = ones(Float64, total_size)
    diag_vec[indices] .= 0.0
    
    return spdiagm(0 => diag_vec)
end

"""
    get_boundary_matrix(config, faces, component) -> SparseMatrixCSC

A diagonal projection that zeroes the entries of `component` on the `faces`, a vector of
`(normal, side)` tuples (see `DomainBoundary`), and keeps all others: `Np × Np` for
`NodalComponent()`, `3Np × 3Np` for `NormalComponent()` and `TangentialComponent()`.
With `R` this matrix and `D = I - R`, `R*K*R + D` replaces the rows and columns of the
boundary entries of a system matrix `K` by the identity: homogeneous Dirichlet
conditions on those entries.

# Example

```julia
faces = [(n, s) for n in (DirX(), DirY(), DirZ()) for s in (Positive(), Negative())]
R = get_boundary_matrix(domain, faces, NodalComponent())    # every wall grounded
D = Diagonal(ones(domain.Np)) - R
Φ = cholesky(D + R*L*R) \\ (R*q)                               # L = Gᵀ M_ε G
```
"""
function get_boundary_matrix(config::FITDomain, faces::Vector{<:DomainBoundary}, component::BoundaryComponent)
    indices = get_custom_boundary_indices(config, faces, component)
    total_size = component isa NodalComponent ? config.Np : 3 * config.Np
    return _create_projection_matrix(indices, total_size)
end

# --- Ghost / Dead Component Matrix ---
"""
    get_ghost_matrix(config) -> SparseMatrixCSC

A diagonal `3Np × 3Np` projection that zeroes the dead entries of an edge vector, those
of `get_ghost_indices`, and keeps all others, e.g. to clean `G * Φ` before reading it as
a field. System matrices assembled with the material matrices do not need it, since
those are zero on the dead edges.
"""
function get_ghost_matrix(config::FITDomain)
    indices = get_ghost_indices(config)
    return _create_projection_matrix(indices, 3*config.Np)
end

"""
    get_ghost_indices(config) -> Vector{Int}

The dead entries of an edge vector, sorted: the edges that would leave the last node
plane of each direction, `Nv Nw + Nu Nw + Nu Nv` of them, i.e. the normal components on
the three positive faces. Their rows of `get_gradient(config, Primal())` difference
unrelated nodes, and `get_permittivity` and `get_conductivity` are exactly zero on them.
For facet vectors the dead entries are the tangential components on the positive faces
instead (see `TangentialComponent`).
"""
function get_ghost_indices(config::FITDomain)
    #indentify dead edges/facets
    idx_X = get_boundary_indices(config, DirX(), Positive(), NormalComponent())[1][1]
    idx_Y = get_boundary_indices(config, DirY(), Positive(), NormalComponent())[1][1]
    idx_Z = get_boundary_indices(config, DirZ(), Positive(), NormalComponent())[1][1]
    
    indices = sort!(unique!(vcat(idx_X, idx_Y, idx_Z)))
        
    return indices
end

"""
    get_all_tangential_boundary_indices(config) -> Vector{Int}

The tangential entries of an edge or facet vector on all six faces of the domain, sorted.
For edges these are the edges in the walls, e.g. to impose `n × E = 0` on a closed PEC
box, or `n × A = 0`, no flux through the walls, in magnetostatics.
"""
function get_all_tangential_boundary_indices(config::FITDomain)
    tx_pos = get_boundary_indices(config, DirX(), Positive(), TangentialComponent())[1]
    tx_neg = get_boundary_indices(config, DirX(), Negative(), TangentialComponent())[1]
    
    ty_pos = get_boundary_indices(config, DirY(), Positive(), TangentialComponent())[1]
    ty_neg = get_boundary_indices(config, DirY(), Negative(), TangentialComponent())[1]
    
    tz_pos = get_boundary_indices(config, DirZ(), Positive(), TangentialComponent())[1]
    tz_neg = get_boundary_indices(config, DirZ(), Negative(), TangentialComponent())[1]
    
    indices = vcat(
        tx_pos..., tx_neg...,
        ty_pos..., ty_neg...,
        tz_pos..., tz_neg...
    )
    
    return sort!(unique!(indices))
end

"""
    get_all_normal_boundary_indices(config) -> Vector{Int}

The normal entries of an edge or facet vector on all six faces of the domain, sorted.
For facets these are the fluxes through the walls; for edges, the edges that leave the
walls (dead on the positive faces).
"""
function get_all_normal_boundary_indices(config::FITDomain)
    nx_pos = get_boundary_indices(config, DirX(), Positive(), NormalComponent())[1]
    nx_neg = get_boundary_indices(config, DirX(), Negative(), NormalComponent())[1]
    
    ny_pos = get_boundary_indices(config, DirY(), Positive(), NormalComponent())[1]
    ny_neg = get_boundary_indices(config, DirY(), Negative(), NormalComponent())[1]
    
    nz_pos = get_boundary_indices(config, DirZ(), Positive(), NormalComponent())[1]
    nz_neg = get_boundary_indices(config, DirZ(), Negative(), NormalComponent())[1]
    
    indices = vcat(
        nx_pos..., nx_neg...,
        ny_pos..., ny_neg...,
        nz_pos..., nz_neg...
    )
    
    return sort!(unique!(indices))
end

"""
    get_boundary_indices(config, normal, side) -> (p, i, j, k)
    get_boundary_indices(config, normal, side, component) -> (indices, i, j, k)

The entries on one face of the domain, `normal` being `DirX()`, `DirY()` or `DirZ()` and
`side` `Positive()` or `Negative()`. The first form gives the indices `p` of the face's
nodes in `1:Np` and the grid indices `i`, `j`, `k` that span the face: a number along
the normal, ranges across it. The second gives, for `component`, a tuple of index
vectors: one into a nodal vector for `NodalComponent()`, one into an edge or facet
vector of length `3Np` for `NormalComponent()`, and two, one per in-plane component, for
`TangentialComponent()`.
"""
function get_boundary_indices(config::FITDomain, normals::Direction, side::DomainFace, component::BoundaryComponent)
    p_base, i, j, k = get_boundary_indices(config, normals, side)
    
    offsets = get_component_offsets(normals, component, config.Np)
    
    p_components = Tuple(p_base .+ offset for offset in offsets)
    
    return p_components, i, j, k
end

"""
    get_custom_boundary_indices(config, faces, component) -> Vector{Int}

The entries of `component` on all the `faces`, a vector of `(normal, side)` tuples (see
`DomainBoundary`), sorted and without repetitions: the indices that
`get_boundary_matrix` zeroes, e.g. the tangential edges of the walls that are PEC.
"""
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

# 1. DirX Boundaries (+u, -u)
function get_boundary_indices(config::FITDomain, ::DirX, ::Positive)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = config.Nu:config.Nu, 1:config.Nv, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, config.Nu, j, k
end

function get_boundary_indices(config::FITDomain, ::DirX, ::Negative)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:1, 1:config.Nv, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, 1, j, k
end

# 2. DirY Boundaries (+v, -v)
function get_boundary_indices(config::FITDomain, ::DirY, ::Positive)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, config.Nv:config.Nv, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, config.Nv, k
end

function get_boundary_indices(config::FITDomain, ::DirY, ::Negative)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, 1:1, 1:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, 1, k
end

# 3. DirZ Boundaries (+w, -w)
function get_boundary_indices(config::FITDomain, ::DirZ, ::Positive)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, 1:config.Nv, config.Nw:config.Nw
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, j, config.Nw
end

function get_boundary_indices(config::FITDomain, ::DirZ, ::Negative)
    Mu, Mv, Mw = 1, config.Nu, config.Nu * config.Nv
    i, j, k = 1:config.Nu, 1:config.Nv, 1:1
    
    p = _generate_boundary_indices(i, j, k, Mu, Mv, Mw)
    return p, i, j, 1
end