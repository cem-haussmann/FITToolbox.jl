#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

abstract type AbstractFITObject end

abstract type AbstractSolid  <: AbstractFITObject end
abstract type AbstractSource <: AbstractFITObject end
# currents on the primal edges; their edges with the current per unit (±1) are _unit_edges
abstract type AbstractEdgeSource <: AbstractSource end
# closed current loops
abstract type AbstractLoopSource <: AbstractEdgeSource end

struct Material
    σ::Float64
    ε_r::Float64
    μ_r::Float64
    color::Union{Nothing,String}   # for plotting only; nothing = choose automatically

    function Material(σ, ε_r, μ_r, color::Union{Nothing,AbstractString}=nothing)
        _check_material(σ, ε_r, μ_r)
        color === nothing || _valid_color(color) ||
            throw(ArgumentError("unknown color \"$color\": use a hex code like \"#B87333\" \
                                 or one of $(join(sort(collect(keys(_NAMED_COLORS))), ", "))"))
        new(σ, ε_r, μ_r, color === nothing      ? nothing :
                         startswith(color, '#') ? String(color) : lowercase(color))
    end
end

const _NAMED_COLORS = Dict(
    "red"       => "#E41A1C", "lightred"   => "#FB9A99", "darkred"   => "#8B0000",
    "green"     => "#2CA02C", "lightgreen" => "#B2DF8A", "darkgreen" => "#006400",
    "blue"      => "#1F78B4", "lightblue"  => "#A6CEE3", "darkblue"  => "#08306B",
    "orange"    => "#FF7F00", "yellow"     => "#FFD92F", "purple"    => "#6A3D9A",
    "pink"      => "#F781BF", "cyan"       => "#17BECF", "brown"     => "#8C564B",
    "gray"      => "#808080", "lightgray"  => "#D3D3D3", "darkgray"  => "#404040",
    "black"     => "#000000", "white"      => "#FFFFFF",
    # typical conductor colours
    "copper"    => "#B87333", "aluminium"  => "#A8A9AD", "gold"      => "#D4AF37",
    # a cut-out (see _solid_material); view_domain draws it transparent or as a hole
    "background" => "#D3D3D3",
)

# Material values the matrices can take: a negative σ makes a lossy medium gain energy, and
# ε_r or μ_r of zero or below leave M_ε or M_ν singular or indefinite, which surfaces much
# later, or never, as a solver that fails. Checked for every solid and the background.
function _check_material(σ, ε_r, μ_r)
    isfinite(σ) && σ >= 0 ||
        throw(ArgumentError("conductivity σ must be finite and not negative, got $σ"))
    isfinite(ε_r) && ε_r > 0 ||
        throw(ArgumentError("relative permittivity ε_r must be finite and positive, got $ε_r"))
    isfinite(μ_r) && μ_r > 0 ||
        throw(ArgumentError("relative permeability μ_r must be finite and positive, got $μ_r"))
    return nothing
end

const _HEX_COLOR =r"^#(?:[0-9A-Fa-f]{3}|[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$"
_valid_color(c) = occursin(_HEX_COLOR, c) || haskey(_NAMED_COLORS, lowercase(c))

# The material of a new solid. color = "background" makes it a cut-out: it takes the
# background material of the domain, whatever σ, ε_r and μ_r say.
function _solid_material(domain, σ, ε_r, μ_r, color)
    (color !== nothing && lowercase(color) == "background") || return Material(σ, ε_r, μ_r, color)
    (σ, ε_r, μ_r) == (0.0, 1.0, 1.0) ||
        @warn "color = \"background\" makes the object a cut-out of the background material; \
               the given σ, ε_r and μ_r are ignored" σ ε_r μ_r
    return Material(domain.σ, domain.ε_r, domain.μ_r, "background")
end

# hex code for plotting; names are resolved through _NAMED_COLORS
_color_hex(c::String) = startswith(c, '#') ? c : _NAMED_COLORS[lowercase(c)]

# History types; they are needed by the FITDomain struct. The functions are in History.jl
abstract type AbstractOperation end

struct CreateObject <: AbstractOperation
    obj_id::Int
    name::String
    object::AbstractFITObject
end

struct DeleteObject <: AbstractOperation
    obj_id::Int
end

mutable struct DomainHistory
    steps::Vector{AbstractOperation}   # step number = position in this vector
    obj_counter::Int                   # last obj_id handed out
end

DomainHistory() = DomainHistory(AbstractOperation[], 0)


function _cell_count(extent, res; rtol=1e-9)
    isfinite(extent) && extent > 0 ||
        throw(ArgumentError("domain size must be positive and finite, got $extent"))
    res > 0 || throw(ArgumentError("resolution must be positive, got $res"))
    n = round(Int, extent / res)
    n >= 1 || throw(ArgumentError("domain size $extent is smaller than resolution $res"))
    isapprox(n * res, extent; rtol=rtol) ||
        @warn "Resolution does not divide domain size" extent res leftover=(extent - n*res)
    return n
end

"""
    FITDomain

The computational domain made by `create_domain`: the primal grid and its dual, the
background material, the material of every cell, and the history of the objects created
in it. Build it with `create_domain`, not directly.

Its fields, all lengths in m:
- `nodes_u`, `nodes_v`, `nodes_w`: the primal nodes along each axis, `Nu`, `Nv` and `Nw` of
  them; `Np = Nu Nv Nw` nodes in all, numbered `p = i + (j-1)Nu + (k-1)NuNv`.
- `edges_u`, `edges_v`, `edges_w`: the primal edge lengths, i.e. the cell widths, `Nu - 1`
  along u; `edges_u_center`, … their midpoints, where the dual nodes lie.
- `dual_edges_u`, `dual_edges_v`, `dual_edges_w`: the dual edge lengths, from cell centre
  to cell centre and halved at the boundary, `Nu` along u.
- `primal_facets_u`, … and `dual_facets_u`, …: the facet areas normal to each axis,
  `(Nv-1) × (Nw-1)` primal and `Nv × Nw` dual ones normal to u.
- `σ`, `ε_r`, `μ_r`: the background material.
- `material`: `(Nu-1) × (Nv-1) × (Nw-1) × 3`, the conductivity, relative permittivity and
  relative permeability of every cell. It is rebuilt from the history whenever that
  changes; values written into it directly are not part of the history, so
  `remove_object!`, `undo!`, `change_resolution` and `save_domain` warn that they are lost.
"""
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
    _history::DomainHistory

        function FITDomain(
            nodes_u::A, nodes_v::A, nodes_w::A,
            edges_u::A, edges_v::A, edges_w::A,
            dual_edges_u::A, dual_edges_v::A, dual_edges_w::A,
            edges_u_center::A, edges_v_center::A, edges_w_center::A,
            primal_facets_u::M, primal_facets_v::M, primal_facets_w::M,
            dual_facets_u::M, dual_facets_v::M, dual_facets_w::M,
            σ::T, ε_r::T, μ_r::T,
            material::Array{T, 4},
            history::DomainHistory = DomainHistory()
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
            edges_u_center, edges_v_center, edges_w_center,
            primal_facets_u, primal_facets_v, primal_facets_w,
            dual_facets_u, dual_facets_v, dual_facets_w,
            Nu, Nv, Nw, Np,
            σ, ε_r, μ_r,
            material,
            history
        )
    end
end

# Slack for comparing positions with the nodes of one axis. The nodes are a running sum
# of the cell widths, so the last one can miss the requested extent by rounding: ten cells
# of 0.1 m end at 0.9999999999999999 m. Without the slack an object reaching exactly to
# the boundary would count as beyond it.
_slack(nodes) = 1e-9 * (nodes[end] - nodes[1])

# true if x lies within lo … hi of this axis, up to the slack
_within(x, nodes, lo = nodes[1], hi = nodes[end]) =
    lo - _slack(nodes) <= x <= hi + _slack(nodes)

# The cells along `nodes` (cell i from nodes[i] to nodes[i+1]) that can hold part of a
# solid within [lo, hi], its extent along that axis: those overlapping it, and one more on
# each side against rounding. Painting a solid scans these alone, so its cost grows with
# its size on the grid, not with the grid.
_cell_range(nodes, lo, hi) =
    max(1, searchsortedlast(nodes, lo) - 1):min(length(nodes) - 1, searchsortedfirst(nodes, hi) + 1)
_cell_ranges(domain, (lo, hi)) =
    ntuple(d -> _cell_range((domain.nodes_u, domain.nodes_v, domain.nodes_w)[d], lo[d], hi[d]), 3)

# Paint a solid into domain.material, the cells of _fill_cells!. Returns their number.
function _paint!(domain, solid::AbstractSolid)
    σ, ε_r, μ_r = solid.material.σ, solid.material.ε_r, solid.material.μ_r
    md = domain.material
    return _fill_cells!(domain, solid) do i, j, k
        @inbounds md[i, j, k, 1] = σ
        @inbounds md[i, j, k, 2] = ε_r
        @inbounds md[i, j, k, 3] = μ_r
    end
end

# Call set!(i, j, k) for every cell a round solid (sphere, cylinder, torus, helix) takes:
# a cell is taken when at least half of it lies inside, judged by 27 sub-voxel centres at
# 1/6, 3/6, 5/6 of the cell in each direction (14 or more inside). Only the cells around
# the solid (_extent) are scanned. Returns the number of cells. Besides the material, the
# views of view_domain fill their labels with it. The brick has a method of its own.
function _fill_cells!(set!, domain, solid::AbstractSolid)
    # node positions relative to the solid's reference point, the frame of _inside_test
    ref = _reference(solid)
    du = domain.nodes_u .- ref[1]
    dv = domain.nodes_v .- ref[2]
    dw = domain.nodes_w .- ref[3]
    inside = _inside_test(solid)

    # k outermost for threading (stride Nu*Nv, so threads write far apart),
    # i innermost because the cell arrays are contiguous in i
    cells = _cell_ranges(domain, _extent(solid))
    filled = zeros(Int, domain.Nw - 1)
    @threads for k in cells[3]
        ww = (dw[k+1] - dw[k]) / 6.0
        w_samples = (dw[k] + ww, dw[k] + 3.0*ww, dw[k] + 5.0*ww)

        @inbounds for j in cells[2]
            vw = (dv[j+1] - dv[j]) / 6.0
            v_samples = (dv[j] + vw, dv[j] + 3.0*vw, dv[j] + 5.0*vw)

            for i in cells[1]
                uw = (du[i+1] - du[i]) / 6.0
                u_samples = (du[i] + uw, du[i] + 3.0*uw, du[i] + 5.0*uw)

                count_inside = 0
                for us in u_samples, vs in v_samples, ws in w_samples
                    inside((us, vs, ws)) && (count_inside += 1)
                end

                if count_inside >= 14
                    set!(i, j, k)
                    filled[k] += 1
                end
            end
        end
    end
    return sum(filled)
end

# Lengths of solids and loops: positive and finite. An infinite radius would pass a test
# for > 0 and then fill or clip against the whole domain without saying so.
_positive_length(x) = isfinite(x) && x > 0

# A solid must reach into the domain. Its reference point (the corner of a brick, the
# centre of a sphere or torus, the base of a cylinder or helix) may lie outside, e.g. for
# a sphere centred just beyond a wall; the part beyond the domain is clipped, with the
# warning of each solid. Throws for a position that is not finite, and for a solid that
# misses the domain altogether: judged by its bounding box (_extent) along each axis, so a
# box that only touches the domain from outside counts as outside.
function _check_reaches_domain(domain::FITDomain, solid::AbstractSolid, kind, reference)
    ref = _reference(solid)
    all(isfinite, ref) || throw(ArgumentError("$kind $reference must be finite, got $ref m"))
    lo, hi = _extent(solid)
    nodes = (domain.nodes_u, domain.nodes_v, domain.nodes_w)
    for (d, name) in zip(1:3, ("u", "v", "w"))
        lo[d] < nodes[d][end] && hi[d] > nodes[d][1] && continue
        throw(ArgumentError("the $kind lies outside the domain: it spans $name = \
                             $(_short(lo[d])) … $(_short(hi[d])) m, the domain $name = \
                             $(_short(nodes[d][1])) … $(_short(nodes[d][end])) m"))
    end
    return nothing
end

"""
    create_domain(domain_size, resolution; units = "m", σ = 0.0, ε_r = 1.0, μ_r = 1.0)
        -> FITDomain
    create_domain(Edges_U, Edges_V, Edges_W; units = "m", σ = 0.0, ε_r = 1.0, μ_r = 1.0)
        -> FITDomain

A computational domain: a box from the origin with the staggered grid pair of the Finite
Integration Technique on it, filled with a homogeneous background material of
conductivity `σ` (S/m), relative permittivity `ε_r` and relative permeability `μ_r`.
Objects and sources go in with the `create_*!` functions, which record them in the
history of the domain (see `list_objects`).

The first form gives an equidistant grid: `domain_size` and `resolution` are three
values each, the extent and the cell width along u, v and w. A resolution that does not
divide the extent gives the nearest whole number of cells, with a warning; the domain
is then that much smaller or larger. The second form takes the cell widths along each
axis, for a graded grid; the nodes are their running sums from 0. All lengths are given
in `units` ("nm", "μm", "mm", "cm", "dm", "m" or "km"); the domain stores them in m.

Throws an `ArgumentError` if `domain_size` or `resolution` does not have three values, or
if a size, a resolution or a cell width is not positive and finite, if an axis has no
cell, or if the background material is not physical: `σ` must be finite and not
negative, `ε_r` and `μ_r` finite and positive. The solids check their material the same
way.

# Examples

```julia
domain = create_domain([1.0, 1.0, 1.0], [0.01, 0.01, 0.01])        # 100³ cells of 1 cm
domain = create_domain([100, 100, 100], [1, 1, 1]; units = "cm")    # the same
edges  = vcat(20.0, 10.0, 5.0, fill(2.0, 50), 5.0, 10.0, 20.0)      # fine in the middle
domain = create_domain(edges, edges, edges; units = "mm", ε_r = 2.0)
```
"""
function create_domain(domain_size, resolution; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r  = Float64(μ_r)

    length(domain_size) == 3 && length(resolution) == 3 ||
        throw(ArgumentError("domain_size and resolution need 3 values each (u, v, w), got \
                             $(length(domain_size)) and $(length(resolution))"))

    # the domain starts at the origin
    _create_domain(0.0, domain_size[1], 0.0, domain_size[2], 0.0, domain_size[3],
                   resolution[1], resolution[2], resolution[3]; units=units, σ=σ, ε_r=ε_r, μ_r=μ_r)
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

# non-equidistant grids: the cell widths along each axis, in `units`
function create_domain(Edges_U, Edges_V, Edges_W; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r  = Float64(μ_r)
    _check_material(σ, ε_r, μ_r)
    unitToMeter = check_units(units)

    # zero, negative or NaN widths would give nodes out of order, and every search on
    # them silently a wrong cell
    for (name, e) in (("Edges_U", Edges_U), ("Edges_V", Edges_V), ("Edges_W", Edges_W))
        !isempty(e) && all(x -> x isa Real && isfinite(x) && x > 0, e) ||
            throw(ArgumentError("$name must hold at least one cell width, all positive \
                                 and finite, got $(_shown(e))"))
    end

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
        Edges_U_Center, Edges_V_Center, Edges_W_Center,
        Primal_Facets_U, Primal_Facets_V, Primal_Facets_W,
        Dual_Facets_U, Dual_Facets_V, Dual_Facets_W,
        σ, ε_r, μ_r,
        material
    )
    return domain
end

"""
    get_cell_count(domain) -> Int

The number of cells of the primal grid, `(Nu-1) (Nv-1) (Nw-1)`, e.g. to judge the size of
a model before solving it. The cells along each axis are `size(domain.material)[1:3]`;
the domain has `domain.Np` nodes, and about `3Np` edges and facets, the unknowns of the
vector quantities.

# Example

```julia
domain = create_domain([1.0, 0.8, 0.6], [0.01, 0.01, 0.01])
get_cell_count(domain)          # 100 × 80 × 60 = 480000
```
"""
get_cell_count(domain::FITDomain) = (domain.Nu - 1) * (domain.Nv - 1) * (domain.Nw - 1)

# Displayed as a summary: the default show of a struct prints every field, the node
# vectors, facet matrices and material array, thousands of characters on one line.
_cells_label(d::FITDomain) = "$(d.Nu - 1)×$(d.Nv - 1)×$(d.Nw - 1) cells"
_n_objects_label(n) = n == 1 ? "1 object" : "$n objects"

function Base.show(io::IO, d::FITDomain)
    print(io, "FITDomain(", _cells_label(d), ", ",
          _n_objects_label(length(_active_ids(d._history))), ")")
end

function Base.show(io::IO, ::MIME"text/plain", d::FITDomain)
    edges = (d.edges_u, d.edges_v, d.edges_w)
    graded = any(e -> !all(x -> isapprox(x, e[1]; rtol = 1e-9), e), edges)
    n = length(_active_ids(d._history))
    println(io, "FITDomain: ", _cells_label(d), " (", get_cell_count(d), "), ", d.Np, " nodes, ",
            graded ? "graded" : "equidistant")
    println(io, "  extent:     ", join((_short(e[end] - e[1]) for e in
                                       (d.nodes_u, d.nodes_v, d.nodes_w)), " × "), " m")
    graded ||
        println(io, "  resolution: ", join((_short(e[1]) for e in edges), " × "), " m")
    println(io, "  background: σ = ", _short(d.σ), ", ε_r = ", _short(d.ε_r), ", μ_r = ", _short(d.μ_r))
    print(io, "  objects:    ", _n_objects_label(n), n == 0 ? "" : " (see list_objects)")
end
