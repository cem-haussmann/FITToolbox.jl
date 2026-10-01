#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025


abstract type AbstractFITObject end

abstract type AbstractSolid  <: AbstractFITObject end
abstract type AbstractSource <: AbstractFITObject end

struct Material
    σ::Float64
    ε_r::Float64
    μ_r::Float64
    color::Union{Nothing,String}   # for plotting only; nothing = choose automatically

    function Material(σ, ε_r, μ_r, color::Union{Nothing,AbstractString}=nothing)
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
)

const _HEX_COLOR = r"^#(?:[0-9A-Fa-f]{3}|[0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$"
_valid_color(c) = occursin(_HEX_COLOR, c) || haskey(_NAMED_COLORS, lowercase(c))

# hex code for plotting; names are resolved through _NAMED_COLORS
_color_hex(c::String) = startswith(c, '#') ? c : _NAMED_COLORS[lowercase(c)]

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
            edges_u_center, edges_v_center, edges_w_center,   # ← move before facets
            primal_facets_u, primal_facets_v, primal_facets_w,
            dual_facets_u, dual_facets_v, dual_facets_w,
            Nu, Nv, Nw, Np,
            σ, ε_r, μ_r,
            material,
            history
        )
    end
end

Base.propertynames(d::FITDomain, private::Bool=false) =
    private ? fieldnames(typeof(d)) : filter(!=(:_history), fieldnames(typeof(d)))

function _record_create!(domain::FITDomain, object::AbstractFITObject; name::String)
    h = domain._history
    h.obj_counter += 1
    push!(h.steps, CreateObject(h.obj_counter, name, object))
    return h.obj_counter
end

function _record_delete!(domain::FITDomain, obj_id::Int)
    h = domain._history
    obj_id in _active_ids(h) ||
        throw(ArgumentError("no object with obj_id $obj_id in the domain"))
    push!(h.steps, DeleteObject(obj_id))
    return nothing
end

function _active_ids(h::DomainHistory)
    ids = Int[]
    for s in h.steps
        s isa CreateObject && push!(ids, s.obj_id)
        s isa DeleteObject && filter!(!=(s.obj_id), ids)
    end
    return ids
end

function _active_entries(h::DomainHistory)
    active = Set(_active_ids(h))
    return CreateObject[s for s in h.steps if s isa CreateObject && s.obj_id in active]
end

const _DEFAULT_NAME = r"^obj\d+$"

# final name for the next object: "obj<id>" if none is given; user names must be unique
function _resolve_name(domain::FITDomain, name::AbstractString)
    h = domain._history
    isempty(name) && return "obj$(h.obj_counter + 1)"
    occursin(_DEFAULT_NAME, name) &&
        throw(ArgumentError("names of the form \"obj<number>\" are reserved for unnamed objects, got \"$name\""))
    any(e -> e.name == name, _active_entries(h)) &&
        throw(ArgumentError("an object named \"$name\" already exists in the domain"))
    return String(name)
end

# the existing object with this name or id
function _entry(domain::FITDomain, name::AbstractString)
    for e in _active_entries(domain._history)
        e.name == name && return e
    end
    throw(ArgumentError("no object named \"$name\" in the domain"))
end

function _entry(domain::FITDomain, obj_id::Int)
    for e in _active_entries(domain._history)
        e.obj_id == obj_id && return e
    end
    throw(ArgumentError("no object with obj_id $obj_id in the domain"))
end

# Rebuild domain.material from the history: background first, then every existing
# solid in creation order, so later objects overwrite earlier ones as before.
# Values written directly into domain.material are lost.
function _replay!(domain::FITDomain; warn::Bool=true)
    warn && @warn "Rebuilding the material from the history: values written directly into \
           domain.material are discarded. Use create_brick!, create_sphere!, … to keep them." maxlog=1
    md = domain.material
    md[:, :, :, 1] .= domain.σ
    md[:, :, :, 2] .= domain.ε_r
    md[:, :, :, 3] .= domain.μ_r

    active = Set(_active_ids(domain._history))
    for step in domain._history.steps
        step isa CreateObject            || continue   # skip DeleteObject steps
        step.obj_id in active            || continue   # skip objects deleted later
        step.object isa AbstractSolid    || continue   # sources do not paint material
        _apply!(domain, step.object)
    end
    return domain
end

"""
    remove_object!(domain, obj_id)

Remove the object `obj_id` from `domain` and rebuild the material from the history.
The removal is recorded as a step of its own, so `undo!` brings the object back.
Throws an `ArgumentError` if no object with that id exists.
"""
function remove_object!(domain::FITDomain, obj_id::Int)
    _record_delete!(domain, obj_id)     # throws for unknown or already deleted ids
    _replay!(domain)
    return nothing
end

remove_object!(domain::FITDomain, name::AbstractString) =
    remove_object!(domain, _entry(domain, name).obj_id)

"""
    undo!(domain)

Revert the last step in the history of `domain`: a created object is removed again,
a removed object comes back. The material is rebuilt from the history. Object ids
are never handed out twice, so the next object still gets a new `obj_id`.

Throws an `ArgumentError` if the history is empty.
"""
function undo!(domain::FITDomain)
    h = domain._history
    isempty(h.steps) &&
        throw(ArgumentError("nothing to undo: the history of this domain is empty"))
    pop!(h.steps)
    _replay!(domain)
    return nothing
end
# true if domain.material is exactly what the history produces,
# i.e. nothing was written into it directly
function _material_matches_history(domain::FITDomain)
    tmp = deepcopy(domain)
    _replay!(tmp; warn=false)
    return tmp.material == domain.material
end

"""
    change_resolution(domain, resolution; units="m") -> FITDomain
    change_resolution(domain, Edges_U, Edges_V, Edges_W; units="m") -> FITDomain

A new domain with the same extent and background material as `domain`, and the
full history of `domain` replayed on it.

The first form gives an equidistant grid: `resolution` is one value for all
directions or three values (u, v, w). The second form takes the cell widths per
direction, as `create_domain(Edges_U, Edges_V, Edges_W)` does, for non-equidistant
grids; each must add up to the extent of `domain` in that direction. All lengths
are given in `units`.

`domain` itself is not changed, and both domains have independent histories from
then on. Values written directly into `domain.material` are not part of the
history; if there are any, a warning says they are not carried over.

Throws an `ArgumentError` if the edges do not match the extent of `domain`, or if
an object cannot be created on the new grid, e.g. a brick thinner than one cell.
Source vectors such as the result of `create_circular_loop_source` are not
updated; compute them again on the new domain.
"""
function change_resolution(domain::FITDomain, Edges_U, Edges_V, Edges_W; units="m")
    unitToMeter = check_units(units)
    for (dir, edges, nodes) in (("u", Edges_U, domain.nodes_u),
                                ("v", Edges_V, domain.nodes_v),
                                ("w", Edges_W, domain.nodes_w))
        L = nodes[end] - nodes[1]
        isapprox(sum(edges) * unitToMeter, L; rtol=1e-9) ||
            throw(ArgumentError("the $dir edges add up to $(sum(edges) * unitToMeter) m, \
                                 but the domain is $L m long in $dir"))
    end
    _material_matches_history(domain) ||
        @warn "domain.material contains values that were not created through create_* \
               functions; they are not carried over to the new grid"

    new = create_domain(Edges_U, Edges_V, Edges_W; units=units,
                        σ=domain.σ, ε_r=domain.ε_r, μ_r=domain.μ_r)
    # copy the history: steps are immutable, so sharing them is safe
    append!(new._history.steps, domain._history.steps)
    new._history.obj_counter = domain._history.obj_counter
    _replay!(new; warn=false)
    return new
end

function change_resolution(domain::FITDomain, resolution; units="m")
    resolution isa Real && (resolution = (resolution, resolution, resolution))
    length(resolution) == 3 ||
        throw(ArgumentError("resolution needs 1 or 3 values, got $(length(resolution))"))
    unitToMeter = check_units(units)
    extents = (domain.nodes_u[end] - domain.nodes_u[1],
               domain.nodes_v[end] - domain.nodes_v[1],
               domain.nodes_w[end] - domain.nodes_w[1])
    edges = [fill(r * unitToMeter, _cell_count(L, r * unitToMeter)) for (L, r) in zip(extents, resolution)]
    return change_resolution(domain, edges...; units="m")
end

"""
    get_source(domain, name)   -> Vector{Float64}
    get_source(domain, obj_id) -> Vector{Float64}

Source vector of one source in `domain`, e.g. from `create_circular_loop_source`,
discretized on the current grid and scaled by its `current`. Length `3Np` in the
edge layout.

Use this after `change_resolution`, `remove_object!` or `undo!`: a vector returned
earlier by `create_circular_loop_source` is not updated by these.

Throws an `ArgumentError` if no such object exists or it is not a source.
"""
get_source(domain::FITDomain, name::AbstractString) = get_source(domain, _entry(domain, name).obj_id)

function get_source(domain::FITDomain, obj_id::Int)
    e = _entry(domain, obj_id)
    e.object isa AbstractSource ||
        throw(ArgumentError("object $obj_id \"$(e.name)\" is a $(nameof(typeof(e.object))), not a source"))
    return _discretize_source(domain, e.object)
end

# Result of list_objects: behaves like a Vector{CreateObject}, prints as a table
struct ObjectList <: AbstractVector{CreateObject}
    entries::Vector{CreateObject}
end

Base.size(l::ObjectList) = size(l.entries)
Base.getindex(l::ObjectList, i::Int) = l.entries[i]

"""
    list_objects(domain) -> ObjectList

All objects that currently exist in `domain`, in creation order; removed objects
are left out. Displays as a table. Index or iterate it like a vector; each entry
has the fields `obj_id`, `name` and `object`.

# Example

```julia
list_objects(domain)                     # show the table
remove_object!(domain, "coil")           # by name
J = get_source(domain, "coil")           # one source by name
```
"""
list_objects(domain::FITDomain) = ObjectList(_active_entries(domain._history))

function Base.show(io::IO, ::MIME"text/plain", l::ObjectList)
    n = length(l)
    print(io, n, n == 1 ? " object" : " objects")
    n == 0 && return
    println(io, ":")
    header = ("obj_id", "type", "name", "geometry (m)", "material")
    rows = [(string(e.obj_id), string(nameof(typeof(e.object))),
             isempty(e.name) ? "-" : e.name,
             _describe_geometry(e.object), _describe_material(e.object)) for e in l]
    w = [maximum(textwidth, (header[c], (r[c] for r in rows)...)) for c in 1:5]
    line(r) = join((rpad(r[c], w[c]) for c in 1:5), "  ")
    print(io, "  ", rstrip(line(header)))
    for r in rows
        print(io, "\n  ", rstrip(line(r)))
    end
end

_describe_material(m::Material) =
    "σ=$(m.σ) ε_r=$(m.ε_r) μ_r=$(m.μ_r)" * (m.color === nothing ? "" : " color=$(m.color)")

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