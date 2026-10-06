#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 14/01/2026

struct Brick <: AbstractSolid
    origin::NTuple{3,Float64}
    lengths::NTuple{3,Float64}
    material::Material
end

_describe_geometry(b::Brick) = "origin=$(_short(b.origin)) lengths=$(_short(b.lengths))"
_describe_material(b::Brick) = _describe_material(b.material)

function create_brick!(domain::FITDomain, u_o, v_o, w_o, u_length, v_length, w_length;
                       units="m", σ=0.0, ε_r=1.0, μ_r=1.0,
                       name::String="", color::Union{Nothing,String}=nothing)
    unitToMeter = check_units(units)
    brick = Brick((u_o*unitToMeter, v_o*unitToMeter, w_o*unitToMeter),
                  (u_length*unitToMeter, v_length*unitToMeter, w_length*unitToMeter),
                  Material(σ, ε_r, μ_r, color))
    name = _resolve_name(domain, name)
    _apply!(domain, brick)
    return _record_create!(domain, brick; name)
end


function _apply!(domain::FITDomain, brick::Brick)
    u_min, v_min, w_min = brick.origin
    u_max, v_max, w_max = brick.origin .+ brick.lengths
    σ, ε_r, μ_r = brick.material.σ, brick.material.ε_r, brick.material.μ_r

    Nodes_U = domain.nodes_u
    Nodes_V = domain.nodes_v
    Nodes_W = domain.nodes_w

    # Check the requested extent against the grid before snapping. After
    # CheckClosest the indices are clamped into range, so testing them cannot
    # detect a brick that lies outside.
    (u_min >= Nodes_U[1] && u_max <= Nodes_U[end] &&
     v_min >= Nodes_V[1] && v_max <= Nodes_V[end] &&
     w_min >= Nodes_W[1] && w_max <= Nodes_W[end]) ||
        throw(ArgumentError("brick extends beyond the domain: requested u $(u_min)–$(u_max), \
                             v $(v_min)–$(v_max), w $(w_min)–$(w_max) m"))
    
    function CheckClosest(index, nodes, pos)
        Nmax = length(nodes)
        index = min(index, Nmax)  # clamp if beyond last node
        if index > 1 && abs(nodes[index-1] - pos) < abs(nodes[index] - pos)
            index -= 1
        end
        return index
    end

    i_min = CheckClosest(searchsortedfirst(Nodes_U, u_min), Nodes_U, u_min)
    i_max = CheckClosest(searchsortedfirst(Nodes_U, u_max), Nodes_U, u_max)
    j_min = CheckClosest(searchsortedfirst(Nodes_V, v_min), Nodes_V, v_min)
    j_max = CheckClosest(searchsortedfirst(Nodes_V, v_max), Nodes_V, v_max)
    k_min = CheckClosest(searchsortedfirst(Nodes_W, w_min), Nodes_W, w_min)
    k_max = CheckClosest(searchsortedfirst(Nodes_W, w_max), Nodes_W, w_max)

    i_min < i_max ||
        throw(ArgumentError("brick has zero extent in u after snapping to the grid — \
                             u_length = $(brick.lengths[1]) m is smaller than one cell"))
    j_min < j_max ||
        throw(ArgumentError("brick has zero extent in v after snapping to the grid — \
                             v_length = $(brick.lengths[2]) m is smaller than one cell"))
    k_min < k_max ||
        throw(ArgumentError("brick has zero extent in w after snapping to the grid — \
                             w_length = $(brick.lengths[3]) m is smaller than one cell"))

    Material_distribution = domain.material

    @debug "Creating brick with snapped grid coordinates:" *
    "\n  u: $(Nodes_U[i_min]) m  to  $(Nodes_U[i_max]) m  (requested: $(u_min) m to $(u_max) m)" *
    "\n  v: $(Nodes_V[j_min]) m  to  $(Nodes_V[j_max]) m  (requested: $(v_min) m to $(v_max) m)" *
    "\n  w: $(Nodes_W[k_min]) m  to  $(Nodes_W[k_max]) m  (requested: $(w_min) m to $(w_max) m)"

    #set the properties to the domain
    md = Material_distribution
    @inbounds for i in i_min:i_max-1, j in j_min:j_max-1, k in k_min:k_max-1
        md[i,j,k,1] = σ
        md[i,j,k,2] = ε_r
        md[i,j,k,3] = μ_r
    end
end

@deprecate create_cube!(args...; kwargs...) create_brick!(args...; kwargs...)

