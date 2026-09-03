#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 14/01/2026

function create_cube!(domain::FITDomain, u_o, v_o, w_o, u_length, v_length, w_length; units="m", σ=0.0, ε_r=1.0, μ_r=1.0)
    σ = Float64(σ)
    ε_r = Float64(ε_r)
    μ_r  = Float64(μ_r)
    unitToMeter = check_units(units)
    
    start_pos_u = u_o*unitToMeter
    start_pos_v = v_o*unitToMeter
    start_pos_w = w_o*unitToMeter

    plate_length_u = u_length*unitToMeter
    plate_length_v = v_length*unitToMeter
    plate_length_w = w_length*unitToMeter

    Nodes_U = domain.nodes_u
    Nodes_V = domain.nodes_v
    Nodes_W = domain.nodes_w

    u_min = start_pos_u 
    u_max = start_pos_u + plate_length_u
    v_min = start_pos_v 
    v_max = start_pos_v + plate_length_v
    w_min = start_pos_w 
    w_max = start_pos_w + plate_length_w

    # Check the requested extent against the grid before snapping. After
    # CheckClosest the indices are clamped into range, so testing them cannot
    # detect a cube that lies outside.
    (u_min >= Nodes_U[1] && u_max <= Nodes_U[end] &&
     v_min >= Nodes_V[1] && v_max <= Nodes_V[end] &&
     w_min >= Nodes_W[1] && w_max <= Nodes_W[end]) ||
        throw(ArgumentError("cube extends beyond the domain: requested u $(u_min)–$(u_max), \
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
        throw(ArgumentError("cube has zero extent in u after snapping to the grid — \
                             u_length = $u_length $units is smaller than one cell"))
    j_min < j_max ||
        throw(ArgumentError("cube has zero extent in v after snapping to the grid — \
                             v_length = $v_length $units is smaller than one cell"))
    k_min < k_max ||
        throw(ArgumentError("cube has zero extent in w after snapping to the grid — \
                             w_length = $w_length $units is smaller than one cell"))

    Material_distribution = domain.material

    @debug "Creating cube with snapped grid coordinates:" *
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