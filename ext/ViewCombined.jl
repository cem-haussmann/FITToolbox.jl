# ViewCombined.jl
# The "combined" view of view_domain: the solids put together as smooth surfaces. Where
# they overlap the later one wins and cut-outs are subtracted, as in the material, but
# from the shapes as defined, not from the grid. Included from FITToolboxMakieExt.jl.
#
# Each solid is a signed distance field (negative inside); "A, minus a later B" is
# max(d_A, -d_B). The surface d = 0 is extracted with surface nets on a uniform sampling
# grid, independent of the FIT grid. Solids of the same colour and material that touch,
# e.g. the pieces of a wire, form a group: one surface, their union, with no seam
# between them, whose faces are shared out among them.

const _COMBINED_CELLS  = 240                # sampling cells along the longest side, at least
const _COMBINED_ACROSS = 3                  # samples across the thinnest part of a solid
const _COMBINED_BUDGET = 4_000_000          # samples around all solids together, at most

# The thinnest part of a solid, which the sampling has to resolve
_thinnest(b::FITToolbox.Brick)    = minimum(b.lengths)
_thinnest(s::FITToolbox.Sphere)   = 2 * s.radius
_thinnest(c::FITToolbox.Cylinder) = min(2 * c.radius, c.height)
_thinnest(t::FITToolbox.Torus)    = 2 * t.minor_radius
_thinnest(h::FITToolbox.Helix)    = 2 * h.minor_radius

# The sampling cells along the longest side of the domain: _COMBINED_CELLS, or more, so
# that the thinnest part of every solid, e.g. a thin plate, is _COMBINED_ACROSS samples
# thick, as far as the samples around all solids, which grow with the cube of the cells,
# stay within _COMBINED_BUDGET. The fields are sampled only around the solids, so a small
# solid in a large domain costs little.
function _combined_cells(domain)
    solids = [e.object for e in FITToolbox.list_objects(domain).entries
              if e.object isa FITToolbox.AbstractSolid && !_is_background(e.object)]
    isempty(solids) && return _COMBINED_CELLS
    L = maximum((domain.nodes_u[end] - domain.nodes_u[1], domain.nodes_v[end] - domain.nodes_v[1],
                 domain.nodes_w[end] - domain.nodes_w[1]))
    h = minimum(_thinnest, solids) / _COMBINED_ACROSS
    volume = sum(o -> ((lo, hi) = FITToolbox._extent(o); prod(hi .- lo)), solids)
    h = max(h, cbrt(volume / _COMBINED_BUDGET))
    return max(_COMBINED_CELLS, round(Int, L / h))
end

# The sampling grid: uniform nodes along each axis, over the domain
function _combined_grid(domain; cells = _combined_cells(domain))
    lo = (domain.nodes_u[1], domain.nodes_v[1], domain.nodes_w[1])
    hi = (domain.nodes_u[end], domain.nodes_v[end], domain.nodes_w[end])
    h = maximum(hi .- lo) / cells
    return ntuple(d -> collect(range(lo[d], hi[d]; length = max(1, round(Int, (hi[d] - lo[d]) / h)) + 1)), 3)
end

# The node indices of the grid within [lo, hi], and one node more on each side
_node_range(nodes, lo, hi) = ntuple(d -> max(1, searchsortedfirst(nodes[d], lo[d]) - 1):
                                         min(length(nodes[d]), searchsortedlast(nodes[d], hi[d]) + 1), 3)

# The field of solid `k` of `entries`: its signed distance minus every later solid that
# is not in its group `group` (obj_ids), on the nodes around it. Returns the node
# ranges and the values.
function _combined_field(nodes, entries, k, group)
    o = entries[k].object
    rng = _node_range(nodes, FITToolbox._extent(o)...)
    F = Array{Float64}(undef, length.(rng))
    fill_max!(F, o, sign, sub) = begin
        ref = FITToolbox._reference(o)
        Threads.@threads for c in collect(CartesianIndices(sub))
            i, j, l = Tuple(c)
            p = (nodes[1][rng[1][i]] - ref[1], nodes[2][rng[2][j]] - ref[2], nodes[3][rng[3][l]] - ref[3])
            v = sign * FITToolbox._signed_distance(o, p)
            F[i, j, l] = sign > 0 ? v : max(F[i, j, l], v)
        end
    end
    fill_max!(F, o, 1, axes(F))
    for later in entries[k+1:end]
        q = later.object
        q isa FITToolbox.AbstractSolid && !(later.obj_id in group) || continue
        lo, hi = FITToolbox._extent(q)
        sub = ntuple(d -> findall(i -> lo[d] - 2 * _step(nodes[d]) <= nodes[d][i] <=
                                       hi[d] + 2 * _step(nodes[d]), rng[d]), 3)
        any(isempty, sub) && continue
        fill_max!(F, q, -1, (first(sub[1]):last(sub[1]), first(sub[2]):last(sub[2]),
                             first(sub[3]):last(sub[3])))
    end
    return rng, F
end

# the step of uniform nodes
_step(v) = length(v) > 1 ? v[2] - v[1] : 1.0

# The surface F = 0 (inside F ≤ 0) by surface nets: a vertex in every cube the surface
# passes, at the mean of its crossings on the cube's edges, and a quad for every grid
# edge it crosses, joining the four cubes around that edge. Beyond the array F counts
# as outside, so the surface is closed. `xs`, `ys`, `zs` are the uniform node positions.
function _surface_nets(xs, ys, zs, F)
    n = size(F)
    h = (_step(xs), _step(ys), _step(zs))
    start = (xs[1], ys[1], zs[1])
    pad = minimum(h)
    val(i, j, k) = (1 <= i <= n[1] && 1 <= j <= n[2] && 1 <= k <= n[3]) ? F[i, j, k] : pad
    pos(i, j, k) = (start[1] + (i - 1) * h[1], start[2] + (j - 1) * h[2], start[3] + (k - 1) * h[3])
    corner(c) = ((c - 1) & 1, (c - 1) >> 1 & 1, (c - 1) >> 2 & 1)       # c = 1 + dx + 2dy + 4dz
    cube_edges = ((1, 2), (3, 4), (5, 6), (7, 8), (1, 3), (2, 4), (5, 7), (6, 8),
                  (1, 5), (2, 6), (3, 7), (4, 8))

    # vertices: cubes 0 … n along each axis, including the padding around the array
    vid = zeros(Int32, n[1] + 1, n[2] + 1, n[3] + 1)
    pts, nrm = Point3f[], Vec3f[]
    for k in 0:n[3], j in 0:n[2], i in 0:n[1]
        v = ntuple(c -> val(i + corner(c)[1], j + corner(c)[2], k + corner(c)[3]), 8)
        inside = count(<=(0), v)
        (inside == 0 || inside == 8) && continue
        acc, m = (0.0, 0.0, 0.0), 0
        for (a, b) in cube_edges
            (v[a] <= 0) == (v[b] <= 0) && continue
            t = v[a] / (v[a] - v[b])
            pa = pos((i, j, k) .+ corner(a)...)
            pb = pos((i, j, k) .+ corner(b)...)
            acc = acc .+ pa .+ t .* (pb .- pa)
            m += 1
        end
        push!(pts, Point3f(acc ./ m))
        # the normal: the gradient of F over the cube
        g = ((v[2] - v[1] + v[4] - v[3] + v[6] - v[5] + v[8] - v[7]) / h[1],
             (v[3] - v[1] + v[4] - v[2] + v[7] - v[5] + v[8] - v[6]) / h[2],
             (v[5] - v[1] + v[6] - v[2] + v[7] - v[3] + v[8] - v[4]) / h[3])
        push!(nrm, Vec3f(g ./ max(sqrt(sum(abs2, g)), eps())))
        vid[i + 1, j + 1, k + 1] = length(pts)
    end

    # faces: one quad per crossed grid edge, wound so that it faces outwards
    faces = _TriFace[]
    cube(i, j, k) = vid[i + 1, j + 1, k + 1]
    for d in 1:3
        e = ntuple(x -> Int(x == d), 3)
        # the two other axes, and the four cubes around an edge in order around it
        p, q = filter(!=(d), (1, 2, 3))
        for k in 0:n[3], j in 0:n[2], i in 0:n[1]
            node = (i, j, k)
            (node[p] >= 1 && node[q] >= 1) || continue
            node[d] <= n[d] || continue
            f0, f1 = val(node...), val((node .+ e)...)
            (f0 <= 0) == (f1 <= 0) && continue
            off(dp, dq) = ntuple(x -> x == p ? -dp : x == q ? -dq : 0, 3)
            c = (cube((node .+ off(1, 1))...), cube((node .+ off(0, 1))...),
                 cube((node .+ off(0, 0))...), cube((node .+ off(1, 0))...))
            any(iszero, c) && continue
            # The four cubes run counter-clockwise in the (p, q) plane, so the quad faces
            # along e_p × e_q, i.e. +d when (p, q, d) is in cyclic order. Outwards is +d
            # when the inside is at the lower node.
            out = f0 <= 0 ? 1 : -1
            v = out * (d == 2 ? -1 : 1) > 0 ? c : (c[1], c[4], c[3], c[2])
            # a quad is rarely flat: split it along the diagonal whose two triangles both
            # face along the surface normal, so that it does not fold over
            navg = nrm[v[1]] + nrm[v[2]] + nrm[v[3]] + nrm[v[4]]
            facing(a, b, cc) = (u = pts[b] - pts[a]; w = pts[cc] - pts[a];
                                sum((u[2]*w[3] - u[3]*w[2], u[3]*w[1] - u[1]*w[3],
                                     u[1]*w[2] - u[2]*w[1]) .* navg))
            if min(facing(v[1], v[2], v[3]), facing(v[1], v[3], v[4])) >=
               min(facing(v[1], v[2], v[4]), facing(v[2], v[3], v[4]))
                push!(faces, _TriFace(v[1], v[2], v[3]), _TriFace(v[1], v[3], v[4]))
            else
                push!(faces, _TriFace(v[1], v[2], v[4]), _TriFace(v[2], v[3], v[4]))
            end
        end
    end
    isempty(faces) && return nothing
    return Makie.GeometryBasics.Mesh(pts, faces; normal = nrm)
end

# The combined surfaces of all solids that are not cut-outs, by obj_id; with a cut, only
# what lies before `pos` along `axis`, closed at the cut. Each group is one surface, the
# union of its fields, shared out among its solids (see _share_faces).
function _combined_meshes(fields, groups, nodes; cut = nothing)
    meshes = Dict{Int,Any}()
    for group in groups
        ids = [id for id in group if haskey(fields, id)]
        isempty(ids) && continue
        # the nodes around all of them; beyond a solid's own nodes it is outside
        rng = ntuple(d -> (minimum(first(fields[id][1][d]) for id in ids):
                           maximum(last(fields[id][1][d]) for id in ids))::UnitRange{Int}, 3)
        G = fill(minimum(_step, nodes), length.(rng))
        owner = zeros(Int32, size(G))            # the solid, by its place in ids, at each node
        for (k, id) in enumerate(ids)
            r, F = fields[id]
            _merge_field!(G, owner, F, ntuple(d -> first(r[d]) - first(rng[d]), 3), k)
        end
        cut === nothing || _cut_field!(G, cut[1], nodes[cut[1]][rng[cut[1]]] .- cut[2])
        m = _surface_nets(nodes[1][rng[1]], nodes[2][rng[2]], nodes[3][rng[3]], G)
        m === nothing && continue
        merge!(meshes, length(ids) == 1 ? Dict(only(ids) => m) :
                                          _share_faces(m, ids, owner, rng, nodes))
    end
    return meshes
end

# G = min(G, F) on the nodes of F, which start at `offset` in G, noting solid k in `owner`
# where it is the smallest
function _merge_field!(G, owner, F, offset, k)
    for c in CartesianIndices(F)
        i = CartesianIndex(Tuple(c) .+ offset)
        F[c] < G[i] && (G[i] = F[c]; owner[i] = k)
    end
    return nothing
end

# Everything beyond the cut outside: G = max(G, x - pos) along axis a, with `beyond` the
# values x - pos at the nodes of G along it
function _cut_field!(G, a, beyond)
    for c in CartesianIndices(G)
        G[c] = max(G[c], beyond[c[a]])
    end
    return nothing
end

# The faces of a group's surface for each of its solids: a face belongs to the solid
# whose field is smallest at the sampling node nearest to its centre (`owner`, over the
# nodes `rng`), the one whose surface it lies on. The solids share the vertices along
# their borders, so together they are the group's surface, without a seam.
function _share_faces(m, ids, owner, rng, nodes)
    pts, nrm = m.position, m.normal
    nearest(p) = CartesianIndex(ntuple(d -> clamp(round(Int, (p[d] - nodes[d][1]) / _step(nodes[d])) + 2 -
                                                  first(rng[d]), 1, size(owner, d)), 3))
    owned = [_TriFace[] for _ in ids]
    for f in Makie.GeometryBasics.faces(m)
        push!(owned[max(1, owner[nearest((pts[f[1]] + pts[f[2]] + pts[f[3]]) / 3)])], f)
    end
    owned = Dict(zip(ids, owned))
    meshes = Dict{Int,Any}()
    for (id, fs) in owned
        isempty(fs) && continue
        # its own vertices, numbered anew
        used = sort!(unique!([convert(Int, i) for f in fs for i in f]))
        new = Dict(v => k for (k, v) in enumerate(used))
        renumber(f) = _TriFace((new[convert(Int, i)] for i in f)...)
        meshes[id] = Makie.GeometryBasics.Mesh(pts[used], renumber.(fs); normal = nrm[used])
    end
    return meshes
end

# The groups of solids drawn as one surface, as obj_ids in their order of creation:
# solids of the same colour and material whose extents touch or overlap, to within two
# sampling steps, and solids connected through such ones. Where they overlap it does
# not matter which one wins, so no seam is drawn between them. Every other solid that
# is not a cut-out forms a group of its own.
function _combined_groups(domain, nodes)
    entries = FITToolbox.list_objects(domain).entries
    colors = _object_colors(entries)
    solids = [e for e in entries if e.object isa FITToolbox.AbstractSolid && !_is_background(e.object)]
    key(e) = (m = e.object.material; (colors[e.obj_id], m.σ, m.ε_r, m.μ_r))
    margin = ntuple(d -> 2 * _step(nodes[d]), 3)
    function touch(a, b)
        (la, ha), (lb, hb) = FITToolbox._extent(a.object), FITToolbox._extent(b.object)
        return all(d -> la[d] - margin[d] <= hb[d] && lb[d] - margin[d] <= ha[d], 1:3)
    end
    # union-find over the solids, by their index in `solids`
    root = collect(eachindex(solids))
    find(k) = root[k] == k ? k : (root[k] = find(root[k]))
    for a in eachindex(solids), b in a+1:lastindex(solids)
        key(solids[a]) == key(solids[b]) && touch(solids[a], solids[b]) &&
            (root[find(b)] = find(a))
    end
    groups = Dict{Int,Vector{Int}}()
    for k in eachindex(solids)
        push!(get!(groups, find(k), Int[]), solids[k].obj_id)
    end
    return sort!(collect(values(groups)); by = first)
end

# The fields of all solids that are not cut-outs, by obj_id, each minus the later solids
# outside its group
function _combined_fields(domain, nodes; groups = _combined_groups(domain, nodes))
    group = Dict(id => Set(g) for g in groups for id in g)
    entries = [e for e in FITToolbox.list_objects(domain).entries if e.object isa FITToolbox.AbstractSolid]
    return Dict(entries[k].obj_id => _combined_field(nodes, entries, k, group[entries[k].obj_id])
                for k in eachindex(entries) if !_is_background(entries[k].object))
end
