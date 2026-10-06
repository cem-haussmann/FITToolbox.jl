# SaveLoad.jl
# Norman Haussmann (haussmann@uni-wuppertal.de)
# Chair of Electromagnetic Theory, University of Wuppertal
# Date: 05/10/2026

"""
    save_domain(path, domain) -> path

Write the grid, the background material and the history of `domain` to the TOML
file `path`. `load_domain(path)` rebuilds the domain from it. `domain.material`
itself is not saved: values written directly into it are lost, with a warning.
"""
function save_domain(path::AbstractString, domain::FITDomain)
    _material_matches_history(domain) ||
        @warn "domain.material contains values that were not created through create_* \
               functions; they are not saved"
    h = domain._history
    d = Dict("format_version" => 1,
             "domain" => Dict("edges_u" => collect(domain.edges_u),
                              "edges_v" => collect(domain.edges_v),
                              "edges_w" => collect(domain.edges_w),
                              "sigma" => domain.σ, "eps_r" => domain.ε_r, "mu_r" => domain.μ_r,
                              "obj_counter" => h.obj_counter),
             "steps" => [_step_to_dict(s) for s in h.steps])
    open(io -> TOML.print(io, d; sorted=true), path, "w")
    return path
end

"""
    load_domain(path) -> FITDomain

Read a file written by `save_domain` and rebuild the domain: same grid, same
background material, same history, so `list_objects` and `undo!` work as before.

Throws an `ArgumentError` if the file is not a domain file or an object no longer
fits the grid.
"""
function load_domain(path::AbstractString)
    d = TOML.parsefile(path)
    get(d, "format_version", 0) == 1 ||
        throw(ArgumentError("$path is not a FITToolbox domain file of format version 1"))
    g = d["domain"]
    steps = [_step_from_dict(s) for s in get(d, "steps", [])]
    _check_steps(steps, g["obj_counter"], path)
    domain = create_domain(g["edges_u"], g["edges_v"], g["edges_w"];
                           σ=g["sigma"], ε_r=g["eps_r"], μ_r=g["mu_r"])
    append!(domain._history.steps, steps)
    domain._history.obj_counter = g["obj_counter"]
    _replay!(domain; warn=false)
    return domain
end

# the rules create_* and remove_object! enforce, for steps read from a file:
# unique obj_ids, unique names among existing objects, "obj<n>" only for obj_id n,
# deletes only of existing objects, and an obj_counter that hands out new ids
function _check_steps(steps, obj_counter, path)
    fail(msg) = throw(ArgumentError("$path: $msg"))
    seen   = Set{Int}()
    active = Dict{Int,String}()                       # obj_id => name
    for s in steps
        if s isa CreateObject
            s.obj_id in seen && fail("obj_id $(s.obj_id) is created twice")
            s.name in values(active) &&
                fail("an object named \"$(s.name)\" already exists (obj_id $(s.obj_id))")
            occursin(_DEFAULT_NAME, s.name) && s.name != "obj$(s.obj_id)" &&
                fail("name \"$(s.name)\" is reserved for obj_id $(s.name[4:end]), not $(s.obj_id)")
            push!(seen, s.obj_id)
            active[s.obj_id] = s.name
        else
            haskey(active, s.obj_id) ||
                fail("a step deletes obj_id $(s.obj_id), which does not exist at that point")
            delete!(active, s.obj_id)
        end
    end
    obj_counter >= maximum(seen; init=0) ||
        fail("obj_counter = $obj_counter is smaller than the largest obj_id $(maximum(seen))")
    return nothing
end

# one step as a flat Dict; a solid's material sits next to its geometry
_step_to_dict(s::DeleteObject) = Dict("op" => "delete", "obj_id" => s.obj_id)

function _step_to_dict(s::CreateObject)
    o = s.object
    d = Dict{String,Any}("op" => "create", "obj_id" => s.obj_id, "name" => s.name,
                         "type" => string(nameof(typeof(o))))
    if o isa Brick
        d["origin"]  = collect(o.origin)
        d["lengths"] = collect(o.lengths)
    elseif o isa Sphere
        d["center"] = collect(o.center)
        d["radius"] = o.radius
    elseif o isa CircularLoop
        d["center"]  = collect(o.center)
        d["radius"]  = o.radius
        d["normal"]  = string(nameof(typeof(o.normal)))
        d["current"] = real(o.current)
        o.current isa Complex && (d["current_im"] = imag(o.current))   # TOML has no complex numbers
    else
        throw(ArgumentError("cannot save objects of type $(nameof(typeof(o)))"))
    end
    if o isa AbstractSolid
        m = o.material
        d["sigma"], d["eps_r"], d["mu_r"] = m.σ, m.ε_r, m.μ_r
        m.color === nothing || (d["color"] = m.color)   # TOML has no null
    end
    return d
end

function _step_from_dict(d)
    op = get(d, "op", "")
    op == "delete" && return DeleteObject(d["obj_id"])
    op == "create" || throw(ArgumentError("unknown step op \"$op\" in domain file"))
    t = d["type"]
    material() = Material(d["sigma"], d["eps_r"], d["mu_r"], get(d, "color", nothing))
    o = if t == "Brick"
        Brick(Tuple(d["origin"]), Tuple(d["lengths"]), material())
    elseif t == "Sphere"
        Sphere(Tuple(d["center"]), d["radius"], material())
    elseif t == "CircularLoop"
        normal = Dict("DirX" => DirX(), "DirY" => DirY(), "DirZ" => DirZ())[d["normal"]]
        current = haskey(d, "current_im") ? complex(d["current"], d["current_im"]) : d["current"]
        CircularLoop(Tuple(d["center"]), d["radius"], normal, current)
    else
        throw(ArgumentError("unknown object type \"$t\" in domain file"))
    end
    return CreateObject(d["obj_id"], d["name"], o)
end
