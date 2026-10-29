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

Throws an `ArgumentError` that names the file and the place in it if the file is not
a domain file of format version 1, if a value is missing or of the wrong kind (e.g. a
position without three coordinates, an unknown axis or colour, a cell width that is
not positive), if the history is one that `create_*` and `remove_object!` could not
have produced (an obj_id or name used twice, an empty name, a delete of an object
that does not exist), or if an object no longer fits the grid. A file that is not
TOML at all throws the `TOML.ParserError` of the parser.
"""
function load_domain(path::AbstractString)
    d = TOML.parsefile(path)
    get(d, "format_version", 0) == 1 ||
        throw(ArgumentError("$path is not a FITToolbox domain file of format version 1"))
    g = _value(d, "domain", path, x -> x isa AbstractDict, "a table, [domain]")
    place = "$path, [domain]"
    edges = [_value(g, key, place, x -> x isa AbstractVector && !isempty(x) &&
                                        all(y -> _is_number(y) && isfinite(y) && y > 0, x),
                    "a list of positive cell widths") for key in ("edges_u", "edges_v", "edges_w")]
    σ, ε_r, μ_r = _number(g, "sigma", place), _number(g, "eps_r", place), _number(g, "mu_r", place)
    obj_counter = _value(g, "obj_counter", place, x -> _is_integer(x) && x >= 0,
                         "a non-negative integer")
    raw = get(d, "steps", [])
    raw isa AbstractVector ||
        throw(ArgumentError("$path: \"steps\" must be an array of tables, [[steps]]"))
    steps = [_step_from_dict(s, "$path, step $n") for (n, s) in enumerate(raw)]
    _check_steps(steps, obj_counter, path)
    domain = create_domain(edges...; σ, ε_r, μ_r)
    append!(domain._history.steps, steps)
    domain._history.obj_counter = obj_counter
    _replay!(domain; warn=false)
    _check_sources(domain)
    return domain
end

# Reading values from a file. TOML has integers, floats, booleans and strings, and
# Bool <: Integer in Julia, so `true` must not pass for a number.
_is_number(x)  = x isa Real && !(x isa Bool)
_is_integer(x) = x isa Integer && !(x isa Bool)

# The value of `key` in the table `t`, if it is `ok`; else an ArgumentError that names
# the `place` in the file and what was `expected`
function _value(t, key, place, ok, expected)
    haskey(t, key) || throw(ArgumentError("$place: \"$key\" is missing"))
    v = t[key]
    ok(v) || throw(ArgumentError("$place: \"$key\" must be $expected, got $(_shown(v))"))
    return v
end
_shown(v) = (s = repr(v); length(s) > 60 ? first(s, 57) * "..." : s)

_number(t, key, place) = _value(t, key, place, _is_number, "a number")
_numbers(t, key, place, n) =
    Tuple(_value(t, key, place, x -> x isa AbstractVector && length(x) == n && all(_is_number, x),
                 "$n numbers"))
_flag(t, key, place) = _value(t, key, place, x -> x isa Bool, "true or false")
_text(t, key, place) = _value(t, key, place, x -> x isa AbstractString, "a string")

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
            isempty(s.name) && fail("obj_id $(s.obj_id) has an empty name")
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
    elseif o isa Cylinder
        d["base"]   = collect(o.base)
        d["radius"] = o.radius
        d["height"] = o.height
        d["axis"]   = string(nameof(typeof(o.axis)))
    elseif o isa Torus
        d["center"]       = collect(o.center)
        d["major_radius"] = o.major_radius
        d["minor_radius"] = o.minor_radius
        d["axis"]         = string(nameof(typeof(o.axis)))
    elseif o isa Helix
        d["base"]         = collect(o.base)
        d["major_radius"] = o.major_radius
        d["minor_radius"] = o.minor_radius
        d["pitch"]        = o.pitch
        d["turns"]        = o.turns
        d["axis"]         = string(nameof(typeof(o.axis)))
        d["start_angle"]  = o.start_angle
        d["left_handed"]  = o.left_handed
    elseif o isa PointCharge
        d["position"] = collect(o.position)
        d["charge"]   = o.charge
    elseif o isa PolylineSource
        d["points"]  = [collect(p) for p in o.points]
        d["closed"]  = o.closed
        d["current"] = real(o.current)
        o.current isa Complex && (d["current_im"] = imag(o.current))   # TOML has no complex numbers
    elseif o isa RectangularLoop
        d["origin"]  = collect(o.origin)
        d["lengths"] = collect(o.lengths)
        d["normal"]  = string(nameof(typeof(o.normal)))
        d["current"] = real(o.current)
        o.current isa Complex && (d["current_im"] = imag(o.current))   # TOML has no complex numbers
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

const _DIRECTIONS = Dict("DirX" => DirX(), "DirY" => DirY(), "DirZ" => DirZ())

# One step from a file, every value checked; `place` names it in the errors
function _step_from_dict(d, place)
    d isa AbstractDict || throw(ArgumentError("$place is not a table"))
    op = get(d, "op", "")
    obj_id() = _value(d, "obj_id", place, x -> _is_integer(x) && x >= 1, "a positive integer")
    op == "delete" && return DeleteObject(obj_id())
    op == "create" ||
        throw(ArgumentError("$place: unknown op $(_shown(op)), neither \"create\" nor \"delete\""))
    t = _text(d, "type", place)
    num(key)   = _number(d, key, place)
    point(key) = _numbers(d, key, place, 3)
    axis(key)  = _DIRECTIONS[_value(d, key, place, x -> x isa AbstractString && haskey(_DIRECTIONS, x),
                                    "\"DirX\", \"DirY\" or \"DirZ\"")]
    function material()
        color = haskey(d, "color") ?
            _value(d, "color", place, c -> c isa AbstractString && _valid_color(c),
                   "a hex code like \"#B87333\" or a colour name") : nothing
        return Material(num("sigma"), num("eps_r"), num("mu_r"), color)
    end
    current() = haskey(d, "current_im") ? complex(num("current"), num("current_im")) : num("current")
    o = if t == "Brick"
        Brick(point("origin"), point("lengths"), material())
    elseif t == "Sphere"
        Sphere(point("center"), num("radius"), material())
    elseif t == "Cylinder"
        Cylinder(point("base"), num("radius"), num("height"), axis("axis"), material())
    elseif t == "Torus"
        Torus(point("center"), num("major_radius"), num("minor_radius"), axis("axis"), material())
    elseif t == "Helix"
        Helix(point("base"), num("major_radius"), num("minor_radius"), num("pitch"), num("turns"),
              axis("axis"), num("start_angle"), _flag(d, "left_handed", place), material())
    elseif t == "PointCharge"
        PointCharge(point("position"), num("charge"))
    elseif t == "PolylineSource"
        points = _value(d, "points", place,
                        x -> x isa AbstractVector &&
                             all(p -> p isa AbstractVector && length(p) == 3 && all(_is_number, p), x),
                        "a list of points of 3 numbers each")
        PolylineSource([Tuple(Float64.(p)) for p in points], _flag(d, "closed", place), current())
    elseif t == "RectangularLoop"
        RectangularLoop(point("origin"), _numbers(d, "lengths", place, 2), axis("normal"), current())
    elseif t == "CircularLoop"
        CircularLoop(point("center"), num("radius"), axis("normal"), current())
    else
        throw(ArgumentError("$place: unknown object type $(_shown(t))"))
    end
    return CreateObject(obj_id(), _text(d, "name", place), o)
end
