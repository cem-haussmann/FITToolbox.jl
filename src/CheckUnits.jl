#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 14/01/2026

function check_units(units)
    units == "nm" && return 1e-9
    (units == "μm" || units == "µm" || units == "um") && return 1e-6
    units == "mm" && return 1e-3
    units == "cm" && return 1e-2
    units == "dm" && return 1e-1
    units == "m"  && return 1.0
    units == "km" && return 1e3
    throw(ArgumentError("unknown unit $(repr(units))"))
end

function convert_to_meter(x, y, z, units)
    unitToMeter = check_units(units)
    return x*unitToMeter, y*unitToMeter, z*unitToMeter
end
