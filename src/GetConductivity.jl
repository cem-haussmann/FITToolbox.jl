#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

function get_conductivity(config)
    return _get_dual_facet_averaged_material(config,1)
end

function M_σ(config)
    return get_conductivity(config)
end