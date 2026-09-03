#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 06/08/2026

function get_divergence(config, ::Primal)
    Pu,Pv,Pw = GetPuPvPw(config)
    return [Pu Pv Pw]  
end

# S̃ = - transpose(G)
function get_divergence(config, ::Dual)
    G = get_gradient(config, Primal())
    return - G'
end

function S(config)
    return get_divergence(config, Primal())
end

function S̃(config)
    return get_divergence(config, Dual())
end