#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

function get_gradient(config)
    @warn "This call is no longer supported, it is returning the gradient for the primal grid!"
    return get_gradient(config, Primal())
end

function get_gradient(config, ::Primal)
    Pu,Pv,Pw = GetPuPvPw(config)
    return [Pu;
            Pv;
            Pw]  
end

# G̃ = - transpose(S)
function get_gradient(config, ::Dual)
    S = get_divergence(config, Primal())
    return - S'
end

function G(config)
    return get_gradient(config, Primal())
end

function G̃(config)
    return get_gradient(config, Dual())
end