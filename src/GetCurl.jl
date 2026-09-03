#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

function get_curl(config)
    @warn "This call is no longer supported, it is returning the curl operator for the primal grid!"
    return get_curl(config, Primal())
end


function get_curl(config, ::Primal)
    Nu = config.Nu
    Nv = config.Nv
    Nw = config.Nw
    Np = Nu*Nv*Nw
 
    Pu,Pv,Pw = GetPuPvPw(config)
    C = [spzeros(Np,Np) -1*Pw Pv;
            Pw spzeros(Np,Np) -1*Pu;
            -1*Pv Pu spzeros(Np,Np)] 
    return C
end

#C̃ = transpose(C)
function get_curl(config, ::Dual)
    C = get_curl(config, Primal())
    return sparse(C')
end

function C(config)
    return get_curl(config, Primal())
end

function C̃(config)
    return get_curl(config, Dual())
end