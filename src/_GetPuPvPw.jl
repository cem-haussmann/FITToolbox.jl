#This file was created by 
#Norman Haussmann (haussmann@uni-wuppertal.de), Chair of Electromagnetic Theory, University of Wuppertal
#Date: 16/09/2025

function GetPuPvPw(config)
    Nu = config.Nu
    Nv = config.Nv
    Nw = config.Nw
    
    Np = Nu*Nv*Nw 
    Mu = 1
    Mv = Nu
    Mw = Nu*Nv

    Pu = spdiagm(0 => -ones(Np), Mu => ones(Np-Mu))
    Pv = spdiagm(0 => -ones(Np), Mv => ones(Np-Mv))
    Pw = spdiagm(0 => -ones(Np), Mw => ones(Np-Mw))
    return Pu,Pv,Pw
end