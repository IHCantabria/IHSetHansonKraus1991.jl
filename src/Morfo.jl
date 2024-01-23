function ADEAN(D50)
   ###########################################################################    
   # Dean parameter; D50 in meters   
   ###########################################################################    
    A=0.51.*wMOORE(D50).^0.44
    
    return A
end

function ALST(Hb,Dirb,hb,ANGbati,K)

    ###########################################################################    
    # Alongshore sediment transport
    #    
    # INPUT:
    # Hb:        wave height.
    # Tb:        wave period.
    # Dirb:      wave direction. Nautical convention.
    # hb:        depth of wave conditions.
    # ANGbati:   bathymetry angle; the normal of the shoreline in nautical degrees.
    # method:    alongshore sediment transport formula. Default is CERQ
    # calp:      calibration parameters
    # K1:        Komar calibration parameter
    # K:         SPM calibration parameter
    #    
    # OUTPUT:
    # q:      alongshore sediment transport relative to the bathymetry angle.
    #
    # DEPENDENCIAS:
    # rel_angle_cartesian; DIRrel
    ###########################################################################

    
    DIRrel = rel_angle_cartesian(nauticalDir2cartesianDir(Dirb),ANGbati)
    PerpRange = abs.(DIRrel) .< 90
    PerpRange = vec(PerpRange)
    q = zeros(size(Hb))
    q0 = zeros(size(Hb))
    # K1 = 0.39; # dimensionless, 0.39 SPM using Hs &  0.77 using Hrms Komar[1970]
    # K1 = calp.K1
    # q[PerpRange] = K[PerpRange].*(Hb[PerpRange].^2).*...
    #     cos(deg2rad[DIRrel[PerpRange]]).*sin(deg2rad[DIRrel[PerpRange]])
    rho = 1025 # #saltwater mass density SPM
    rhos = 2650 # # sand mass density SPM
    
    p = 0.4 # # porosity SPM
    gammab = Hb[PerpRange]./hb[PerpRange]
    gammab[isnan.(gammab)].=Inf
    cnts = rho.*sqrt.(9.81)./(16. .*sqrt.(gammab).*(rhos.-rho).*(1.0.-p));      
    q0[PerpRange] .= K[PerpRange].*cnts.*(Hb[PerpRange].^(5. /2.))
    q[PerpRange] .= q0[PerpRange].*sin.(2. .*deg2rad.(DIRrel[PerpRange]))
    #     elif method .== "spm":
    # #        K = 0.0098
    #         rho = 1025
    #         rhos = 2650
    #         lbda = 0.4
    # #        K = calp.K
    #         q[PerpRange] = rho/8.*K[PerpRange]/((rhos-rho)*lbda)*Hb[PerpRange]*np.sqrt(9.81*hb[PerpRange])*np.cos(np.radians[DIRrel[PerpRange]])*np.sin(np.radians[DIRrel[PerpRange]])
    #         
    # q[1] = q[2] + (q[3] - q[2])

    # if q[3]-q[2] < 0
    #     q[1] = q[2] + (q[3] - q[2])
    # else
    #     q[1] = q[2] - (q[3] - q[2])
    # end

    # q[end] = q[end-1] + (q[end-2] - q[end-1])
    # if q[end-2]-q[end-1] < 0
    #     q[end] = q[end-1] + (q[end-2] - q[end-1])
    # else
    #     q[end] = q[end-1] - (q[end-2] - q[end-1])
    # end


    q[1]=q[2]
    q[end]=q[end-1]
    # q[1]=0
    # q[end]=0
    

    return q, q0
end

function BruunRule(hc,D50,Hberm,slr)
    #    ###########################################################################    
    #    # Bruun Rule
    #    # INPUT:
    #    # hc:     depth of closure
    #    # D50:      Mean sediment grain size (m)
    #    # Hberm:    Berm Height [m]
    #    # slr:      Expected Sea Level Rise [m]
    #    # OUTPUT:
    #    # r:        expected progradation/recession [m]
    #
    #    ###########################################################################    
    Wc = wast(hc,D50)
    
    r=slr.*Wc./(Hberm.+hc)
    
    return r
end

function depthOfClosure(Hs12,Ts12)
    #    ###########################################################################    
    #    # Closure depth, Birkemeier[1985]
    #    # Hs12:     Significant wave height exceed 12 hours in a year.
    #    # Ts12:     Significant wave period exceed 12 hours in a year.
    #    # hs_ecdf = ECDF[hsDOWserie]
    #    # f_hs = interpolate.interp1d[hs_ecdf.y,hs_ecdf.x]
    #    # Hs12 = f_hs[1-12./365./24.]
    #    # Ts12 = 5.7*np.sqrt(Hs12)
    #    ###########################################################################
            
    dc = 1.75.*Hs12 - 57.9.*(Hs12.^2. /(9.81.*Ts12.^2))
        
    return dc
end

function Hs12Calc(Hs,Tp)

    ###########################################################################    
    # Significant Wave Height exceed 12 hours a year
    #
    # INPUT:
    # Hs:     Significant wave height.
    # Tp:     Wave Peak period.
    #
    # OUTPUT:
    # Hs12:     Significant wave height exceed 12 hours in a year.
    # Ts12:     Significant wave period exceed 12 hours in a year.
    ###########################################################################   
    
    Hs12=zeros(size(Hs))
    Ts12=zeros(size(Tp))
    # for i=1:size(Hs,2)
    for i in eachcol(Hs)
        Hs12calc=prctile[Hs[:,i],((365*24-12)/(365*24))*100]
        buscHS12=Hs[:,i].>=(Hs12calc-0.1) & Hs[:,i].<=(Hs12calc+0.1)
        f, xi = ksdensity[Tp[buscHS12,i]]
        _ , ii = maximum(f)
        Ts12[:,i]=xi[ii]
        Hs12[:,i]=Hs12calc
    end
    return Hs12,Ts12
end

function wast(hb,D50)

   ###########################################################################    
   # Width of the active surf zone
   # hb:   depth of closure
   # D50:  mean sediment grain size (m)
   ###########################################################################    
  #    hb = hb+CM ????? see why is introducing the tidal range in the width of the active surf zone    
    
    wsf=(hb./ADEAN(D50)).^(3.0/2.)
    
    return wsf
end

function wMOORE(D50)
    ###########################################################################    
    # Fall velocity; D50 in meters. Moore 1982
    ###########################################################################    
    
    ws = zeros(size(D50))
    for i in eachindex(D50)
        if D50[i].<=0.1*1e-3
            ws[i] = 1.1*1e6*D50[i].^2
        elseif D50[i].>0.1*1e-3 && D50[i].<1.0*1e-3
            ws[i] = 273.0.*D50[i].^1.1
        elseif D50[i].>1*1e-3
            ws[i] = 4.36.*D50[i].^0.5
        end
    end
    
    return ws
end