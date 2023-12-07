#%% Variables
    # ynew  : new shoreline position
    # y     : last shoreline position
    # r     : dune retreat
    # dev     : dune eroded volume
    # dt    : time step
    # dx    : distance between transect
    # ti    : index for the current step
    # wavec : wave climate at breaking during the new time step and 2 previous conditions
    #   hb  : wave height at breaking
    #   tb  : wave period at breaking
    #   θb: wave direction at breaking
    #   depthb  : breaking depth
    #   surge: storm surge
    #   tide: astronomical tide
    #   mmsl: monthly mean sea level NOT anomaly
    #   slr : mean sea level rise
    # q     : alongshore sediment transport from previous steps
    # yeq   : cross-shore equilibrium position from previous steps
    # ay0   : Miller&Dean maximum progadation position without forcings from previous steps
    # cnst  : model constants
    #   dc  : depth of closure DCT
    #   D50 : mean sediment grain size
    #   Hberm: berm height
    #   Hs12: wave height overcome 12hours in a year
    #   Ts12: representative period asociated to the wave height overcome 12 hours in a year
    # calp  : calibration parameters
    #   K1  : CERQ alongshore sediment transport calibration parameter kal
    #   AY0 : Miller&Dean maximum progadation position without forcings
    #   kero: Miller&Dean erosion rate kacr
    #   kacr: Miller&Dean accretion rate kero
    # trs   : transects position and angles
    #   X0  : x coordinate, origin of the transect
    #   Y0  : y coordinate, origin of the transect
    #   phi : transect orientation
    # bdc   : boundary conditions
    #   ql  : alongshore sediment transport sl
    #   qc  : cross shore sediment transport sc
    # slope : slope of the dry beach where swash
    # dunezz : dune toe and crest elevations
    # model : kind of processes to turn on ['crosshore','alongshore','mmsl','slr','dune']
    # theta : implicitness of the numerical scheme

    #    hs = wavec.hs
    #    tp = wavec.tp
    #    θe = wavec.dir
    #    depth = wavec.depth

    #    hb = wavec.hb
    #    tb = wavec.tb
    #    θb = wavec.dirb
    #    depthb = wavec.depthb

    #    surge = wavec.ss
    #    tide = wavec.tide
    #    mmsl = wavec.mmsl
    #    slr = wavec.msl




#

function run_OneLine()
    
    println("Loading libraries...")
    wrkDir = pwd()
    mods = wrkDir*"Modules\\"
    dats = wrkDir*"Data\\"

    # mods = wrkDir*"\\Modules\\"
    # dats = wrkDir*"\\Data\\"

    println("Loading datasets...")

    wavF = NCDataset(dats*"wav.nc")
    ensF = NCDataset(dats*"ens.nc")
    trsF = NCDataset(dats*"trs.nc")
    # slF = NCDataset(dats*"sl.nc")
    configF = NCDataset(dats*"config.nc")

    println("Unpacking datasets...")

    dt, dx, bctype, Ndata,  NS = configF["dt"][:][1], configF["dx"][:][1], configF["bctype"][1], configF["Ndata"][:], configF["NS"][:]
    
    kal = collect(skipmissing(parF["kal"][:]))

    Y0 = collect(skipmissing(trsF["Y0"][:]))
    X0 = collect(skipmissing(trsF["X0"][:]))
    Yf = collect(skipmissing(trsF["Yf"][:]))
    Xf = collect(skipmissing(trsF["Xf"][:]))
    phi = collect(skipmissing(trsF["phi"][:]))

    Hs = collect(skipmissing(wavF["hs"][:]))
    Tp = collect(skipmissing(wavF["tp"][:]))
    θ = collect(skipmissing(wavF["dir"][:]))
    depth = collect(skipmissing(wavF["depth"][:]))
    doc = collect(skipmissing(wavF["doc"][:]))

    yi = collect(skipmissing(configF["yi"][:]))

    n = length(Y0) + 1
    m = Int(length(Hs)/(n))

    close(wavF)
    close(ensF)
    close(configF)
    close(trsF)

    println("Datasets closed...")

    Hs = convert(Array{Float64},transpose(reshape(Hs, (m, n))))
    Tp = convert(Array{Float64},transpose(reshape(Tp, (m, n))))
    θ = convert(Array{Float64},transpose(reshape(θ, (m, n))))
    doc = convert(Array{Float64},transpose(reshape(doc, (m, n))))

    yi = convert(Array{Float64},yi)

    ##########START HERE#############

    println("Starting COCOONED - Longshore Only...")

    y_tot, q, hb, depthb, θb, q0 = OneLine(yi, dt, dx, Hs, Tp, θ, depth, 
                                            doc, kal, X0, Y0, phi, bctype)

    # cd("..")
    # cd("Results")

    println("\n\n****************Finished****************\n\n")
    # txt=['Tiempo para terminar:',num2str(floor(TT/60/60)),' h ',
    #     num2str(floor(60*(TT/60/60-floor(TT/60/60)))),' min ',
    #     num2str(60*(TT/60-floor(TT/60))),' s'];
    # disp(['Tiempo total de simulación: ', num2str((now-start_time)*24*60*60,'#.0f'), ' s']);

    #     # get the new positions

    XN,YN = abs_pos(X0,Y0,deg2rad.(phi),yi)
    return XN, YN, y_tot, q
end

############################################################################################
#################One-line#################One-line#################One-line#################
#################One-line#################One-line#################One-line#################
#################One-line#################One-line#################One-line#################
#################One-line#################One-line#################One-line#################
#################One-line#################One-line#################One-line#################
############################################################################################

function OneLine(yi, dt, dx, hs::Matrix{Float64}, tp::Matrix{Float64}, θ::Matrix{Float64}, depth, doc::Matrix{Float64}, kal, X0, Y0, phi, bctype)

    nti = size(hs,2)
    # nti = 1000
    tii = 2
    desl = 1

    n1 = length(X0)
    n2 = size(tp,1)
    mt = size(tp,2)

    ysol = convert(Array{Float64},zeros(n1, mt))
    ysol[:,1] .= yi

    # sedbgtal = zeros((n1, mt))
    # dQdx = zeros((n1, mt))
    hb = convert(Array{Float64},zeros(n2, mt))
    θb = convert(Array{Float64},zeros(n2, mt))
    depthb = convert(Array{Float64},zeros(n2, mt))
    q = convert(Array{Float64},zeros(n2, mt))
    q0 = convert(Array{Float64},zeros(n2, mt))
    
    for pos = 1:(nti-2)

        ti = pos+desl

        p1=ti-desl
        p2=ti+desl

        ynew, hb[:,ti],θb[:,ti],depthb[:,ti],
            q[:,ti], q0[:,ti] = ydir_L(ysol[:,ti-1],dt,dx,tii,hs[:,p1:p2],
            tp[:,p1:p2],θ[:,p1:p2],depth,hb[:,p1:p2],
            θb[:,p1:p2],depthb[:,p1:p2],q[:,p1:p2],
            doc[:,p1:p2],kal, X0, Y0, phi, bctype)

        # resFun(x) = residualL(x,ysol[:,ti-1],dt,dx,tii,hs[:,p1:p2],
        # tp[:,p1:p2],θ[:,p1:p2],depth,hb[:,p1:p2],
        # θb[:,p1:p2],depthb[:,p1:p2],q[:,p1:p2],
        # doc[:,p1:p2],kal, X0, Y0, phi,theta)

        # ynew = optimize.newton_krylov(resFun, ysol[:,ti-1] ; method="minres")

        # lb = ysol[:,ti-1] .- 1
        # ub = ysol[:,ti-1] .+ 1
        
        # ynew = nlboxsolve(resFun,ysol[:,ti-1],lb,ub, method= :jfnk, xtol=1e-1,ftol=1e-1).zero

        # ynew[isnan.(ynew)] .= ysol[isnan.(ynew),ti-1]
        ysol[:,ti] = ynew

        p1=p1+1
        p2=p2+1

        # hb[:,ti],θb[:,ti],depthb[:,ti], q[:,ti], sedbgtal[:,ti],
        # dQdx[:,ti]= residualL_vol(ysol[:,ti],ysol[:,ti-1],dt,
        # dx,tii,hs[:,p1:p2],tp[:,p1:p2],θ[:,p1:p2],depth,hb[:,p1:p2],θb[:,p1:p2],
        # depthb[:,p1:p2],q[:,p1:p2],doc[:,p1:p2],kal, X0, Y0, phi,theta)

        if pos % 100 == 0
            @printf("\n Progress of %.2f %%",pos/(nti-1) .* 100)
        end
        # println(pos)
    #     #    res[pos+1,:], WAVEC.hb[pos+1,:], WAVEC.θb[pos+1,:], WAVEC.depthb[pos+1,:], q[pos+1,:], yeq[pos+1,:], CALP.AY0[pos+1,:], sedbgt[pos+1,:], sedbgtal[pos+1,:], sedbgtcr[pos+1,:], sedbgtbc[pos+1,:] = residual(ynew,y,dt,dx,pos+1,WAVEC,q,yeq,CNST,CALP,TRS,BDC,model,theta)
    #     #    yold = y    
        # tiempo_iter(pos)=toc
    end

    return ysol, q, hb, depthb, θb, q0
end

function ydir_L(y,dt,dx,ti,hs,tp,θe,depth,hb,θb,depthb,q,doc,kal,X0, Y0, phi, bctype)

    XN, YN = abs_pos(X0,Y0,nauticalDir2cartesianDir(phi).*pi./180.,y)
    
    alfas = zeros(size(hs,1))
    alfas[2:end-1] = shore_angle(XN,YN,θe[:,ti])
    # alfas = GM.cartesianDir2nauticalDir(alfas)
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    # println(alfas)
    try
        hb[:,ti], θb[:,ti], depthb[:,ti] = BreakingPropagation(hs[:,ti],tp[:,ti],θe[:,ti],depth,alfas .+ 90, "spectral")
    catch
        # println("\n")
        # println(alfas .- angRel)
        # println("\n")
        # println(hs[:,ti])
        # println("\n")
        # println(tp[:,ti])
        # println("\n")
        # println(dire[:,ti])
        # println("\n") 
        println("\n Waves diverged -- Q_tot = 0")     
    end
    # println(hb[1,ti])
    #modificado LFP
    # 
    # kal = kal[:,ti]
    
    dc = 0.5.*(doc[2:end,ti-1:ti+1]+doc[1:end-1,ti-1:ti+1]); # # trapezoidal rule for the closure depth in m+1/2

    q[:,ti], q0 = ALST(hb[:,ti], θb[:,ti], depthb[:,ti], alfas.+ 90, kal) # # for Tairua is -90 according to how the reference system was defined
    if bctype == "Dirichlet"
        q[1,ti] = 0
        q[end,ti] = 0
    elseif bctype == "Neumann"
        q[1,ti] = q[2,ti]
        q[end,ti] = q[end-1,ti]
    end

    if (dx^2  * minimum(dc) / (4*maximum(q0))) < dt
        println("WARNING: COURANT CONDITION VIOLATED")
    end

    # println("qs=", size(q))
    # println("ysol=", size(y))
    # println("kal=", size(kal))
    return y .- (dt*60*60) ./ dc[:,ti] .* (q[2:end,ti].-q[1:end-1,ti])./dx, hb[:,ti], θb[:,ti], depthb[:,ti], q[:,ti], q0

end

function residualL(ynew,y,dt,dx,ti,hs,tp,θe,depth,hb,θb,depthb,q,doc,kal,X0, Y0, phi,theta, bctype, algulo_rel)

    # println(sum(isreal.(ynew)))
    XN, YN = abs_pos(X0,Y0,nauticalDir2cartesianDir(phi).*pi./180.,ynew)
    
    alfas = zeros(size(hs,1))
    alfas[2:end-1] = shore_angle(XN,YN,θe[:,ti])
    alfas[1] = alfas[2]; alfas[end] = alfas[end-1]; # # ghost condition for the relative angle()
    
    hb[:,ti], θb[:,ti], depthb[:,ti] .= BreakingPropagation(hs[:,ti],tp[:,ti],θe[:,ti],depth,alfas .- angulo_rel, "spectral")

    # println(hb[1,ti])
    #modificado LFP
    # 
    # kal = kal[:,ti]
    
    dc = 0.5.*(doc[2:end,ti-1:ti+1]+doc[1:end-1,ti-1:ti+1]); # # trapezoidal rule for the closure depth in m+1/2

    q[:,ti], _ = ALST(hb[:,ti],θb[:,ti],depthb[:,ti],alfas .- angulo_rel,kal) # # for Tairua is -90 according to how the reference system was defined

    if bctype == "Dirichlet"
        q[1,ti] = 0
        q[end,ti] = 0
    elseif bctype == "Neumann"
        q[1,ti] = q[2,ti]
        q[end,ti] = q[end-1,ti]
    end

    return (ynew.-y)./dt .+ (theta.*(q[2:end,ti].-q[1:end-1,ti])./dx./dc[:,ti] .+ (1-theta).*(q[2:end,ti-1].-q[1:end-1,ti-1])./dx./dc[:,ti-1]) 

end
