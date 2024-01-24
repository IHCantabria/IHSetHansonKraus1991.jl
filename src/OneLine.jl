"""
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
"""

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
    hb = convert(Array{Float64}, zeros(n2, mt))
    θb = convert(Array{Float64}, zeros(n2, mt))
    depthb = convert(Array{Float64}, zeros(n2, mt))
    q = convert(Array{Float64}, zeros(n2, mt))
    q0 = convert(Array{Float64}, zeros(n2, mt))
    

    time_init = now()
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

        ysol[:,ti] = ynew

        p1=p1+1
        p2=p2+1


        if pos % 1000 == 0
            elp_t = (now() - time_init).value
            @printf("\n Progress of %.2f %% - ",pos/(nti-1) .* 100)
            @printf("Average time per step: %.2f [ms] - ",(elp_t/pos))
            @printf("Estimated time to finish: %.2f [s] - ",(elp_t/1000/pos)*(nti-2-pos))
            @printf("Elapsed time: %.2f [s]", elp_t/1000)
        end
    end
    
    println("\n***************************************************************")
    println("End of simulation")
    println("***************************************************************")
    @printf("\nElapsed simulation time: %.2f seconds \n",(now() - time_init).value/1000)
    println("***************************************************************")
    
    return ysol, q
end

function ydir_L(y,dt,dx,ti,hs,tp,θe,depth,hb,θb,depthb,q,doc,kal,X0, Y0, phi, bctype)

    XN, YN = abs_pos(X0,Y0,nauticalDir2cartesianDir(phi).*pi./180.,y)
    
    alfas = zeros(size(hs,1))
    alfas[2:end-1] = shore_angle(XN,YN,θe[:,ti])
    # alfas = cartesianDir2nauticalDir(alfas)
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

    q[:,ti], q0 = ALST(hb[:,ti], tp[:,ti],θb[:,ti], depthb[:,ti], alfas.+ 90, kal) # # for Tairua is -90 according to how the reference system was defined
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

function run_OneLine()
    
    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data_1L/"

    println("Loading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"
    trsF = dats*"par.nc"
    parF = dats*"par.nc"
    # slF = NCDataset(dats*"sl.nc")
    configF = dats*"config.nc"

    println("Unpacking datasets...")

    dt, dx, bctype = ncread(configF, "dt"), ncread(configF, "dx"), ncread(configF, "bctype")

    kal = ncread(parF, "kal")
    Y0 = ncread(trsF, "Y0")
    X0 = ncread(trsF, "X0")
    Yf = ncread(trsF, "Yf")
    Xf = ncread(trsF, "Xf")
    phi = ncread(trsF, "phi")

    Hs = ncread(wavF, "hs")
    Tp = ncread(wavF, "tp")
    θ = ncread(wavF, "dir")
    depth = ncread(wavF, "depth")
    doc = ncread(wavF, "doc")

    yi = ncread(configF, "yi")

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    time_w = DateTime.(YY, MM, DD, HH)

    ##########START HERE#############

    println("Starting COCOONED - Longshore Only...")

    y_tot, q = OneLine(yi, dt, dx, Hs, Tp, θ, depth, doc, kal, X0, Y0, phi, bctype)

    # cd("..")
    # cd("Results")

    println("\n\n****************Finished****************\n\n")
    # txt=['Tiempo para terminar:',num2str(floor(TT/60/60)),' h ',
    #     num2str(floor(60*(TT/60/60-floor(TT/60/60)))),' min ',
    #     num2str(60*(TT/60-floor(TT/60))),' s'];
    # disp(['Tiempo total de simulación: ', num2str((now-start_time)*24*60*60,'#.0f'), ' s']);

    #     # get the new positions

    XN,YN = abs_pos(X0, Y0, deg2rad.(phi), yi)
    
    return XN, YN, y_tot, q

end

function cal_OneLine()
    
    println("Loading libraries...")
    wrkDir = pwd()
    dats = wrkDir*"/data_1L/"

    println("Loading datasets...")

    wavF = dats*"wav.nc"
    ensF = dats*"ens.nc"
    trsF = dats*"trs.nc"
    # slF = NCDataset(dats*"sl.nc")
    configF = dats*"config.nc"

    println("Unpacking datasets...")

    dt, dx, bctype, MetObj = ncread(configF, "dt")[1], ncread(configF, "dx")[1], ncread(configF, "bctype")[1], ncread(configF, "MetObj")[1]

    Y0 = ncread(trsF, "Y0")
    X0 = ncread(trsF, "X0")
    Yf = ncread(trsF, "Yf")
    Xf = ncread(trsF, "Xf")
    phi = ncread(trsF, "phi")

    Hs = ncread(wavF, "hs")
    Tp = ncread(wavF, "tp")
    θ = ncread(wavF, "dir")
    depth = ncread(wavF, "depth")
    doc = ncread(wavF, "doc")

    yi = ncread(configF, "yi")

    YY, MM, DD, HH = ncread(wavF, "Y"), ncread(wavF, "M"), ncread(wavF, "D"), ncread(wavF, "h")

    YYo, MMo, DDo, HHo = ncread(ensF, "Y"), ncread(ensF, "M"), ncread(ensF, "D"), ncread(ensF, "h")
    
    Y_obs = ncread(ensF, "Obs")

    t_obs = DateTime.(YYo, MMo, DDo, HHo)

    t_wav = DateTime.(YY, MM, DD, HH)

    ii =  t_obs .<= t_wav[end] .&& t_obs .>= t_wav[1]

    t_obs, Y_obs = t_obs[ii], Y_obs[:,ii]

    ii =  t_wav .<= t_obs[end] .&& t_wav .>= t_obs[1]

    t_wav, Hs, Tp, θ = t_wav[ii], Hs[:,ii], Tp[:,ii], θ[:,ii]

    idx_obs = zeros(length(t_obs))

    for i in eachindex(t_obs)
        idx_obs[i] = argmin(abs.(t_wav .- t_obs[i]))
    end

    idx_obs = convert(Array{Int64},idx_obs)    

    ##########START HERE#############

    time_start = now()

    println("Starting COCOONED - Longshore Only...")

    function Calibra_(Χ)
        # kal = fill(exp(Χ[1]), size(Hs,1))
        kal = fill(Χ[1], size(Hs,1))
        Ymd, _ = OneLine(yi, dt, dx, Hs, Tp, θ, depth, doc, kal, X0, Y0, phi, bctype)
        YYsl = Ymd[:,idx_obs]
        if MetObj == "Pearson"
            rp = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                rp[i] = 1 -  abs(sum((YYsl[:,i].-mean(YYsl[:,i])).*(Y_obs[:,i] .- mean(Y_obs[:,i])))/(std(YYsl[:,i])*std(Y_obs[:,i])*length(YYsl[:,i])))
            end
            return mean(rp)
        elseif MetObj == "RMSE"
            rmse = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                rmse[i] = abs(sqrt(mean((YYsl[:,i] .- Y_obs[:,i]).^2))/5)
            end
            return mean(rmse)
        elseif MetObj == "MSS"
            mss = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                mss[i] = sum((YYsl[:,i] .- Y_obs[:,i]).^2)/length(YYsl[:,i])/(var(YYsl[:,i])+var(Y_obs[:,i])+(mean(YYsl[:,i])-mean(Y_obs[:,i]))^2)
            end
            return mean(mss)
        elseif MetObj == "BSS"
            bss = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                bss[i] = (mean((YYsl[:,i] .- Y_obs[:,i]).^2) - mean((YYref[:,i] .- Y_obs[:,i]).^2))/mean((YYref[:,i] .- Y_obs[:,i]).^2)
            end
            return mean(bss)
        elseif MetObj == "Double"
            mss = zeros(size(YYsl,2))
            rmse = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                mss[i] = sum((YYsl[:,i] .- Y_obs[:,i]).^2)/length(YYsl[:,i])/(var(YYsl[:,i])+var(Y_obs[:,i])+(mean(YYsl[:,i])-mean(Y_obs[:,i]))^2)
                rmse[i] = abs(sqrt(mean((YYsl[:,i] .- Y_obs[:,i]).^2))/5)
            end
            return (mean(mss), mean(rmse))
        elseif MetObj == "Triple"
            mss = zeros(size(YYsl,2))
            rmse = zeros(size(YYsl,2))
            rp = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                mss[i] = sum((YYsl[:,i] .- Y_obs[:,i]).^2)/length(YYsl[:,i])/(var(YYsl[:,i])+var(Y_obs[:,i])+(mean(YYsl[:,i])-mean(Y_obs[:,i]))^2)
                rmse[i] = abs(sqrt(mean((YYsl[:,i] .- Y_obs[:,i]).^2))/5)
                rp[i] = 1 -  abs(sum((YYsl[:,i].-mean(YYsl[:,i])).*(Y_obs[:,i] .- mean(Y_obs[:,i])))/(std(YYsl[:,i])*std(Y_obs[:,i])*length(YYsl[:,i])))
            end
            return (mean(mss), mean(rmse), mean(rp))
        elseif MetObj == "Double2"
            mss = zeros(size(YYsl,2))
            rp = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                mss[i] = sum((YYsl[:,i] .- Y_obs[:,i]).^2)/length(YYsl[:,i])/(var(YYsl[:,i])+var(Y_obs[:,i])+(mean(YYsl[:,i])-mean(Y_obs[:,i]))^2)
                rp[i] = 1 -  abs(sum((YYsl[:,i].-mean(YYsl[:,i])).*(Y_obs[:,i] .- mean(Y_obs[:,i])))/(std(YYsl[:,i])*std(Y_obs[:,i])*length(YYsl[:,i])))
            end
            return (mean(mss), mean(rp))
        elseif MetObj == "Double3"
            rmse = zeros(size(YYsl,2))
            rp = zeros(size(YYsl,2))
            for i in keys(YYsl[1,:])
                rmse[i] = abs(sqrt(mean((YYsl[:,i] .- Y_obs[:,i]).^2))/5)
                rp[i] = 1 -  abs(sum((YYsl[:,i].-mean(YYsl[:,i])).*(Y_obs[:,i] .- mean(Y_obs[:,i])))/(std(YYsl[:,i])*std(Y_obs[:,i])*length(YYsl[:,i])))
            end
            return (mean(rmse), mean(rp))
        end
    end

    PS = 10
    ME = 10

    boundsr = [(1e-2, 5)]

    if MetObj == "Double" || MetObj == "Double2" || MetObj == "Double3"
        resr = bboptimize(Calibra_; 
                        # Method = :simultaneous_perturbation_stochastic_approximation,
                        SearchRange = boundsr,
                        NumDimensions = 1,
                        PopulationSize = PS,
                        MaxSteps = ME,
                        FitnessTolerance = 1e-6,
                        FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                        TraceMode=:compact,
                        ϵ=0.1,
                        τ = 0.05,
                        MaxStepsWithoutEpsProgress = 100,
                        Method=:borg_moea)
    elseif MetObj == "Triple"
        resr = bboptimize(Calibra_; 
                        # Method = :simultaneous_perturbation_stochastic_approximation,
                        SearchRange = boundsr,
                        NumDimensions = 1,
                        PopulationSize = PS,
                        MaxSteps = ME,
                        FitnessTolerance = 1e-6,
                        FitnessScheme=ParetoFitnessScheme{3}(is_minimizing=true),
                        TraceMode=:compact,
                        ϵ=0.1,
                        τ = 0.05,
                        MaxStepsWithoutEpsProgress = 100,
                        Method=:borg_moea)
    else
        resr = bboptimize(Calibra_; 
                        Method = :adaptive_de_rand_1_bin,
                        SearchRange = boundsr,
                        NumDimensions = 1,
                        PopulationSize = PS,
                        MaxSteps = ME,
                        FitnessTolerance = 1e-6,
                        TraceMode=:compact,
                        ϵ=0.1,
                        τ = 0.05,
                        MaxStepsWithoutEpsProgress = 100)
    end

    popr = best_candidate(resr)

    Ymdr, q_tot = OneLine(yi, dt, dx, Hs, Tp, θ, depth, doc, exp(popr[1]), X0, Y0, phi, bctype)

    Ysl = Ymdr[idx_obs,:]
    aRP = sum((Ysl.-mean(Ysl)).*(Y_obs .- mean(Y_obs)))/(std(Ysl)*std(Y_obs)*length(Ysl))
    aRMSE = sqrt(mean((Ysl .- Y_obs).^2))
    aMSS = 1 - sum((Ysl .- Y_obs).^2)/length(Ysl)/(var(Ysl)+var(Y_obs)+(mean(Ysl)-mean(Y_obs))^2)

    XN,YN = abs_pos(X0, Y0, deg2rad.(phi), yi)

    total_time = (now() - time_init).value/1000
    minutes = floor(total_time/60)
    seconds = total_time - minutes*60

    println("\n****************************************************************")
    println("End of Calibration")
    println("***************************************************************")
    @printf("\nElapsed calibration time: %.0f min and %.2f\n", minutes, seconds)
    println("****************************************************************")

    println("\n\n****************Writing output****************\n\n")

    year_atts = Dict("long_name" => "Year")
    month_atts = Dict("long_name" => "Month")
    day_atts = Dict("long_name" => "Day")
    hour_atts = Dict("long_name" => "Hour")
    trs_atts = Dict("long_name" => "Transect number")

    output = wrkDir*"/results/Shoreline_HansonKraus1989.nc"
    nccreate(output, "year",
                "dim", length(YY),
                atts = year_atts)
    ncwrite(YY, output, "year")
    nccreate(output, "month",
                "dim", length(MM),
                atts = month_atts)
    ncwrite(MM, output, "month")
    nccreate(output, "day",
                "dim", length(DD),
                atts = day_atts)
    ncwrite(DD, output, "day")
    nccreate(output, "hour",
                "dim", length(HH),
                atts = hour_atts)
    ncwrite(HH, output, "hour")  
    nccreate(output, "trs",
                "dim", length(X0),
                atts = trs_atts)
    ncwrite(1:length(X0), output, "trs")

    Y_atts = Dict("units" => "m",
        "long_name" => "Shoreline position",
        "standard_name" => "Y")
    Yi_atts = Dict("units" => "m",
        "long_name" => "Initial shoreline position",
        "standard_name" => "Yi")
    XN_atts = Dict("units" => "m",
        "long_name" => "X coordinate",
        "standard_name" => "XN")
    YN_atts = Dict("units" => "m",
        "long_name" => "Y coordinate",
        "standard_name" => "YN")
    q_tot_atts = Dict("units" => "m3/s",
        "long_name" => "Total alongshore sediment transport",
        "standard_name" => "q_tot")
    kal_atts = Dict("units" => "-",
        "long_name" => "K parameter",
        "standard_name" => "K")
    RP_atts = Dict("units" => "-",
        "long_name" => "Pearson correlation coefficient",
        "standard_name" => "RP")
    RMSE_atts = Dict("units" => "m",
        "long_name" => "Root mean square error",
        "standard_name" => "RMSE")
    MSS_atts = Dict("units" => "-",
        "long_name" => "Mielke Skill Score",
        "standard_name" => "MSS")


    nccreate(output, "Y",
                "dim", (length(X0), length(YY)),
                atts = Y_atts)
    ncwrite(Ymdr, output, "Y")
    nccreate(output, "XN",
                "dim", (length(X0), length(YY)),
                atts = XN_atts)
    ncwrite(XN, output, "XN")
    nccreate(output, "YN",
                "dim", (length(X0), length(YY)),
                atts = YN_atts)
    ncwrite(YN, output, "YN")
    nccreate(output, "q_tot",
                "dim", (length(X0), length(YY)),
                atts = q_tot_atts)
    ncwrite(q_tot, output, "q_tot")
    nccreate(output, "Yi",
                "len", 1,
                atts = Yi_atts)
    ncwrite(yi, output, "Yi")
    nccreate(output, "K",
                "len", 1,
                atts = K_atts)
    # ncwrite([exp(popr[1])], output, "K")
    ncwrite([popr[1]], output, "K")
    nccreate(output, "RP",
                "len", 1,
                atts = RP_atts)
    ncwrite([aRP], output, "RP")
    nccreate(output, "RMSE",
                "len", 1,
                atts = RMSE_atts)
    ncwrite([aRMSE], output, "RMSE")
    nccreate(output, "MSS",
                "len", 1,
                atts = MSS_atts)
    ncwrite([aMSS], output, "MSS")



    println("\n\n****************Finished****************\n\n")
    # txt=['Tiempo para terminar:',num2str(floor(TT/60/60)),' h ',
    #     num2str(floor(60*(TT/60/60-floor(TT/60/60)))),' min ',
    #     num2str(60*(TT/60-floor(TT/60))),' s'];
    # disp(['Tiempo total de simulación: ', num2str((now-start_time)*24*60*60,'#.0f'), ' s']);

    #     # get the new positions

    

end