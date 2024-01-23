function BreakingPropagation(H1,T1,DIR1,h1,ANGbati,breakType)
    ###########################################################################    
    # Propagation of waves using linear theory assuming rectilinear & parallel bathymetry
    #    
    # INPUT:
    # H1:        wave height.
    # T1:        wave period.
    # DIR1:      wave direction. Nautical convention.
    # h1:        depth of wave conditions.
    # ANGbati:   bathymetry angle; the normal of the shoreline. Cartesian notation
    # breakType: type of breaking condition. Spectral | monochromatic.
    #    
    # OUTPUT:
    # H2:        wave height during breaking. Wave period is assumed invariant due to linear theory
    # DIR2:      wave direction during breaking. Nautical convention.
    # h2:        depth of breaking
    ###########################################################################    

    if breakType == "mono"
        Bcoef=0.78
    elseif breakType == "spectral"
        Bcoef = 0.45
    end

    DIRrel = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1),ANGbati)

    h2l0 = H1./Bcoef; # # initial condition for breaking depth
        
    
    H2 = zeros(length(H1),1); 
    DIR2 = zeros(length(DIR1),1); 
    h2 = zeros(length(H1),1)
    
    
    H2[h2l0.>=h1] = H1[h2l0.>=h1]
    DIR2[h2l0.>=h1] = DIR1[h2l0.>=h1]
    h2[h2l0.>=h1] = h2l0[h2l0.>=h1]  # # check that the initial depth is deeper than the breaking value
    
    H2[H1.<=0.1] = H1[H1.<=0.1]
    DIR2[H1.<=0.1] = DIR1[H1.<=0.1]
    h2[H1.<=0.1] = h2l0[H1.<=0.1]  
    
    propProf = (abs.(DIRrel).<=90) .&& (H1.>0.1) .&& (h2l0.<h1)
    propProf = vec(propProf)
    
    if sum(propProf)>0
        myFun(x) = LinearShoalBreak_Residual(x, H1[propProf], T1[propProf], DIR1[propProf], h1[propProf], ANGbati[propProf], Bcoef)
        
        # lb = zeros(size(h2l0[propProf])) .+ 0.1
        # ub = zeros(size(h2l0[propProf])) .+ 30

        # try
        #     nlboxsolve(myFun,h2l0[propProf], lb, ub, xtol=1e-1,ftol=1e-1).zero
        # catch
        #     println("\n\n", H1, "\n\n")
        #     println(T1, "\n\n")
        #     println(DIR1, "\n\n")
        #     println(h1, "\n\n")
        #     println(ANGbati, "\n\n")
        #     println(DIRrel, "\n\n")
        # end
        # h2l = nlboxsolve(myFun,h2l0[propProf], lb, ub, xtol=1e-1,ftol=1e-1).zero
        # h2l = nlsolve(myFun,h2l0[propProf]).zero

        # println("h2l = ",h2l)
        h2l = optimize.newton_krylov(myFun,h2l0[propProf]; method="minres")
        H2l, DIR2l = LinearShoalBreak_ResidualVOL(h2l, H1[propProf],T1[propProf], DIR1[propProf], h1[propProf], ANGbati[propProf], Bcoef);                
        H2[propProf] = H2l
        DIR2[propProf] = DIR2l
        h2[propProf] = h2l
    end
        
    return H2, DIR2, h2
end

function GroupCelerity(L,T,h)
    ###########################################################################    
    # CELERITY GROUP
    # L: wave lenght.
    # T: wave period.
    # h: depth of wave conditions.
    ###########################################################################       
    
    c = L./T
    k = 2 .*pi./L
    N = 1.0.+ 2.0 .*k.*h./sinh.(2.0 .*k.*h)
    Cg = c ./ 2.0 .* N
    
    return Cg
end

function hunt(T,d)

    ###########################################################################    
    # Wave lenght from Hunt's approximation
    #
    # INPUT:
    # T:     Wave Peak period.
    # d:     Local depth.
    # OUTPUT:
    # L:     Wave length.
    ###########################################################################   
   
    
    g=9.81; #[m/s^2]

    G= (2. .* pi ./T ) .^2 .*(d./g)
    
    # p=Polynomial([.067,.0864,.4622,.6522,1])
    p=Polynomial([1.0,.6522,.4622,.0864,.0675])
    
    F = G .+ 1.0 ./p.(G)

    L= T.*(g.*d./F).^.5
    
    return L
end

function LinearShoal(H1, T1, DIR1, h1, h2, ANGbati)
    ###########################################################################    
    # Wave shoaling & refraction applying linear theory with parallel; rectilinear bathymetry.
    #    
    # INPUT:
    # H1:        initial wave height.
    # T1:        wave period.
    # DIR1:      initial wave direction. Nautical convention.
    # h1:        initial depth of wave conditions.
    # h2:        final depth of wave conditions.
    # ANGbati:   bathymetry angle; the normal of the shoreline. Cartesian convention
    #
    # OUTPUT:
    # H2:        wave height during breaking. Wave period is assumed invariant due to linear theory.
    # DIR2:      wave direction during breaking. Nautical convention.
    ###########################################################################

    
    relDir1 = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1),ANGbati)

    # L1, _ =RelDisp(h1,T1)
    # L2, _ =RelDisp(h2,T1)
    L1 = hunt(T1,h1)
    L2 = hunt(T1,h2)
    CG1 = GroupCelerity(L1,T1,h1)
    CG2 = GroupCelerity(L2,T1,h2)
    relDir2 = Snell_Law(L1,L2,relDir1)
    KS = sqrt.(CG1./CG2)
    KR = sqrt.(cos.(relDir1.*pi./180.)./cos.(relDir2.*pi./180.))
    H2 = H1.*KS.*KR
    DIR2 = cartesianDir2nauticalDir(abs_angle_cartesian(relDir2,ANGbati))
    
    return H2, DIR2
end

function LinearShoalBreak_Residual(h2l, H1, T1, DIR1, h1, ANGbati, Bcoef)

    H2l, _ = LinearShoal(H1, T1, DIR1, h1, h2l, ANGbati)
    H2comp = h2l.*Bcoef
    res = H2l-H2comp

    return res
end

function LinearShoalBreak_ResidualVOL(h2l, H1, T1, DIR1, h1, ANGbati, Bcoef)

    H2l, DIR2l = LinearShoal(H1, T1, DIR1, h1, h2l, ANGbati)
    H2comp = h2l.*Bcoef
    res = H2l-H2comp

    return H2l, DIR2l
end

function RelDisp(h, T)
    
    g=9.81
    
    L=hunt(T,h)
    
    Li=hunt[T,h]
    error =1
    while error.>1E-6
        L=g.*T.^2.0./(2*pi).*tanh(2*pi.*h./Li)
        error = sum(sum(abs(L-Li)))
        Li=L
    end
    
    C=g.*T./(2*pi).*tanh.(2*pi.*h./L)
    
    return L,C
end

function RU2_Stockdon2006(slope,hs0,tp)
    ###########################################################################    
    # Run up 2# STOCKDON 2006
    #
    # INPUT:
    # slope:  Beach Slope in swash zone H:V. if slope is V:H all the slope terms multiply
    # hs0:     Significant wave height in deep water.
    # tp:     Peak period.
    #
    # OUTPUT:
    # runup2:   Run-up exceed 2%.
    ###########################################################################
    
    g=9.81
    hs0 = 0.5.*(hs0[2:end]+hs0[1:end-1])
    tp = 0.5.*(tp[2:end]+tp[1:end-1])
    L0 = g.*tp.*tp./2.0/pi
    slope = 1.0./slope; # # V:H
    setup = 0.35.*slope.*(hs0.*L0).^(1.0/2.)
    infgr = (hs0.*L0.*(0.563.*slope.*slope.+0.004)).^(1.0./2.)./2.
    runup2 = 1.1 .* (setup.+infgr) # # eq 19 Stockdon 2006
    return runup2
end

function Snell_Law(L1,L2,alpha1)
    ###########################################################################    
    # Wave refraction using snell law.
    #    
    # INPUT:
    # L1:     initial wave length.
    # L1:     final wave length.
    # alpha1: initial wave dir. Cartesian notation.
    #
    # OUTPUT:
    # alpha1: final wave dir. Cartesian notation.
    ###########################################################################    
    alpha=asin.(L2.*sin.(alpha1.*pi./180.0)./L1).*180.0./pi
    
    return alpha
end