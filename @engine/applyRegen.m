function eng = applyRegen(eng, T_i, dm, h, w, t, coolant, n)
    %{ 
    This function applies a regenerative cooling system to the engine
    
    T_i = [K] coolant channel input temp
    dm_i = [kg/s] coolant massflow rate
    h = [m] channel heights: [inlet, throat, exit] 
    w = [m] channel widths: [inlet, throat, exit]
    t = [m] [fin thickness, jacket thickness]
    coolant = [coolant thermal conductivity, vaporization enthalpy]
    n = number of cooling channels
    
    %}
    %%%%% CHANNEL STUFF
    eng.regen = true;

    eng.dm_c = dm;    
    eng.T_c(1)      = T_i;
    eng.lambda_c    = coolant(1);
    eng.Q_cVap      = coolant(2);
    eng.cp_c        = coolant(3);
    mu_c            = coolant(4);
    k_c             = coolant(5);                                           % 
    
    eng.t_f         = t(1);
    eng.t_j         = t(2);
    
    eng.n_cc        = n;
    
    
    %%%%% INTERPOLATE CHANNEL HEIGHT AND WIDTH 
    dw1 = w(2) - w(1);
    dw2 = w(3) - w(2);
    dC1 = 2*pi*(min(eng.r) - eng.r(1));
    dC2 = 2*piI(eng.r(end) - min(eng.r));

    %%% INLET
    eng.h_cc(1) = h(1);
    eng.w_cc(1) = w(1);
    
    %%% NOZZLE
    j = 2;
    while eng.r(j) > min(eng.r)
        dC              = 2*pi*(eng.r(j) - eng.r(j-1));d
        dw              = dw1*dC/dC1;
        
        eng.w_cc(j)     = eng.w_cc(j-1) + dw;
        eng.h_cc(j)     = A/eng.w_cc(j);

        j = j + 1;
    end
    
    %%% THROAT
    eng.h_cc(j) = h(2);
    eng.w_cc(j) = w(2);
    
    %%% COMBUSTION CHAMBER
    j = j + 1;
    while j < eng.NS
        dC              = 2*pi*(eng.r(j) - eng.r(j-1));d
        dw              = dw2*dC/dC2;
        
        eng.w_cc(j)     = eng.w_cc(j-1) + dw;
        eng.h_cc(j)     = A/eng.w_cc(j);
        
        j = j + 1;
    end
    
    %%% OUTLET
    eng.h_cc(j) = h(3);
    eng.w_cc(j) = w(3);
    
    %%%%% OTHER STUFF
    eng.Re_c        = 2*eng.mDot/(mu_c.*(eng.h_cc+eng.w_cc));
    eng.Pr_c        = eng.cp_c*mu_c/k_c;
    

end