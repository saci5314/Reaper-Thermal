function eng = Bartz(eng)
    %{ 
    This function performs a 1D steady state heat tranfer analysis across
    the nozzle wall via the Bartz method. The method has been modified to
    consider arbitrarily complex regen, film, and tbc cooling.
    
    Sources:
    
    [1] Ponomarenko, A., "Thermal Analysis of Thrust Chambers", RPA: Tool
        for Rocket Propulsion Analyis. https://www.rocket-propulsion.com/d
        ownloads/pub/RPA_ThermalAnalysis.pdf.

    [2] Bartz, D. R. , "A Simple Equation for Rapid Estimation of Rocket
        Nozzle Convective Heat Transfer Coefficients", Jet Propulsion
        Laboratory, California Insititute of Technology.

    [3] Grisson, W. M., "Liquid Film Cooling in Rocket Engines", Air Force 
        Office of Scientific Research.
    
    [4] Huang, D. H., Huzel, D. K., "Modern Engineering for Design of Liquid
        Propellant Rocket Engines", Rocketdyne Division of Rockwell
        International.

    [5] Boysan, M. E., "Anlaysis of Regenerative Cooling in Liquid
        Propellant Rocket Engines".

    %}



    %%%                             %%%
    %%%      ENGINE PARAMETERS      %%%
    %%%                             %%%
    
    
    
    %%% CONSTANTS
    g = 9.8;                                                                                % [m/s^2] gravitational acceleration
    sigma_sb = 5.67e-8;                                                                     % [W/(m^2*K^4)] Stefan-Boltzman constant
    
    %%% FILM STUFF
    et = .15;                                                                               % !!! make this a function of axial position !!!
    K_t = 1 + 10.2*et;                                                                      % boundary layer turbulence correction factor
    K_m = (eng.Mc/eng.Mg)^.14;                                                              % coolant-gas molecular weight ratio correction factor

    %%% OTHER STUFF
    D_t = 2*min(eng.r);                                                                     % [m] throat diameter
    A_t = pi*min(eng.r)^2;                                                                  % [m^2] throat area
        
    cStar = eng.p_c*A_t*g/eng.wDot;                                                         % [m/s] characteristic velocity
    
    ec_r_T_0 = emissivity(eng, eng.T_0);                                                    % emissivity of reaction products at hot gas temp

    
    
    %%% PREALLOCATION
    T_wg = zeros(1, eng.NS);
    
    
    
    
    
    %%%                 %%%
    %%%     EXECUTE     %%%
    %%%                 %%%
    
    
    
    tic; fin = 0;
    while fin ~= 1
        %%% RADIATIVE FLUX
        updateT_aw();
        update_q_r();
 
        %%% CONVECTIVE FLUX
        update_q_w();
        
        %%% NET GAS SIDE FLUX
        eng.q = eng.q_w + eng.q_r;
        
        %%% FILM EVOLUTION
        if eng.film
            updateT_bl();
            correct_q_w();                               
        end
        
        %%% COOLANT SIDE 
        if eng.regen
            updateT_c(); 
            updateT_wc();
        end

        %%% GAS SIDE WALL TEMP 
        updateT_wg();

        %%% CONVERGENCE CRITERIA 
        for j = 1:eng.NS
            delta = abs((eng.T_wg(j)-T_wg(j))/ebng.T_wg(j)); 
            
            if delta > 0.05*T_wg(j); break
            elseif j == eng.NS; fin = 1;
            end
        end
 
        eng.T_wg = T_wg; 
    end
    toc;
    
    

    %%%                  %%%
    %%%     FUNCTIONS    %%%
    %%%                  %%%
    
    
    
    %%%%% UPDATE RADIATIVE HEAT TRANSFER (GO FOR TESTING)
    
    function update_q_r()    
        for j = 1:eng.NS
            % via [1]: 
            ec_r_T_aw   = emissivity(eng.T_aw(j));                                          % emissivity of reaction products at wall temp
            ec_aw       = eng.ec_w/(1-(1-eng.ec_w)*(1-ec_r_T_aw));                          % emissivity of wall

            eng.q_r(j)  = ec_aw*sigma_sb*(ec_r_T_0*T_0^4 - ec_r_T_aw*eng.T_aw(j)^4);        % [W/m^2] radiative heat flux 

            % via [5]:
            %{
            P = P_cr/(1+C1)^(gamma/(gamma-1));                                              % [kg/cm^2] what the fuck      
            p_CO2 = P*molFrac_CO2;                                                          % [kg/cm^2] partial pressure maybe???
            p_H2O = P*molFrac_H2O;                                                          % [kg/cm^2] ^^^

            L_e = 1.2*n.r(j);

            q_r_CO2 = 3.5*(p_CO2*L_e)^(1/3)*((n.T_aw(j)/100)^3.5-(n.T_wg(j)/100)^3.5);      % [W/m^2] radiative gas side heat flux due to CO2([5] eq 2.10)
            q_r_H2O = 3.5*p_H2O^.8*L_e^.6*((n.T_aw(j)/100)^3-(n.T_wg(j)/100)^3);            % [W/m^2] radiative flux due to H20  ([5] eq 2.11)
            n.q_r(j) = q_r_CO2+q_r_H2O;                                                     % [W/m^2] radiative heat flux ([5] eq 2.9)
            %}
        end
    end



    %%%%% UPDATE CONVECTIVE HEAT TRANSEFR (GO FOR TESTING)
    
    function update_q_w()
        for j = 1:eng.NS
            C           = 1 + eng.M(j)^2*(gamma-1)/2;
            sig_h_g     = (.5*eng.T_wg(j)*C/eng.T_0 + .5)^-.68*C^-.12;                      % h_g correction factor
            h_g         = (.026/D_t^.2)*(eng.mu_g^.2*eng.cp_g/eng.Pr_g^.6)...
                          *(eng.p_c*g/cStar)^.8*(D_t*eng.R_t)^.1...
                          *(A_t/(pi*eng.r(j)^2))^.9*sig_h_g;                                % [W/m^2] heat transfer coefficient ([1], [2] eq 7, [4] eq 4-13, [5]) 

            eng.q_w(j)  = h_g*(eng.T_aw(j) - eng.T_wg); 
        end
    end

    

    %%%%% FILM-CORRECT CONVECTIVE HEAT TRANSFER (GO FOR TESTING)
    
    function correct_q_w()
        %%% CORRECT FOR NEW BOUNDARY LAYER TEMP OVER GAS FILM LENGTH
        update_q_w();
        
        %%% OMIT FOR LIQUID FILM LENGTH
        for j = 1:eng.NS
            if eng.T_f ~= 0; eng.q_w(j) = 0; end
        end
        
        %%% CORRECTED NET HEAT FLUX
        eng.q = eng.q_w + eng.q_r;
    end
    
    


    %%%%% UPDATE FILM COOLANT TEMPS (GO FOR TESTING)
    
    function updateT_bl()
        eng.T_f(eng.NS)     = eng.T_c(eng.NS);
        dmDot_bl            = 0;
        j                   = eng.liqStart;
        
        %%% LIQUID FILM FLOW
        while eng.mDot_bl(j) >= 0 && eng.r(j) ~= min(eng.r)
            j = j - 1;

            % BOUNDARY LAYER TEMPS
            dT_f                = 2*pi*eng.r(j)*eng.q(j)/eng.mDot_bl(j)*eng.dx(j);              % [K] increase due to core gas flow
            eng.T_f(j)          = eng.T_f(j-1) + dT_f;
            
            % CHECK FOR VAPORIZATION
            if eng.T_f(j) >= eng.T_cVap
                dmDot_bl           = 2*pi*eng.r(j)*eng.q(j)/eng.Q_cVap*eng.dx(j);               % [kg/s] massflow of evaporation
                eng.mDot_bl(j)     = eng.mDot_bl(j) - dmDot_bl;                                 % [kg/s] liquid boundary layer massflow 
            end
        end
        
        eng.T_f(1:j)         = 0;
        eng.liqEnd           = j;
        eng.T_aw(j)          = eng.T_f(j+1);
      
        %%% GASEOUS FILM FLOW
        while j > 1
            j = j - 1;
            
            % MASSFLOW
            dmDot_1             = .1963*K_t*G*(eng.mu_g*eng.mDot_bl(j))^-.25;                   % [kg/s] boundary layer massflow increase due to free stream gas entrainment
            
            dmDot_bl            = dmDot_bl + dmDot_1;                                           % [kg/s] projected boundary layer massflow
            
            dmDot_2             = .1963*K_t*G*(eng.mu_g*eng.mDot_bl(j-1))^-.25;                 % [kg/s] projected dmDot_1 by 2nd order Runge Kutta
            dmDot_e             = .5*(dmDot_1 + dmDot_2);                                       % [kg/s] 2nd order Runge Kutta-corrected dmDot_1
            dmDot_c             = -eng.mDot_bl(j)*(eng.r(j) - eng.r(j+1))/eng.r(j);             % [kg/s] increase in boundary layer massflow due to chamber contraction
            
            eng.mDot_bl(j-1)    = eng.mDot_bl(j-1) + dmDot_e + dmDot_c;                         % [kg/s] gaseous boundary layer massflow
            
            % BOUNDARY LAYER TEMPS            
            
            dT_e                = eng.dmDot_e*(T_r - T)/...
                                (eng.dmDot_bl(j) + dmDot_c*((eng.cp_c/eng.cp_g)/K_m - 1));      % [K/m] dT due to free stream gas entrainment
            dT_r                = eng.q_r(j)/(eng.cp_g*eng.mDot_bl(j))*eng.dx(j);               % [K/m] dT due to thermal radiation 
            
            eng.T_aw(j)         = eng.T_aw(j-1) + dT_e + dT_r;                                  % [K] boundary layer temp
        end 
    end


    
    %%%%% UPDATE REGEN COOLANT TEMPS (GO FOR TESTING)
    
    function updateT_c() 
        for j = eng.NS-1:1 
            dT_c            = 2*pi*eng.r(j)*eng.q(j)*enj.dx(j)/(eng.mDot_c*eng.cp_c);           % [K/m] coolant dT/dx                 
            
            eng.T_c(j)      = eng.T_c(j-1)+dT_c;                                                % [K] new regen coolant temps
        end   
    end



    %%%%% UPDATE ADIABATIC WALL TEMPS (GO FOR TESTING)
    
    function updateT_aw()
        for j = 1:eng.NS
            C               = eng.M(j)^2*(eng.gamma-1)/2;                                                          
            eng.T_aw(j)     = eng.T_0*(1+eng.Pr^.33*C)/(1+C);                                   % [K] adiabatic wall temp
        end
    end

    

    %%%%% UPDATE COOLANT-SIDE WALL TEMPS (GO FOR TESTING)
    
    function updateT_wc()
        for j = 1:eng.NS
            Nu              = 0.0185*eng.Re_c(j)^.8*eng.Pr_c*.4*(eng.T_c(j)-eng.T_wc(j))^.1;    % Nusselt number for methane [1]
            A_wet           = eng.n_cc*eng.w_cc(j)*eng.h_cc(j);                                 % [m^2] average total flow area of cooling passages
            P_wet           = 2*eng.n_cc*(eng.w_cc(j) + eng.h_cc(j));                           % [m]average total perimeter of cooling passages
            d_e             = 4*A_wet/P_wet;                                                    % [m] cooling channel equivalent diameter at station
            h_c             = Nu*eng.lambda_c/d_e;                                              % coolant-side heat transfer coefficient
            
            eng.T_wc(j)     = eng.q(j)/h_c+eng.T_c(j);                                          % [K] new coolant side wall temps    
        end
    end



    %%% UPDATE GAS SIDE WALL TEMPS (GO FOR TESTING)
    function updateT_wg()    
        C_TBC = eng.t_w/eng.lambda_w;
        if eng.t_tbc ~= 0                                                                       % check for TBC layers
            for k = 1:length(eng.t_tbc)
                C_TBC = C_TBC + eng.t_tbc(k)/eng.lambda_tbc(k);                            
            end
        end
        for j = 1:eng.NS
            T_wg(j) = eng.q(j)*C_TBC+eng.T_wc(j);                                               % [K] new gas side wall temps 
        end
    end


end