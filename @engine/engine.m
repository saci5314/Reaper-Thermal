classdef engine
    properties (Access = public)
        %%% NOZZLE CHARACTERISTICS
        NS;                                                                 % # stations
        r;                                                                  % [m] radial position
        l;                                                                  % [m] axial position
        x;                                                                  % [m] spline position
        t_w;                                                                % [m] chamber wall thickness
        R_t;                                                                % [m] throat radius of curvature
        
        lambda_w;                                                           % [] wall thermal conductivity
        ec_w;                                                               % emissivity
        E;                                                                  % [N/m^2] elastic modulus
        a;                                                                  % [W/(m*K)] thermal conductivity
        v;                                                                  % Poisson's ratio
        
        %%% PERFORMANCE CHARACTERISTICS
        p_c;                                                                % [N/m^2] chamber pressure
        T_0;                                                                % [K] stagnation Temp
        gamma;                                                              % specific heat ratio
        cp_g;                                                               % chamber specific heat
        wDot;                                                               % [N/s] prop weight flowrate 
        Pr;                                                                 % Prandtl No.
        mu_g;                                                                 % [] flow viscocity
        
        %%% PLANE STATE     
        q_w;                                                                % [W/m^2] convective heat flux
        q_r;                                                                % [W/m^2] radiative heat flux
        q;                                                                  % [W/m^2] gas-side heat flux
        T_aw;                                                               % [K] adiabatic wall temp
        T_wg;                                                               % [K] gas-side wall temp
        T_wc;                                                               % [K] cool-side wall temp
        p_wg;                                                               % [N/m^2] gas pressure
        
        T_c;                                                                % [K] regen coolant temp
        T_f;                                                               % [K] liquid film temp
        mDot_c;                                                             % [kg/m] regen coolant massflow
        mdot_bl;                                                            % [kg/m] boundary layer/film massflow
        
        M;                                                                  % local mach numbers
        
        sig1;                                                               % [N/m^2] plane stress
        sig2c;                                                              % [N/m^2] cool-side stress
        sig2g;                                                              % [N/m^2] gas-side stress   
        
        %%% COOLING SYSTEMS
        regen = false;
        film = false;
        
        lambda_c;                                                           % [W/(m^2*K)] coolant thermal conductivity
        Q_cVap;                                                             % [J/mol?] coolant vaporization enthalpy
        T_cVap;                                                             % [K] coolant vaporization temp
        Re_c;                                                               % Reynold's number 
        Pr_c;                                                               % Prandtl number 
        cp_c;                                                               % average specific heat
     
        n_cc;                                                               % number of regen channels
        w_cc;                                                               % [m] regen channel width
        h_cc;                                                               % [m] regen channel heigth
        
        liqStart;                                                           % film injection index 
        liqEnd;                                                             % liquid film evaporation point
        
        t_tbc;                                                              % [m] tbc layer thicknesses
        lambda_tbc;                                                         % [W/(m^2*K)] tbc layer thermal conductivities
                
    end
    methods
        function eng = engine(params, nozzle, type, nozInput, NS)             %%% CRITICAL TODO
            if nargin > 0              
                
                %%% NOZZLE CHARACTERISTICS
                eng.t_w         = nozzle(1);
                eng.lambda_w    = nozzle(2);
                eng.E           = nozzle(3);
                eng.a           = nozzle(4);
                eng.v           = nozzle(5);
                
                %%% GENERATE NOZZLE SPLINE
                eng.NS = NS;    
                
                if      type == "rpa";   eng = rpaNozzle(eng, nozInput);
                elseif  type == "moc";   eng = mocNozzle(eng, nozInput);
                end
                
                %%% PERFORMANCE CHARACTERISTICS
                eng.p_c         = params(1);
                eng.T_0         = params(2);
                eng.cp_g        = params(3);
                eng.gamma       = params(4);
                eng.wDot        = params(5);
                
                eng.Pr          = 4*eng.gamma/(9*eng.gamma - 5);            % (Modern 4-15)
                eng.mu          = 46.6e-10*M*eng.T_0;                       % (Modern 4-16)
                eng.Re          = 2*eng.wDot/(pi*eng.r.*eng.mu);                        
                eng.M = mach(eng);
            end    
        end
        
        eng = rpaNozzle(eng, rpaFile);                                      %%% GO FOR TESTING
        eng = mocNozzle(eng, geometry);                                     %%% NOT CRITICAL TODO
        
        M   = mach(eng);                                                    %%% GO FOR TESTING
        eps = emissivity(eng, T);                                      %%% GO FOR TESTING
            
        eng = applyRegen(eng, T_i, dm_i, h, w, t, coolant, n);              %%% GO FOR TESTING
        eng = applyFilm(eng, dm, injLoc);                                   %%% GO FOR TESTING
        eng = applyTBC(eng, t, lambda);                                     %%% GO FOR TESTING
        
        eng = Bartz(eng);                                                   %%% GO FOR TESTING
        eng = nozStresses(eng, material);                                   %%% GO FOR TESTING
        

    end
end