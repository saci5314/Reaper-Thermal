classdef engine
    properties (Access = public)
        %%% GEOMETRY
        NS;                                                                 % # stations
        r;                                                                  % [m] radial position
        l;                                                                  % [m] axial position
        x;                                                                  % [m] spline position

        %%% STATION/PLANE STATE     
        q_w;                                                                % [W/m^2] convective heat flux
        q_r;                                                                % [W/m^2] radiative heat flux
        q;                                                                  % [W/m^2] gas-side heat flux
        T_aw;                                                               % [K] adiabatic wall temp
        T_wg;                                                               % [K] gas-side wall temp
        T_wc;                                                               % [K] cool-side wall temp
        p_wg;                                                               % [N/m^2] gas pressure
        
        T_c;                                                                % [K] regen coolant temp
        T_f;                                                                % [K] film coolant temp
        T_fvap;                                                             % [K] vaporization temp of 
        dm_c;                                                               % [kg/] regen coolant massflow
        dm_f;                                                               % [kg/] film massflow
        
        M;                                                                  % local mach no.
        
        sig1;                                                               % [N/m^2] plane stress
        sig2c;                                                              % [N/m^2] cool-side stress
        sig2g;                                                              % [N/m^2] gas-side stress

        %%% NOZZLE GEOMETRY
        t_w;
        D_t;
        A_t;
        
        V_chamber;
        A_chamber;
        
        p_c
        
        %%% COOLING SYSTEMS
        n_cc;                                                               % number of regen channels
        w_cc;                                                               % [m] regen channel width
        h_cc;                                                               % [m] regen channel heigth
        l_f;                                                                % [m] film injection position
        t_tbc;                                                              % [m] tbc layer thickness(es)
        lambda_tbc;                                                         % tbc layer thermal conductivities

    end
    methods
        function n = nozzle(NS, chamber, wall)
            if nargin < 0 
                n.NS = NS;
                n.t_w = t_w;
                n.lambda_w = lambda_w;
            end    
        end
        
        n = rpaNozzle(n, rpaFile);                                          %%% DONE
        n = mocNozzle(n, l_1, r_1, r_2, D_t, D_ch, theta_c, theta_d, NS);   %%% NOT DONE NOT CRITICAL
        
        M = mach(r, gamma, A_t);                                            %%% DONE 
        eps = emissivity(n, T, p_c);                                        %%% DONE
            
        n = applyRegen(n);                                                  %%% NOT DONE
        n = applyFilm(n, type, dm, injLoc);                                 %%% DONE
        n = applyTBC(n, t, lambda);                                         %%% DONE
        
        n = Bartz(n);                                                       %%% NOT DONE
        n = nozStresses(n, material, t_j, t_f);                             %%% NOT DONE
    end
end