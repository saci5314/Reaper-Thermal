function eps = emissivity(T, p_c)
    %%% This fucntion finds emissivities of reaction products at T
    
    eps = 0;
    
    %    H2O     OH      CO2     CO      H2     
    w = [.41654  .03004  .18721  .34028  .01640];                           % mass fractions
    M = [.018    .017    .044    .028    .002];                             % [kg/M] molar masses     
    i = [5       2       3       2       1];                                % # radiation strips

    
    
    %%%%% By [3] in Bartz.m
    
    
    epsH2O = .825*(1 + rhoH2O)*
    
    eps = epsH20 + epsCO2 - delta_eps;
    
    
    
    %%%%% By [1] in Bartz.m
    
    %%% REACTION PRODUCT CHARACTERIZATON
    rpData = readmatrix("ecData.csv");
    
    Mbar = 1/sum(w./M);                                                     % [kg/M] avg molar mass of mixture
    
    a = 1; b = 1;
    for j = 1:5
        if j == 1; B = 1.09; else B = 1; end
        
        a = b;
        b = a + i(j);
        
        x = (w(j)/m(j))*MBar;                                               % mole fraction 
        
        m(j).pp = x*p_c;
        m(j).B = B;
        m(j).i = i(j);
        m(j).a_3000 = rpData(a:b, 1);
        m(j).omega_0 = rpData(a:b, 2);
        m(j).domega_s = rpData(a:b, 3);
        m(j).a = rpData(a:b, 4);
        m(j).delta = rpData(a);                         
    end
    
    
    
    %%% FIND EMISSIVITY
    
    eps_c = A_c/A_t;                                                        % chamber area contraction ratio
    V_cc = A_t*(L_c*eps_c + sqrt(A_t/pi)*cotd(theta)*(eps_c^(1/3) - 1)/3);  % [m^3] combustion chamber volume ([2] eq 4-5)
    A_cc = 2*L_c*sqrt(pi*eps_c*A_t) + cscd(theta)*(eps_c - a)*A_t;          % [m^2] combuston chamber area ([2] eq 4-6) 
    L = 3.6*V_c/A_c;
    
    ec = 0;                                                                 
    for j = 1:5                                                             % sorry no comments here ... 
        for i = 1:m(j).i                                                    % this math is hand-waivy... just see [1]
            domega = m(j).domga_s(i)*(1-m(j).delta) + m(j).a(i)*(1+m(j).delta)*sqrt(T/1000);
            omega1 = 0; omega2 = 0;
            switch m(j).delta
                case 1 
                    omega1 = m(j).omega_0(i) - domega/2; 
                    omega2 = m(j).omega_0(i) + m(j).domega/2;
                case 0
                    omega1 = m(j).omega_0(i) - (domega - m(j).domega_s(i)); 
                    omega2 = m(j).omega_0(i) + m(j).domega_s(i);            
            end
            
            R = int((3.742*omega*10^-12)/(exp(1.44*omega/T)-1), omega1, omega2, omega); 
            K = 3000*m(j).a_3000(i)/(domega*T);
            
            eps = eps + 10^4*m(j).B*(1-exp(-K)*m(j).pp*L)*R/(sig_sb*T^4);     % emissivity 
        end
    end
    
end