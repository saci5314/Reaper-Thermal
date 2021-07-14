function eng = nozStresses(eng)
%{ 
This function calculates nozzle wall stresses

Source: Martinez-Sancches, M., "Lecture 9: Liquid Cooling", 16.512 - Rocket 
        Propulsion, https://ocw.mit.edu/courses/aeronautics-and-astronautic
        s/16-512-rocket-propulsion-fall-2005/lecture-notes/lecture_9.pdf.
%}

%%% EVALUATE STRESSES
for j = 1:eng.NS
    
    % via linear system of equations (eq. 5, 8, and 13 from source)
    % DOUBLE CHECK MATHS
    b = [eng.E*eng.a*(eng.T_wg(j) - eng.T_wc(j))/(1 - eng.v),
         eng.E*eng.a*(eng.T_wc(j) - eng.T_c(j))/(1 - eng.v),
         eng.p_c*2*eng.r(j) + 2*eng.p_c*eng.w_cc(j)];
         
    A = [0, 1, -1;
         1, -1, 0;
         2*eng.t_j, eng.t_w, 0];                            
    
    sigStation = linsolve(A, b);
    
    n.sig1(j) = sigStation(1); % plane stress
    n.sig2c(j) = sigStation(2); % cold side stress
    n.sig2g(j) = sigStation(3); % hot side stress
end





