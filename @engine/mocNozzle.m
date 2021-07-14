function n = mocNozzle(eng, geometry, material, t_w, NS)
    %{
    This function characterizes combustion chamber + nozzle spline using
    method of characteristics via a few input parameters.

    Sources:

    [1] Grisson, W. M., "Liquid Film Cooling in Rocket Engines", Air Force 
    Office of Scientific Research.

    [2] moc source

    l_1 = [m]
    r_1 = [m]
    r_2 = [m] 
    D_t = [m] throat diameter
    D_ch = [m] combustion chamber diameter
    theta_c = [deg]
    theta_d = [deg]

    %}

    eng.t_w = t_w;
    eng.NS = NS;
       
    %%% NOZZLE MATERIAL PROPERTIES
    eng.lambda_w = material(1);
    eng.E = material(2);
    eng.a = material(3);
    eng.v = material(4);
        

    %%% COMBUSTION CHAMBER CONTOUR
    
    l_1 = geometry(1);
    r_1 = geometry(2);
    r_2 = geometry(3);
    D_t = geometry(4);
    D_ch = geometry(5);
    theta_c = geometry(6);
    theta_d = geometry(7);

    % see pg 59 and figure 6 in [1]
    l_2 = l_1 + r_1*sind(theta_c);
    D_2 = D_ch - 2*r_1*(1-cosd(theta_c));
    D_3 = D_t + 2*r_2*(1-cosd(theta_c));
    l_3 = l_2 + .5*(D_2-D_3)/tand(theta_c);
    l_t = l_3 + r_2*sind(theta_c);
    l_5 = l_t + r_2*sind(theta_d);
    D_5 = D_t + 2*r_2*(1-cosd(theta_d));
    
    eng.r = zeros(1, NS); n.r(1) = r_1;
    eng.l = zeros(1, NS);
    eng.x = zeros(1, NS);
    
    for i = 2:NS/2
        eng.l(i) = eng.l(i-1) + 2*l_5/NS;
        if eng.l(i) <= l_1
        elseif (eng.l(i) > l_1) && (eng.l(i) <= l_2)
            eng.r(i) = .5*D_ch - r_1 + sqrt(r_1^2 - (eng.l(i) - l_1)^2);
        elseif (eng.l(i) > l_2) && (eng.l(i) <= l_3)
            eng.r(i) = .5*D_2 - (eng.l(i) - l_2)*tand(theta_c);
        elseif (eng.l(i) > l_3) && (eng.l(i) <= l_5)
            eng.r(i) = .5*D_t + r_2 - sqrt(r_2^2 - (eng.l(i) - l_5)^2);
        else 
            disp("Something is wrong with your chamber contour calc");
        end
    end
    

    %%% NOZZLE CONTOUR 
    % TODO find NS/2 rs and ls using MOC
    
    

    %%% DOWNSTREAM WALL LOCATIONS
    for i = 2:NS
        eng.x(i) = eng.x(i-1) + ...
                   sqrt((eng.l(i)-eng.(i-1))^2 + (eng.r(i)-eng.r(i-1))^2);
    end
end