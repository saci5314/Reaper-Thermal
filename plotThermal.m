function plotThermal(noz)
    %%% This function plots nozzle wall and coolant temps and heat fluxes

    r = noz.r; l = noz.l;
    
    %%% TEMPERATURES
    figure
    hold on

    % nozzle contour with coolant temp
    yyaxis right
    ylabel('Radius (in)'); ylim([0 10]);
    plot(l, r, '.', 'MarkerEdgeColor', 'p', 'MarkerSize', 5); 
    surface([l;l], [r, r], [noz.T_c;noz.T_c], 'LineWidth', 4);

    % film coolant (UNFINISHED)
   

    % wall temps
    yyaxis left
    ylabel('Temperautre (K)'); ylim([0 5000]);
    plot(l, noz.T_wg, 'LineWidth', 3, 'LineColor', [1, 0.5, 0.5])
    plot(l, noz.T_wc, 'LineWidth', 3, 'LineColor', [1, 0.65, 0.35])

    title('Steady State Temperatures'); 
    legend('Nozzle radius', 'Hot-gas side', 'Coolant side');
    hold off

    %%% HEAT FLUX
    figure
    hold on

    % nozzle contour with coolant temp
    yyaxis right
    ylabel('Radius (in)'); ylim([0 10]);
    plot(l, r, '.', 'MarkerEdgeColor', 'p', 'MarkerSize', 5); 
    surface([l;l], [r, r], [T_c;T_c], 'LineWidth', 4);

    % fluxes
    yyaxis left
    ylabel('Temperautre (K)'); ylim([0 5000]);
    plot(l, noz.q, 'LineWidth', 3, 'LineColor', [.8, 0, 0.4])
    plot(l, noz.q_r, 'LineWidth', 3, 'LineColor', [0.4, 0, 0.8])
    plot(l, noz.q_w, 'LineWidth', 3, 'LineColor', [1, 0, 1])

    title('Steady State Temperatures'); 
    legend('Nozzle radius', 'Convective Heat Flux', 'Radiative Heat Flux', 'Net Heat Flux');
    hold off

end