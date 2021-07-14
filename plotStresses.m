function plotStresses(noz)
    %%% This function plots nozzle wall stresses
    
    r = noz.r; l = noz.l;
    
    % nozzle contour
    yyaxis right
    ylabel('Radius (in)'); ylim([0 10]);
    plot(l, r, '.', 'MarkerEdgeColor', 'p', 'MarkerSize', 5); 
    surface([l;l], [r, r], [T_c;T_c], 'LineWidth', 4);

    % stresses
    yyaxis left
    ylabel('Temperautre (K)'); ylim([0 5000]);
    plot(l, noz.sig1, 'LineWidth', 3, 'LineColor', [.8, 0, 0.4])
    plot(l, noz.sig2c, 'LineWidth', 3, 'LineColor', [0.4, 0, 0.8])
    plot(l, noz.sig2g, 'LineWidth', 3, 'LineColor', [1, 0, 1])

    title('Steady State Temperatures'); 
    legend('Nozzle radius', 'Convective Heat Flux', 'Radiative Heat Flux', 'Net Heat Flux');
    hold off
end
