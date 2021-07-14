function n = applyTBC(n, t, lambda)
    %{ 
    This function models thermal barrier coating(s) in your nozzle.

    Source:

    Ponomarenko, A., "Thermal Analysis of Thrust Chambers", RPA: Tool
    for Rocket Propulsion Analyis. https://www.rocket-propulsion.com/d
    ownloads/pub/RPA_ThermalAnalysis.pdf.

    %}

    n.t_tbc = t;
    n.lambda_tbc = lambda;
end
