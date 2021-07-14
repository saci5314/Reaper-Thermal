# Reaper-Thermal

This is a reverse engieered version of RPA designed for running trade studies on the thermal and structural performance of liquid bipropellant rocket engines with complex cooling systems. The program inputs are:
- High-level engine geometry (chamber and throat radii,
throat radius of curvature, chamber length)
- Nozzle material properties
- High-level engine performance parameters (fuel/LOX
massflow, chamber pressure, stagnation temp, etc.)
- Reaction product mass fractions from NASA’s Chemical
Equilibrium and Applications (CEA)

Under the assumption of axisymmetric flow, an iterative solver is used to find steady-state wall temps along the nozzle spline. Those values are then used with core gas pressures to calculate principal stresses and strains in the nozzle wall. This analysis coupled with the object-oriented nature of the program will allow us to optimize the effectiveness of our cooling configuration, prolong the engines life cycle, and avoid rapid unscheduled disassembly.

## Main sources:
- Ponomarenko, A., “Thermal Analysis of Thrust Chambers”
- Bartz, D. R., “A Simple Equation for the Rapid Estimation
of Rocket Nozzle Convective Heat Transfer Coefficients”
- Grisson, D. R., “Liquid Film Cooling in Rocket Engines”
- Huang, D. H., Huzel, D. K., “Modern Engineering for the
Design of Liquid Propellant Rocket Engines”
- Boysan, M. E., “Analysis of Regenerative Cooling in Liquid
Propellant Rocket Engines”
