LORAAX: LOw Reynolds number Aircraft Aerodynamics with Xfoil

LORAAX is a 3D panel code fully coupled with a 2D integral boundary layer
method (Xfoil). Unlike traditional panel methods, It can predict nonlinearities
in the lift and pitching moment curves due to the onset of stall, separation
bubbles, and laminar-turbulent transition. Xfoil boundary layer calculations are
performed on the fly via the libxfoil library. Subsonic compressibility effects
are accounted for through the Prandtl-Glauert transformation. LORAAX is not
limited to a planar wake; the wake can be optionally deformed to the correct
force-free shape as part of the iterative solution. Since the boundary layer
method assumes 2D flow, only 0-sideslip conditions are currently supported, and
accuracy will suffer for very low aspect ratios and/or highly swept wings. The
method is applicable for model-scale Reynolds numbers all the way through
full-scale aircraft in subsonic (not transonic) flight.

LORAAX is designed to analyze wings, of which any number can be modeled in any
combination (for example, the traditional stabilizer would be modeled as a
second wing). Geometry is defined in a simple XML format. Output files include
forces and moments (totals and broken out for each wing), sectional quantities
(e.g., lift and drag distributions across the span), and visualizations of the
surface and wake. These files are designed to be read in ParaView.

Documentation and validation cases will be added. For now, you can find some
initial results and discussion on RCGroups:

https://www.rcgroups.com/forums/showthread.php?3037870-Numerical-and-experimental-results-for-low-Re-airfoils-and-wings/page6
