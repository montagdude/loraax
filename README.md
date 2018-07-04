LORAAX: LOw Reynolds number Aircraft Aerodynamics with Xfoil

LORAAX is a 3D panel code fully coupled with a 2D integral boundary layer
method (Xfoil). Unlike traditional panel methods, LORAAX can predict
nonlinearities in the lift and pitching moment curves due to the onset of
stall, separation bubbles, and laminar-turbulent transition. Xfoil boundary
layer calculations are performed on the fly via the libxfoil library, which is
very fast, robust, and does not require the user to manually generate polars.
Subsonic compressibility effects are accounted for through a Prandtl-Glauert
correction. Whereas many panel and vortex lattice codes assume a planar wake,
LORAAX starts with a planar wake and deforms it to its correct shape as part
of the iterative solution. Since the boundary layer method assumes 2D flow, only
0-sideslip conditions are currently supported, and accuracy of the method will
suffer for very low aspect ratios and/or highly swept wings. The method is
applicable for model-scale Reynolds numbers all the way through subsonic
full-scale aircraft.

LORAAX is designed to analyze wings, of which any number can be modeled in any
combination (for example, the traditional stabilizer would be modeled as a
second wing). Geometry is defined in a simple XML format. Output files include
forces and moments (totals and broken out for each wing), sectional quantities
(e.g., lift and drag distributions across the span), and visualizations of the
surface and wake. These files are designed to be read in ParaView.
