# Ice-Stream-Flowline-Model-Simple

# This is a numerical implementation of the flowline model described in Robel et al., JGR, 2014 for MATLAB. It solves for ice thickness,
# grounding line position, velocity, till void ratio, unfrozen till thickness and ice temperature, given a bed topography (defined in
# Base.m and dBasedx.m), and other parameters. You can find a detailed mathematical description of the numerics in Robel et al., JGR,
# 2014. 

# Plug and play instructions:
# The code runs in MATLAB, and should work simply by running the driver file, which is called Flowline_init.m. The first few code blocks in the driver file will set parameters and model configuration, pre-allocate storage for saving output and initialize variables before beginning forward integration of the model itself. All the model parameters and settings are stored in the "parameters" structure (with the syntax parameters.name). The most commonly changed parameters and model settings are set in Flowline_init.m, and the rest are all set in setparams_init.m.

# The model calculates ice thickness, grounding line position, velocity, temperature, till void ratio and unfrozen till thickness over all grounded ice down the center line of an ice stream (from the ice divide to the grounding line). It does this using an operator-splitting where ice thickness, grounding line position, velocity and temperature are calculated used a few different implicit methods and the other variables are calculated explicitly. The modeling strategy and numerical details can be found in the supplement to Robel et al. 2014.

# There is no ice shelf in this model. The version here assumes that ice stream width is constant throughout the flowline. The model currently has the capability of adaptive time stepping. If something starts changing really fast, then the model will chop the time step in half. This time step length will always be tracked by parameters.dtau, so if that is critical for the sea level calculation, be careful that you take it into account.

# As the model is written, it needs the bed topography to be provided as an equation (that is a differentiable function of the form b = b(x) where b is the bed elevation as a function of distance from the ice divide). The bed topography function is specified in Base.m and the x-derivative of the bed topography is specified in dBasedx.m. Currently, the specified bed topography is simply linear, like:
b(x) = b_0 + b_x * x
where b_0 is the bed elevation at the ice divide and b_x is the bed slope. Thus, the spatial derivative is easy, as it is just b_x everywhere. The bed topography is referenced to sea level, so at any given time b=0 is sea level. These functions specifying bed topography and its derivative are used to calculate the driving stress in the velocity solver and (more critically) the grounding line position in the thickness/GL solver. 

# Right now, if you just hit "run" on Flowline_init.m, it should run and continuously plot a bunch of diagnostics as the model evolves. The most immediately interesting variables are plotted in the left column, where you will see (from the top) the basal velocity in the ice stream as a function of distance along the flowline (where the grounding line is at sigma=1), a side-view of the current ice stream thickness and bed along the flowline, and the grounding line position as a function of time. The most recent profile of a given quantity is plotted in red, with all the old profiles plotted in black.

# If you just look at the grounding line position, you will see it undergo oscillations with period of approx. 1500 years. You can play around with the period of these oscillations or turn them off entirely by changing (among other parameters) the surface temperature (parameters.T_s), geothermal heat flux (parameters.G), or accumulation rate (parameters.a_c). The dependence of oscillation characteristics on the various parameters is explored in depth in Robel et al. 2013.
