
### LoopAmpl_iterate.py 
- loops over planning horizon (discount rate x-axis) and initial mean size at maturity (y-axis)
- passes those parameters to the run script

### model_ps_iterate.run
- uses discount rate and initial maturity
- optimises over a time horizon of T (from minT.csv - this has numerical reasons and has nothing to do with the planning horizon. With T=1000 (default) there are numerical errors, which is why we iterate over rounds and only store the first result. details see supplementary material)
- uses model_ps.mod and model.dat data which include all parameters from calibration

outputs are stored in the folder out
