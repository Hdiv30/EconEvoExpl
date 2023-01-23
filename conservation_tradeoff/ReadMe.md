
# LoopAmpl.py
- loops over time point where target needs to be reached and actual target (x-axis)
- passes those parameters to the .run script
- for both 31 (2050) and 81 (2100) it is enough to numerically optimise until T=150 (this T has nothing to do with the planning horizon)

# model_ps.run
- uses model_ps.mod and model.dat data which include all parameters from calibration
- ps = perfect size selection

# model_ps_sm.run
- uses model_ps_sm.mod and model.dat data which include all parameters from calibration
- ps_sm = perfect size and maturity selection

outputs are stored in the folder out
