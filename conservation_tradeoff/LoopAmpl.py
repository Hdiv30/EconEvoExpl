import os
import pandas as pd
import numpy as np


targets = [0]  # t31 (max 53.49 = recovery from start if zero harvest) (min 50.61 = start) 
targets = [59]  # t81 (max 60.18 = recovery from start if zero harvest) (min 50.61 = start) 


T_target = 81  # year 2100 
# T_target = 31  # year 2050 

p = 0.2358  # --> mean maturation length as in data

for T_target in [31, 81]: 
    for disc_rate in [0.02]:
        delta = 1 / (1 + disc_rate)
        for target_lm in targets:
            delta_round2 = np.round(delta, 2) * 100
            T = 150

            T_str = str(int(T))
            delta_str = str(int(delta * 100))
            fileT = open('T.csv', 'w'); fileT.write(str(T)); fileT.close()
            filedelta = open('delta.csv', 'w'); filedelta.write(str(delta)); filedelta.close()
            fileT_target = open('T_target.csv', 'w'); fileT_target.write(str(T_target)); fileT_target.close()
            filetarget_lm = open('target_lm.csv', 'w'); filetarget_lm.write(str(target_lm)); filetarget_lm.close()
            filediscrate = open('disc_rate.csv', 'w'); filediscrate.write(str(disc_rate)); filediscrate.close()
            filep = open('p.csv', 'w'); filep.write(str(p)); filep.close()

            # os.system('ampl model_ps.run')
            os.system('ampl model_ps_sm.run')
