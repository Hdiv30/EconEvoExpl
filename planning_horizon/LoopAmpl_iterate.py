import os
import pandas as pd
import numpy as np

# previously calculated time (not planning horizon) of numerical optimisation so that first output is the same as with optimisation over long time horizon (T>=1000)
findTfile = pd.read_csv("minT.csv", sep='\t', header=0, index_col=None)


timehorizon = np.array([10, 100])  # time horizon
ps = np.array([0.3, 0.3])  # parameter for initial mean maturation
deltas = 0.05 ** (1 / timehorizon)  # discount factor. until 5 % of value
disc_rates = (1 - deltas) / deltas  # discount rate

for disc_rate, delta, p in zip(disc_rates, deltas, ps):

    # find time length in file:
    delta_round = np.round(delta, 2)
    while True:
        try:
            T = np.array(findTfile.loc[findTfile['delta'] == delta_round, 'T'])[0]
            print('T', T)
            break
        except IndexError:
            delta_round = np.round(delta_round + 0.01, 2)
        print('delta (rounded)', delta_round)


    T_str = str(int(T))
    delta_str = str(int(delta * 100))
    fileT = open('T.csv', 'w'); fileT.write(str(T)); fileT.close()
    filedelta = open('delta.csv', 'w'); filedelta.write(str(delta)); filedelta.close()
    filediscrate = open('disc_rate.csv', 'w'); filediscrate.write(str(disc_rate)); filediscrate.close()
    filep = open('p.csv', 'w'); filep.write(str(p)); filep.close()

    # run optimisation in AMPL/KNITRO:
    os.system('ampl model_ps_iterate.run')
