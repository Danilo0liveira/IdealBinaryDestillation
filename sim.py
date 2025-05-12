import numpy as np
import matplotlib as plt

from column_model import column_model
from model_parameters import num_stages

save_data = True # option to save simulation data
tmax = 4800 # simulation duration time
ts = 1 # sampling time
sim_speed = 0.01 # simulation speed

state = np.array([0.98999601, 0.9850687,  0.97890511, 0.97123079, 0.96173054, 0.95005375,
 0.93582728, 0.91867875, 0.89827239, 0.87435803, 0.84683005, 0.81578811,
 0.78158589, 0.74485101, 0.70646153, 0.66747386, 0.62901176, 0.59213975,
 0.55775021, 0.52648835, 0.49872553, 0.47417105, 0.44554726, 0.41300468,
 0.37704318, 0.33853007, 0.29864467, 0.25874731, 0.22019995, 0.1841863,
 0.15157975, 0.1228865,  0.09826259, 0.07758188, 0.06052529, 0.04666705,
 0.03554394, 0.02670323, 0.01973116, 0.0142665,  0.01000397])

time = np.arange(0, tmax, ts)
listinha_pv1 = list()
listinha_pv2 = list()


for k in range(time.size):
    PV, state = column_model(state, time[k], ts)
    print(PV)
    listinha_pv1.append(PV[0])
    listinha_pv2.append(PV[1])

plt.plot(PV[0])