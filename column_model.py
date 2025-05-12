import numpy as np
from model_parameters import *
from scipy.integrate import solve_ivp

def column_model(x0, t, ts,lr):
    x = x0
    solution = solve_ivp(fun = out_column(lr), t_span=(t, t+ts), y0=x, method='RK45', args=(), dense_output=True)
    x = solution.y[:, -1]

    PV = x[0], x[-1]
    
    return PV, x

def out_column(Lr):
    def column(t,x):

        for i in range(num_stages):
            y[i] = (a*x[i])/(1 + (a-1)*x[i])

        #overhead reciever
            dxi[0] = 1/Md * (Vr*y[1] - (D+R)*x[0])

        #rectifying section (2 to nf-1)
        for i in range(1, feed_position):
            dxi[i] = 1/Mt* (Lr*x[i-1] + Vr*y[i+1] - Lr*x[i] - Vr*y[i])
        
        #feed stage
        dxi[feed_position] = 1/Mt * (Lr*x[feed_position-1] + Vs*y[feed_position+1] + F*zf - Ls*x[feed_position] - Vr*y[feed_position])
        
        #rectifying section (nf+1 to ns-1)
        for i in range(feed_position+1, num_stages-1):
            dxi[i] = 1/Mt* (Ls*x[i-1] + Vs*y[i+1] - Ls*x[i] - Vs*y[i])

        #reboiler 
        dxi[num_stages-1] = 1/Mb * (Ls*x[num_stages-2] - B*x[num_stages-1] - Vs*y[num_stages-1])

        return dxi
    return column