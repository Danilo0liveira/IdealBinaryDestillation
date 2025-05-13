import numpy as np
from model_parameters import *
from scipy.integrate import solve_ivp
from scipy.signal import StateSpace, lsim
import matplotlib.pyplot as plt

def column_model(x0, t, ts,lr):
    x = x0
    solution = solve_ivp(fun = ode, t_span=(t, t+ts), y0=x, method='RK45', args=(lr), dense_output=True)
    x = solution.y[:, -1]

    PV = x[0], x[-1]
    
    return PV, x

def ode(x=0,t=0,Lr=Lr):
    def column(x=0,t=0):

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
    return column(t,x)

def ss_column(x, Lr=0):
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

def lstate_model(xs, ys):
    ki = np.zeros(num_stages)
    Am = np.zeros(shape=(num_stages, num_stages))
    Bm = np.zeros(shape=(num_stages, 2))
    Cm = np.zeros(shape=(2, num_stages))
    Dm = np.zeros(shape=(2,2))

    for i in range(num_stages):
        ki[i] = a/(1 + (a-1)*xs[i])**2

    # overhead reciever
    Am[0,0] = -Vr/Md
    Am[0,1] = Vr/Md*ki[1]

    # top section
    for i in range(1,feed_position):
        Am[i,i] = -(Lr+Vr*ki[i])/Mt
        Am[i,i-1] = Lr/Mt
        Am[i,i+1] = Vr/Mt * ki[i+1]

    # feed stage
    Am[feed_position,feed_position-1] = Lr/Mt
    Am[feed_position,feed_position] = -(Ls+Vr*ki[feed_position])/Mt
    Am[feed_position,feed_position+1] = Vs/Mt * ki[i+1]

    # bottom section
    for i in range(feed_position+1, num_stages-1):
        Am[i,i] = -(Ls+Vs*ki[i])/Mt
        Am[i,i-1] = Ls/Mt
        Am[i,i+1] = Vs/Mt * ki[i+1]

    # Reboiler stage
    Am[num_stages-1,num_stages-2] = Ls/Mb
    Am[num_stages-1,num_stages-1] = -(B+Vs*ki[num_stages-1])/Mb

    Bm[0,0] = 0
    Bm[0,1] = 0
    
    for i in range(1,num_stages-1):
        Bm[i,0] = (xs[i-1]-xs[i])/Mt 
        Bm[i,1] = (ys[i+1]-ys[i])/Mt

    Bm[num_stages-1, 0] = (xs[i-1]-xs[i])/Mb 
    Bm[num_stages-1, 1] = (xs[num_stages-1]-ys[num_stages-1])/Mb

    Cm[0,0] = 1
    Cm[1,num_stages-1] = 1

    return Am, Bm, Cm, Dm

xs = np.array([0.98999601, 0.9850687,  0.97890511, 0.97123079, 0.96173054, 0.95005375,
 0.93582728, 0.91867875, 0.89827239, 0.87435803, 0.84683005, 0.81578811,
 0.78158589, 0.74485101, 0.70646153, 0.66747386, 0.62901176, 0.59213975,
 0.55775021, 0.52648835, 0.49872553, 0.47417105, 0.44554726, 0.41300468,
 0.37704318, 0.33853007, 0.29864467, 0.25874731, 0.22019995, 0.1841863,
 0.15157975, 0.1228865,  0.09826259, 0.07758188, 0.06052529, 0.04666705,
 0.03554394, 0.02670323, 0.01973116, 0.0142665,  0.01000397])

ys = np.zeros(shape=xs.shape)

for i in range(num_stages):
    ys[i] = (a*xs[i])/(1 + (a-1)*xs[i])

Am, Bm, Cm, Dm = lstate_model(xs, ys)
sys = StateSpace(Am, Bm, Cm, Dm)

t = np.linspace(0, 400, 400)
u = np.column_stack([-0.01*Lr*np.ones_like(t), 0.0*np.ones_like(t)])
x0 = np.zeros(num_stages)

t_out, y_out, x_out = lsim(sys, U=u, T=t, X0=x0)

plt.plot(t_out, y_out[:,0]+xs[0])
plt.plot(t_out, y_out[:,1]+xs[-1])
plt.xlabel("Tempo")
plt.ylabel("Saída")
plt.grid()
plt.show()