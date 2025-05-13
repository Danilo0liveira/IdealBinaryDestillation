import numpy as np

num_stages=41 # número de estágios
feed_position=20 # estágio de alimentação
F=1 # fluxo molar de alimentação (mol/min)
zf=0.5 # composição de componete leve na alimentação (mol)
qf = 1 # qualidade do stream de alimentação

Vs=3.206 # fluxo molar de vapor na parte superior
Vr= Vs + F*(1-qf) # fluxo molar de vapor na retificação

R=2.706 # fluxo molar no refluxo
a = 1.5 # constante de equilibrio

Ls=R+F*qf # fluxo molar de líquido na seção inferior
Lr=R # fluxo molar de líquido na seção superior

B=Ls-Vs # fluxo molar de material no fundo
D=0.5 # fluxo molar na destilação


Mb=5 # retenção molar no fundo
Md=5 # retenção molar no receptor
Mt=0.5 # retenção molar nas trilhas

dxi = np.zeros(num_stages)
y = np.zeros(num_stages)