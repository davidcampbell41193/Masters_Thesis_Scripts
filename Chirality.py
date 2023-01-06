#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 12:45:21 2018

@author: david
"""

import qutip as q 
import numpy as np
import matplotlib.pyplot as plt
 
n=5 # Number of cavity levels

g1=10.0
g2=g1
kappa=60.0
Delta=0.0
Delta_1=1.0
Delta_2=-Delta_1
Omega=46.0

Gamma = 4.0 * g1 * np.conj(g1) / kappa
G=g2/g1
J = -1.0j * Gamma * G / 2 

P=q.bell_state(state='00') # |B00> = 1 / sqrt(2)*[|0>|0>+|1>|1>]
M=q.bell_state(state='01') # |B01> = 1 / sqrt(2)*[|0>|0>-|1>|1>]
T=q.bell_state(state='10') # |B10> = 1 / sqrt(2)*[|0>|1>+|1>|0>]
S=q.bell_state(state='11') # |B11> = 1 / sqrt(2)*[|0>|1>-|1>|0>]

sm = [ q.tensor( q.sigmam() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmam() ) ]
sp = [ q.tensor( q.sigmap() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmap() ) ]
sx = [ q.tensor( q.sigmax() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmax() ) ]
sz = [ q.tensor( q.sigmaz() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmaz() ) ]

c_ops = [ np.sqrt(Gamma) * ( sm[0] + G*sm[1] ) ]
#c_ops=[2*g1*q.tensor( q.sigmam(), q.identity(2) )/np.sqrt(kappa) + 2*g2*q.tensor( q.identity(2) , q.sigmam() )/np.sqrt(kappa) ]
expect=[ S*S.dag(), T*T.dag() , M*M.dag() , P*P.dag() ]

q_coup = J * sm[0]*sp[1]
q_coup = q_coup + q_coup.dag()
H_qubits = Delta_1*sp[0]*sm[0] + Delta_2*sp[1]*sm[1] + Omega * ( sx[0] + sx[1] ) / 2 + q_coup

psi_0 = q.tensor( q.identity(2) , q.identity(2) )/4

tlist=np.linspace(0 , 300 , 1500)

output=q.mesolve(H_qubits, psi_0 , tlist , c_ops , expect )

states =output.states
#states_reduced = output.states

sing=[]
trip=[]
minus=[]
plus=[]
for y in states:
    sing.append( q.fidelity( S*S.dag() , y ) )
    trip.append( q.fidelity( T*T.dag() , y ) )
    minus.append( q.fidelity( M*M.dag() , y ) )
    plus.append( q.fidelity( P*P.dag() , y ) )
    

H_A = -Delta * q.tensor( q.identity(2) , q.identity(2) , q.num(n) )
H_q = Delta_1*q.tensor( q.sigmap()*q.sigmam() , q.identity(2) , q.identity(n) ) + Delta_2*q.tensor( q.identity(2) , q.sigmap()*q.sigmam() , q.identity(n) ) + Omega * ( q.tensor( q.sigmax() , q.identity(2) , q.identity(n) ) + q.tensor( q.identity(2) , q.sigmax() , q.identity(n) ) )/2
H_I = g1 * ( q.tensor( q.sigmam() , q.identity(2) , q.create(n) ) + q.tensor( q.sigmap(), q.identity(2) , q.destroy(n) ) ) + g2 * ( q.tensor( q.identity(2) , q.sigmam() , q.create(n) ) + q.tensor( q.identity(2) , q.sigmap() , q.destroy(n) ) ) + q.tensor( q_coup , q.identity(n) )

H=H_A+H_q+H_I

result=q.mesolve( H , q.tensor( psi_0 , q.fock_dm(n,0) ) , tlist, np.sqrt(kappa)*q.tensor(q.identity(2), q.identity(2), q.destroy(n)) , [] )

states=result.states

sing=[]
trip=[]
minus=[]
plus=[]

#P1=[]
#P2=[]
#P12=[]
#
for state in states:
    x=state.ptrace([0,1])
    sing.append( q.expect( S*S.dag() , x ) )
    trip.append( q.expect( T*T.dag() , x ) )
    minus.append( q.expect( M*M.dag() , x ) )
    plus.append( q.expect( P*P.dag() , x ) )

#fig, ax=plt.subplots()    
#ax.plot(tlist, P12 , label=r'$P_{12}$')
#ax.plot(tlist, P1 , label=r'$P_1$')
#ax.plot(tlist, P2 , label=r'$P_2$')
#ax.set_xlabel('time')
#ax.set_ylabel('Purity')
#plt.title( r'Full System $g=$%d $\kappa=$%d $\Delta=$%d $\Delta_1=-\Delta_2=$%d $\Omega=$%d' %(g1 , kappa , Delta , Delta_1 , Omega ) )
#plt.legend()
 
 
R = [ output.expect[0] , output.expect[1] , output.expect[2] , output.expect[3] ]
full = [sing , trip , minus , plus ]
    
import itertools
 
colors = [ 'k' , 'c' , 'g' , 'b' ]
cc = itertools.cycle(colors)
plot_lines = []

parameters = [ r"$|S\rangle$" , r"$|T\rangle$" , r"$|\Phi_- \rangle$" , r"$|\Phi_+\rangle$"    ]

for i in range(len(full)):
    
    d1 = full[i]
    d2 = R[i]
    
    plt.hold(True)
    c=next(cc)
    l1, = plt.plot( (g1**2) * tlist / kappa , d1 , '--' , color = c)
    l2, = plt.plot( (g1**2) * tlist / kappa , d2 , color = c )
    
    plot_lines.append( [ l1 , l2 ] )
    
legend1 = plt.legend( plot_lines[0] , ["Full" , "Reduced"] ,  )
plt.legend( [ l[1] for l in plot_lines ] , parameters , loc=4  )
plt.gca().add_artist(legend1)
plt.xlabel(r'$g^2 t/ \kappa$')
plt.ylabel('Fidelity')
plt.title( r'$g=$%d $\kappa=$%d $\delta_c=$%d $\delta_1=-\delta_2=$%d $\Omega=$%d' %(g1 , kappa , Delta , Delta_1 , Omega ) )
plt.savefig('/Users/davidcampbell633/Desktop/Asym_stabilization.pdf')