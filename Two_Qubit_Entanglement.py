#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 10:53:21 2018

@author: davidcampbell63
"""

"""
Spyder Editor

This is a temporary script file.
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

P=q.bell_state(state='00') # |B00> = 1 / sqrt(2)*[|0>|0>+|1>|1>]
M=q.bell_state(state='01') # |B01> = 1 / sqrt(2)*[|0>|0>-|1>|1>]
T=q.bell_state(state='10') # |B10> = 1 / sqrt(2)*[|0>|1>+|1>|0>]
S=q.bell_state(state='11') # |B11> = 1 / sqrt(2)*[|0>|1>-|1>|0>]

c_ops=[2*g1*q.tensor( q.sigmam(), q.identity(2) )/np.sqrt(kappa) + 2*g2*q.tensor( q.identity(2) , q.sigmam() )/np.sqrt(kappa) ]
expect=[ S*S.dag() , T*T.dag() , M*M.dag() , P*P.dag() ]

H_qubits = Delta_1/2*q.tensor( q.sigmaz() , q.identity(2) ) + Delta_2/2*q.tensor( q.identity(2) , q.sigmaz() ) + Omega * ( q.tensor( q.sigmax() , q.identity(2) ) + q.tensor( q.identity(2) , q.sigmax() ) ) / 2

psi_0 = q.tensor( q.qeye(2) , q.qeye(2) )/4

tlist=np.linspace(0 , 2500 , 2000)

output=q.mesolve( H_qubits, psi_0 , tlist , c_ops , expect ) 

H_A = Delta * q.tensor( q.identity(2) , q.identity(2) , q.num(n) )
H_q = Delta_1/2*q.tensor( q.sigmaz() , q.identity(2) , q.identity(n) ) + Delta_2/2*q.tensor( q.identity(2) , q.sigmaz() , q.identity(n) ) + Omega * ( q.tensor( q.sigmax() , q.identity(2) , q.identity(n) ) + q.tensor( q.identity(2) , q.sigmax() , q.identity(n) ) ) / 2
H_I = g1 * ( q.tensor( q.sigmam() , q.identity(2) , q.create(n) ) + q.tensor( q.sigmap(), q.identity(2) , q.destroy(n) ) ) + g2 * ( q.tensor( q.identity(2) , q.sigmam() , q.create(n) ) + q.tensor( q.identity(2) , q.sigmap() , q.destroy(n) ) )    

H=H_A+H_q+H_I

result=q.mesolve( H , q.tensor( psi_0 , q.fock_dm(n,0) )  , tlist, np.sqrt(kappa)* ( q.tensor(q.identity(2), q.identity(2), q.destroy(n)) ) , [] )

states=result.states

sing=[]
trip=[]
minus=[]
plus=[]
cavity = []

for state in states:
    cavity.append( q.expect( q.tensor( q.qeye(2) , q.qeye(2) , q.create(n) ) , state ) )
    x=q.ptrace( state , [0,1] )
    sing.append( q.expect( S*S.dag() , x ) )
    trip.append( q.expect( T*T.dag() , x ) )
    minus.append( q.expect( M*M.dag() , x ) )
    plus.append( q.expect( P*P.dag() , x ) )
    #sing.append( q.expect( q.tensor(S*S.dag(), q.coherent_dm(n,alpha))  , state ) )
    #trip.append( q.expect( q.tensor(T*T.dag(), q.coherent_dm(n,alpha))  , state ) )
    #minus.append( q.expect( q.tensor(M*M.dag(), q.coherent_dm(n,alpha)) , state ) )
    #plus.append( q.expect( q.tensor(P*P.dag(), q.coherent_dm(n,alpha)) , state ) )
 
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
#plt.savefig('/Users/davidcampbell633/Desktop/Sym_stabilization.pdf')
 
#ax.plot(tlist, sing , color='c' , label=r'$|S\rangle|0\rangle_c$' )
#ax.plot(tlist, trip , color='g' , label=r'$|T\rangle|0\rangle_c$')
#ax.plot(tlist, minus , color='r' , label=r'$|\Phi_-\rangle|0\rangle_c$')
#ax.plot(tlist, plus , color='b' , label=r'$|\Phi_+\rangle|0\rangle_c$')
#ax.set_xlabel('time')
#ax.set_ylabel('Fidelity')
#plt.title( r'Full System $g=$%d $\kappa=$%d $\Delta=$%d $\Delta_1=-\Delta_2=$%d $\Omega=$%d' %(g1 , kappa , Delta , Delta_1 , Omega ) )
#plt.legend()
#plt.savefig('/Users/davidcampbell633/Desktop/sing.pdf')

#plt.savefig('Gap_Difference.pdf')
#q.matrix_histogram( q.steadystate(H_qubits , c_ops) )
#rho_ss = q.steadystate( H , [np.sqrt(kappa)*q.tensor(q.identity(2), q.identity(2), q.destroy(n))] )
#q.matrix_histogram( rho_ss.ptrace([0,1]))
#q.matrix_histogram(rho_ss.ptrace(2)) 
#

# USE THIS FOR popultation plots

#for j in range( n ):
#    environ=[]
#    for state in states:
#        x=state.ptrace(2)
#        environ.append(q.expect(q.fock_dm(n,j), x ))
##    fig2 , ax=plt.subplots()
#    plt.plot(tlist, environ , label='%d' %j)
#plt.xlabel('time')
#plt.ylabel('Fidelity')
#plt.legend()
#plt.show()


# 
#beta_1= 2 * Omega * Delta_1
#beta_2= 2 * Omega * Delta_2 

