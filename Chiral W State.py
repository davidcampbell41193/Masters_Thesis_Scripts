#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 04:04:34 2019

@author: davidcampbell633
"""


import qutip as q 
import numpy as np
import matplotlib.pyplot as plt
 
n=10 # Number of cavity levels

g1=4.0
g2=g1
kappa=1.0
Delta=0.0
Delta_1=0.1
Delta_2=1.1
Delta_3=-0.1
Delta_4=-1.1
Omega=60.0

Gamma = 8.0 * g1 * np.conj(g1) / kappa
G=g2/g1
J = -1.0j * Gamma * G / 2

P=q.bell_state(state='00') # |B00> = 1 / sqrt(2)*[|0>|0>+|1>|1>]
M=q.bell_state(state='01') # |B01> = 1 / sqrt(2)*[|0>|0>-|1>|1>]
T=q.bell_state(state='10') # |B10> = 1 / sqrt(2)*[|0>|1>+|1>|0>]
S=q.bell_state(state='11') # |B11> = 1 / sqrt(2)*[|0>|1>-|1>|0>]

sm = [q.tensor( q.sigmam() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmam() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmam() ) ]
sp = [q.tensor( q.sigmap() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmap() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmap() ) ] 
sx = [q.tensor( q.sigmax() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmax() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmax() ) ]
sz = [q.tensor( q.sigmaz() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmaz() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmaz() ) ]

c_ops = [ np.sqrt(Gamma) * ( sm[0] + sm[1] + sm[2] ) ]
#c_ops = [ np.sqrt(Gamma) * ( sm[0] + G*sm[1] ) , np.sqrt(Gamma) * ( sm[1] + sm[2] ) , np.sqrt(Gamma) * ( sm[2] + sm[3] ) ]
#c_ops=[2*g1*q.tensor( q.sigmam(), q.identity(2) )/np.sqrt(kappa) + 2*g2*q.tensor( q.identity(2) , q.sigmam() )/np.sqrt(kappa) ]
#expect=[ P*P.dag() , M*M.dag() , T*T.dag() , S*S.dag() ]

q_coup = J*sm[0]*(sp[1]+sp[2]) + J*sm[1]*(sp[2])
#q_coup = J*sm[0]*sp[1] + J*sm[1]*sp[2] + J*sm[2]*sp[3]
q_coup = q_coup + q_coup.dag()
H_qubits = Delta_1*sp[0]*sm[0] + Delta_2*sp[1]*sm[1] + Delta_3*sp[2]*sm[2] + Delta_4*sp[3]*sm[3] + Omega * ( sx[0] + sx[1] + sx[2] + sx[3] ) + q_coup

psi_0 = q.tensor( q.qeye(2) , q.qeye(2) , q.qeye(2) , q.qeye(2) )/16

tlist=np.linspace(0 , 0.5 , 1000)

output=q.mesolve(H_qubits, psi_0 , tlist , c_ops , [] )

states =output.states
#states_reduced = output.states
#
#singr=[]
#tripr=[]
#minusr=[]
#plusr=[]
#for y in states_reduced:
#    singr.append( q.fidelity( S*S.dag() , y ) )
#    tripr.append( q.fidelity( T*T.dag() , y ) )
#    minusr.append( q.fidelity( M*M.dag() , y ) )
#    plusr.append( q.fidelity( P*P.dag() , y ) )
    

#H_A = -Delta * q.tensor( q.identity(2) , q.identity(2) , q.num(n) )
#H_q = Delta_1*q.tensor( q.sigmap()*q.sigmam() , q.identity(2) , q.identity(n) ) - Delta_2*q.tensor( q.identity(2) , q.sigmap()*q.sigmam() , q.identity(n) ) + Omega * ( q.tensor( q.sigmax() , q.identity(2) , q.identity(n) ) + q.tensor( q.identity(2) , q.sigmax() , q.identity(n) ) )
#H_I = g1 * ( q.tensor( q.sigmam() , q.identity(2) , q.create(n) ) + q.tensor( q.sigmap(), q.identity(2) , q.destroy(n) ) ) + g2 * ( q.tensor( q.identity(2) , q.sigmam() , q.create(n) ) + q.tensor( q.identity(2) , q.sigmap() , q.destroy(n) ) ) + q.tensor( q_coup , q.identity(n) )
#
#H=H_A+H_q+H_I
#
#result=q.mesolve( H , q.tensor( psi_0 , q.fock_dm(n,0) ) , tlist, np.sqrt(kappa)*q.tensor(q.identity(2), q.identity(2), q.destroy(n)) , [] )
#
#states=result.states
#
#sing=[]
#trip=[]
#minus=[]
#plus=[]

P=[]

P1=[]
P2=[]
P3=[]
P4=[]

P12=[]
P13=[]
P14=[]
P23=[]
P24=[]
P34=[]

P123=[]
P124=[]
P134=[]
P234=[]


S=[]

S1=[]
S2=[]
S3=[]
S4=[]

S12=[]
S13=[]
S14=[]
S23=[]
S24=[]
S34=[]

S123=[]
S124=[]
S134=[]
S234=[]

for state in states:
    P1.append( ( q.ptrace( state , 0 )**2 ).tr() )
    P2.append( ( q.ptrace( state , 1 )**2 ).tr() )
    P3.append( ( q.ptrace( state , 2 )**2 ).tr() )
    P4.append( ( q.ptrace( state , 3 )**2 ).tr() )
    
    P12.append( ( q.ptrace( state, [0,1] )**2 ).tr() )
    P13.append( ( q.ptrace( state, [0,2] )**2 ).tr() )
    P14.append( ( q.ptrace( state, [0,3] )**2 ).tr() )
    P23.append( ( q.ptrace( state, [1,2] )**2 ).tr() )
    P24.append( ( q.ptrace( state, [1,3] )**2 ).tr() )
    P34.append( ( q.ptrace( state, [2,3] )**2 ).tr() )
    
    P123.append( ( q.ptrace( state, [0,1,2] )**2 ).tr() )
    P124.append( ( q.ptrace( state, [0,1,3] )**2 ).tr() )
    P134.append( ( q.ptrace( state, [0,2,3] )**2 ).tr() )
    P234.append( ( q.ptrace( state, [1,2,3] )**2 ).tr() )
    
    P.append( ( state**2 ).tr() )
    
    S1.append( q.entropy_vn( q.ptrace( state , 0 ) ) )
    S2.append( q.entropy_vn( q.ptrace( state , 1 ) ) )
    S3.append( q.entropy_vn( q.ptrace( state , 2 ) ) )
    S4.append( q.entropy_vn( q.ptrace( state , 3 ) ) )
    
    S12.append( q.entropy_vn( q.ptrace( state, [0,1] ) ) )
    S13.append( q.entropy_vn( q.ptrace( state, [0,2] ) ) )
    S14.append( q.entropy_vn( q.ptrace( state, [0,3] ) ) )
    S23.append( q.entropy_vn( q.ptrace( state, [1,2] ) ) )
    S24.append( q.entropy_vn( q.ptrace( state, [1,3] ) ) )
    S34.append( q.entropy_vn( q.ptrace( state, [2,3] ) ) )
    
    S123.append( q.entropy_vn( q.ptrace( state, [0,1,2] ) ) )
    S124.append( q.entropy_vn( q.ptrace( state, [0,1,3] ) ) )
    S134.append( q.entropy_vn( q.ptrace( state, [0,2,3] ) ) )
    S234.append( q.entropy_vn( q.ptrace( state, [1,2,3] ) ) )
    
    S.append( q.entropy_vn( state ) )

    #sing.append( q.expect( S*S.dag() , x ) )
    #trip.append( q.expect( T*T.dag() , x ) )
    #minus.append( q.expect( M*M.dag() , x ) )
    #plus.append( q.expect( P*P.dag() , x ) )

fig, ax1=plt.subplots()    
ax1.plot( Gamma*tlist , P1 , label=r'$P_1$')
ax1.plot( Gamma*tlist , P2 , label=r'$P_{2}$')
ax1.plot( Gamma*tlist , P3 , label=r'$P_{3}$')
ax1.plot( Gamma*tlist , P4 , label=r'$P_4$')
ax1.set_xlabel(r'$\Gamma t$')
ax1.set_ylabel('Purity')
plt.title( r'$\Gamma=8\frac{|g|^2}{\kappa}$%d $\Delta_1=$%d $\Delta_2=$%d $\Delta_3=$%d $\Delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) )
plt.legend()
#plt.savefig('/Users/davidcampbell633/Dropbox/David/Notes_Latex/Tetramer_Figs/A_Sigle_Qubit_Populations.pdf')

fig, ax2=plt.subplots()    
ax2.plot( Gamma*tlist , P12 , label=r'$P_{12}$')
ax2.plot( Gamma*tlist , P13 , label=r'$P_{13}$')
ax2.plot( Gamma*tlist , P14 , label=r'$P_{14}$')
ax2.plot( Gamma*tlist , P23 , label=r'$P_{23}$')
ax2.plot( Gamma*tlist , P24 , label=r'$P_{24}$')
ax2.plot( Gamma*tlist , P34 , label=r'$P_{34}$')
ax2.set_xlabel(r'$\Gamma t$')
ax2.set_ylabel('Purity')
plt.title( r'$\Gamma=8\frac{|g|^2}{\kappa}$%d $\Delta_1=$%d $\Delta_2=$%d $\Delta_3=$%d $\Delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) )
plt.legend()
#plt.savefig('/Users/davidcampbell633/Dropbox/David/Notes_Latex/Tetramer_Figs/A_Double_Qubit_Populations.pdf')

fig, ax3=plt.subplots()    
ax3.plot( Gamma*tlist , P123 , label=r'$P_{123}$')
ax3.plot( Gamma*tlist , P124 , label=r'$P_{124}$')
ax3.plot( Gamma*tlist , P134 , label=r'$P_{134}$')
ax3.plot( Gamma*tlist , P234 , label=r'$P_{234}$')
ax3.plot( Gamma*tlist , P , label=r'$P$')
ax3.set_xlabel(r'$\Gamma t$')
ax3.set_ylabel('Purity')
plt.title( r'$\Gamma=8\frac{|g|^2}{\kappa}$%d $\Delta_1=$%d $\Delta_2=$%d $\Delta_3=$%d $\Delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) )
plt.legend()
#plt.savefig('/Users/davidcampbell633/Dropbox/David/Notes_Latex/Tetramer_Figs/A_Triple_Qubit_Populations.pdf')

#fig, ax1=plt.subplots()    
#ax1.plot( Gamma*tlist , S1 , label=r'$S_1$')
#ax1.plot( Gamma*tlist , S2 , label=r'$S_{2}$')
#ax1.plot( Gamma*tlist , S3 , label=r'$S_{3}$')
#ax1.plot( Gamma*tlist , S4 , label=r'$S_4$')
#ax1.set_xlabel(r'$\Gamma t$')
#ax1.set_ylabel('Entropy')
#plt.title( r'$\Gamma=8\frac{|g|^2}{\kappa}$%d $\Delta_1=$%d $\Delta_2=$%d $\Delta_3=$%d $\Delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) )
#plt.legend()
##plt.savefig('/Users/davidcampbell633/Dropbox/David/Notes_Latex/Tetramer_Figs/A_Sigle_Qubit_Populations.pdf')
#
#fig, ax2=plt.subplots()    
#ax2.plot( Gamma*tlist , S12 , label=r'$S_{12}$')
#ax2.plot( Gamma*tlist , S13 , label=r'$S_{13}$')
#ax2.plot( Gamma*tlist , S14 , label=r'$S_{14}$')
#ax2.plot( Gamma*tlist , S23 , label=r'$S_{23}$')
#ax2.plot( Gamma*tlist , S24 , label=r'$S_{24}$')
#ax2.plot( Gamma*tlist , S34 , label=r'$S_{34}$')
#ax2.set_xlabel(r'$\Gamma t$')
#ax2.set_ylabel('Entropy')
#plt.title( r'$\Gamma=8\frac{|g|^2}{\kappa}$%d $\Delta_1=$%d $\Delta_2=$%d $\Delta_3=$%d $\Delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) )
#plt.legend()
##plt.savefig('/Users/davidcampbell633/Dropbox/David/Notes_Latex/Tetramer_Figs/A_Double_Qubit_Populations.pdf')
#
#fig, ax3=plt.subplots()    
#ax3.plot( Gamma*tlist , S123 , label=r'$S_{123}$')
#ax3.plot( Gamma*tlist , S124 , label=r'$S_{124}$')
#ax3.plot( Gamma*tlist , S134 , label=r'$S_{134}$')
#ax3.plot( Gamma*tlist , S234 , label=r'$S_{234}$')
#ax3.plot( Gamma*tlist , S , label=r'$S$')
#ax3.set_xlabel(r'$\Gamma t$')
#ax3.set_ylabel('Entropy')
#plt.title( r'$\Gamma=8\frac{|g|^2}{\kappa}$%d $\Delta_1=$%d $\Delta_2=$%d $\Delta_3=$%d $\Delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) )
#plt.legend()
#plt.savefig('/Users/davidcampbell633/Dropbox/David/Notes_Latex/Tetramer_Figs/A_Triple_Qubit_Populations.pdf')
 


#full = [ output.expect[0] , output.expect[1] , output.expect[2] , output.expect[3] ]
#R = [plus , minus , trip , sing ]
#    
#import itertools
# 
#colors = ['b' , 'g' , 'c' , 'r' ]
#cc = itertools.cycle(colors)
#plot_lines = []
#
#parameters = [ r"$|\Phi_+\rangle$" , r"$|\Phi_\rangle$" , r"$|T\rangle$" , r"$|S\rangle$"  ]
#
#for i in range(len(full)):
#    
#    d1 = full[i]
#    d2 = R[i]
#    
#    plt.hold(True)
#    c=next(cc)
#    l1, = plt.plot( tlist , d1 , '--' , color = c)
#    l2, = plt.plot( tlist , d2 , color = c )
#    
#    plot_lines.append( [ l1 , l2 ] )
#    
#legend1 = plt.legend( plot_lines[0] , ["Full" , "Reduced"] )
#plt.legend( [ l[1] for l in plot_lines ] , parameters , loc=4  )
#plt.gca().add_artist(legend1)
#
#plt.xlabel('time')
#plt.ylabel('Fidelity')
#plt.title('$=$%f$j$ $g=$%d $\kappa=$%d $\Delta=$%d $\Delta_1=-\Delta_2=$%d $\Omega=$%d' %(J.imag ,  g1 , kappa , Delta , Delta_1 , Omega ))
##plt.savefig('/Users/davidcampbell633/Desktop/sing.pdf')
#plt.show()