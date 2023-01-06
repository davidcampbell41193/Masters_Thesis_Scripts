# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import qutip as q 
import numpy as np
import matplotlib.pyplot as plt

Gamma = 10.0 # 8.0 * g1 * np.conj(g1) / kappa
G=1.0 #g2/g1
J = 0.0 #-1j *  Gamma * G / 2

Delta_1=1.0#0.6*Gamma
Delta_2=1.0#0.4*Gamma
Delta_3=-1.0#-0.6*Gamma
Delta_4=-1.0#-0.4*Gamma

Omega=8.0

P=q.bell_state(state='00') # |B00> = 1 / sqrt(2)*[|0>|0>+|1>|1>]
M=q.bell_state(state='01') # |B01> = 1 / sqrt(2)*[|0>|0>-|1>|1>]
T=q.bell_state(state='10') # |B10> = 1 / sqrt(2)*[|0>|1>+|1>|0>]
S=q.bell_state(state='11') # |B11> = 1 / sqrt(2)*[|0>|1>-|1>|0>]

sm = [q.tensor( q.sigmam() , q.qeye(2) , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmam() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmam() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.qeye(2) , q.sigmam() ) ]
sp = [q.tensor( q.sigmap() , q.qeye(2) , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmap() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmap() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.qeye(2) , q.sigmap() ) ] 
sx = [q.tensor( q.sigmax() , q.qeye(2) , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmax() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmax() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.qeye(2) , q.sigmax() ) ]
sz = [q.tensor( q.sigmaz() , q.qeye(2) , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.sigmaz() , q.qeye(2) , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.sigmaz() , q.qeye(2) ) , q.tensor( q.qeye(2) , q.qeye(2) , q.qeye(2) , q.sigmaz() ) ]

c_ops = [ np.sqrt(Gamma) * ( sm[0] + sm[1] + sm[2] + sm[3] ) ]
#c_ops = [ np.sqrt(Gamma) * ( sm[0] + G*sm[1] ) , np.sqrt(Gamma) * ( sm[1] + sm[2] ) , np.sqrt(Gamma) * ( sm[2] + sm[3] ) ]
#c_ops=[2*g1*q.tensor( q.sigmam(), q.identity(2) )/np.sqrt(kappa) + 2*g2*q.tensor( q.identity(2) , q.sigmam() )/np.sqrt(kappa) ]
#expect=[ P*P.dag() , M*M.dag() , T*T.dag() , S*S.dag() ]

q_coup = J*sm[0]*(sp[1]+sp[2]+sp[3]) + J*sm[1]*(sp[2]+sp[3]) + J*sm[2]*sp[3]
#q_coup = J*sm[0]*sp[1] + J*sm[1]*sp[2] + J*sm[2]*sp[3] + J*sm[3]*sp[0]
q_coup = q_coup + q_coup.dag()
H_qubits = Delta_1*sp[0]*sm[0] + Delta_2*sp[1]*sm[1] + Delta_3*sp[2]*sm[2] + Delta_4*sp[3]*sm[3] + Omega * ( sx[0] + sx[1] + sx[2] + sx[3] ) + q_coup

e = q.fock(2,0)
g = q.fock(2,1)
psi_0 = q.tensor(g,g,g,g) # q.tensor( q.qeye(2) , q.qeye(2) , q.qeye(2) , q.qeye(2) )/16 #

Npts=1000
tlist=np.linspace(0 , 100.0 , Npts)

output=q.mesolve(H_qubits, psi_0 , tlist , c_ops , [] )

states =output.states

sing=[]
trip=[]
minus=[]
plus=[]

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




#========== Plot =========================#
fig, (ax1, ax2, ax3)=plt.subplots(1,3,sharey=True)
#my_suptitle=fig.suptitle( r'$\Gamma=4\frac{g^2}{\kappa}=$%d $\delta_a=$%d $\delta_b=$%d $\delta_3=$%d $\delta_4=$%d $\Omega=$%d' %(Gamma , Delta_1 , Delta_2 , Delta_3, Delta_4 , Omega ) , fontsize=14 , verticalalignment='bottom' )
#my_suptitle=fig.suptitle( r'$\Gamma=4\frac{g^2}{\kappa}=$%d $\delta=1$ $\Omega=$%d' %(Gamma , Omega ) , fontsize=14 , verticalalignment='bottom' )


ax1.plot( Gamma*tlist , P1 , label=r'$P_1$')
ax1.plot( Gamma*tlist , P2 , label=r'$P_{2}$')
ax1.plot( Gamma*tlist , P3 , label=r'$P_{3}$')
ax1.plot( Gamma*tlist , P4 , label=r'$P_4$')
ax1.set_xlabel(r'$\Gamma t$')
ax1.set_ylabel('Purity')
ax1.legend(loc=4)
 
ax2.plot( Gamma*tlist , P12 , label=r'$P_{12}$', color='k')
ax2.plot( Gamma*tlist , P13 , ':', label=r'$P_{13}$')
ax2.plot( Gamma*tlist , P14 ,'--' , label=r'$P_{14}$')
ax2.plot( Gamma*tlist , P23 , '--' , label=r'$P_{23}$')
ax2.plot( Gamma*tlist , P24 , ':' , label=r'$P_{24}$')
ax2.plot( Gamma*tlist , P34 ,  label=r'$P_{34}$' , color='r')
#ax2.plot( Gamma*tlist , P2 , ':' , label=r'P', color='k' )
ax2.set_xlabel(r'$\Gamma t$')
#ax2.set_ylabel(r'Purity')
#ax2.set_title( r'$\delta_a \;  -\delta_a \; \delta_b \; -\delta_b$ '  )
ax2.set_title( r'Symmetric'  )
ax2.legend(loc=4)
   
  
ax3.plot( Gamma*tlist , P123 , label=r'$P_{123}$')
ax3.plot( Gamma*tlist , P124 , label=r'$P_{124}$')
ax3.plot( Gamma*tlist , P134 , label=r'$P_{134}$')
ax3.plot( Gamma*tlist , P234 , label=r'$P_{234}$')
ax3.plot( Gamma*tlist , P , ':' , label=r'$P$' , color='k')
ax3.set_xlabel(r'$\Gamma t$')
#ax3.set_ylabel('Purity')
plt.legend(loc=4)
plt.tight_layout(w_pad=0)
#fig.savefig('/Users/davidcampbell633/Desktop/Sym_Homo.pdf', dpi=fig.dpi, bbox_inches='tight',bbox_extra_artists=[my_suptitle])
fig.savefig('/Users/davidcampbell633/Desktop/Sym_Homo.pdf')


#rho_ss = states[Npts-1]
#print( q.expect( S12_S34*S12_S34.dag() , rho_ss )  ) 
#ello = S13_S24 + S14_S23
#ello = ello.unit()
#print( q.expect( ello*ello.dag() , rho_ss )  ) 

#q.matrix_histogram(rho_ss)
