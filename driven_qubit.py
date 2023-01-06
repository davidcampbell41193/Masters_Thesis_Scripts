from qutip import *
import numpy as np
import matplotlib.pyplot as plt

gamma=1.0 ;
#omega=9.0 ; # Driving frequency
nu_prime=0.0 ; # Qubit frequency detuning
Omega=20.0 ; # Drive strength
N= 1.0; # Average number of excitations in the bath 

psi0=fock_dm(2,0) ;   
a=sigmam() ;

times=np.linspace(0,3,1000)          

H0=nu_prime/2*sigmaz() ;
H1=Omega/2*sigmax() ;

def H1_coeff(times, args):
    return np.cos(omega*times)

H=H0+H1 ;

e = fock_dm(2,0)
g = fock_dm(2,1)

output=mesolve(H,psi0,times,[np.sqrt(gamma*(N+1))*a, np.sqrt(gamma*N)*a.dag()], [ e*e.dag() , g*g.dag() , q.sigmaz() , q.sigmax() , q.sigmay() ]) ;


    
left  = 0.125  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.2 
              
plt.subplots_adjust(left=3, bottom=None, right=1000, top=None, wspace=None, hspace=None)

fig=plt.figure()
fig, (ax1)=plt.subplots(1,1) ;
#plt.title(r'$\gamma=$%d $\epsilon=$%d $\Omega=$%d $N=$%d' %(gamma, nu_prime, Omega , N ) , loc='left')
 
my_suptitle=fig.suptitle(r'$\gamma=$%d $\epsilon=$%d $\Omega=$%d $N=$%d' %(gamma, nu_prime, Omega , N ) , verticalalignment='bottom', fontsize=18 )
                
ax1.plot(gamma*output.times, np.transpose(output.expect[0]), color='r' , label=r'$|e\rangle$') ;
ax1.plot(gamma*output.times, np.transpose(output.expect[1]), color='b' , label=r'$|g\rangle$') ;
ax1.set_xlabel(r'$\gamma t$') ;
ax1.set_ylabel('Fidelity') ;
ax1.legend() ; 


ax2.plot(gamma*output.times, np.transpose(output.expect[2]) , color='g' , label=r'$\langle \sigma_z \rangle$') ;
ax2.plot(gamma*output.times, np.transpose(output.expect[3]) , color='c' , label=r'$\langle \sigma_x \rangle$') ;
ax2.plot(gamma*output.times, np.transpose(output.expect[4]) , color='m' , label=r'$\langle \sigma_y \rangle$') ;
ax2.set_xlabel(r'$\gamma t$') ;
ax2.set_ylabel('Expectation Value') ;
ax2.legend()
plt.tight_layout( w_pad=3.0  )
fig.suptitle(r'$\gamma=$%d $\epsilon=$%d $\Omega=$%d $N=$%d' %(gamma, nu_prime, Omega , N ) , verticalalignment='bottom', fontsize=18 )
plt.savefig('/Users/davidcampbell633/Desktop/decoherence.pdf', dpi=fig.dpi, bbox_inches='tight',bbox_extra_artists=[my_suptitle])

