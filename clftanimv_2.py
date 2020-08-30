"""
Simulation for ring of N classical harmonic oscillators
Note, this program as written assumes k=1, m=1 for each oscillator
"""

from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all')

def delta(i1,i2): #define Kronecker delta
    
    if i1 == i2:
        
        return 1
    
    else:
        
        return 0

def Theta(x): #define Heaviside step function
    
    if x > 0:
        
        return 1
    
    else:
        
        return 0
    
def fn(n, N, a, b): #Function describing the initial configuration of compressed/extended springs
    
    x = 2*n*np.pi/N
    
    func = b*(Theta(a-(x+a/2))*Theta(x+a/2)*np.sin(np.pi*(x+a/2)/a)**2
                +Theta((x-a/2)-(2*np.pi-a))*Theta(2*np.pi-(x-a/2))*np.sin(np.pi*((x-a/2)-2*np.pi)/a)**2)
    
    return func
    
def kmat(N): #generates NxN matrix whose eigenvalues give normal mode frequencies
    
    K = [[((N/(2*np.pi))**2)*(2*delta(i,j)-delta((i+1) % N, j)-delta(i, (j+1) % N)) 
          for i in range(N)] for j in range(N)]
    
    #Here, the factor of (N/(2*pi))^2 ensures that the speed of sound does not
    #change with the number of oscillators
    
    return K

def eigen(array):
    
    w,v = LA.eig(array)
    
    return [w,v]

def qn(N, n, t, w, A, alph, beta):
    
    # qn(t) = \sum_{j=0}^{N-1} alph[n]A[n,j]cos(w[j]t) + beta[n]A[n,j]sin(w[j]t)
    
    lmat = []
    
    i = 0
    
    while i < N:
        
        lmat.append(A[n,i]*alph[i]*np.cos(w[i]*t) 
                    + A[n,i]*beta[i]*np.sin(w[i]*t))
        
        i += 1
        
    return sum(lmat)

def qvec(t, N): #vector of positions of each oscillator at given time

    d0 = np.pi/(2*N)
    
    a = 1/6

    K = kmat(N)

    esyst = eigen(K)

    wsq = esyst[0]
    A = esyst[1]
    
    n = 0
    
    while n < N: #This step ensures that all frequencies which are effectively zero
                 #remain very small after taking square root in next step. Note, the way
                 #that LA solves systems of equations means that if we just take these
                 #frequencies to be zero, it will result in a singular matrix which
                 #cannot be solved.
        wsq[n] = abs(wsq[n])
        
        if wsq[n] < 1e-15:
            wsq[n] = wsq[n]**2
        
        n += 1
        
    w = np.sqrt(wsq)

    q0 = []
    p0 = []
    
    if N < 100:

        q0.append(d0)
        p0.append(0)
        
        n = 1
        
        while n < N:
            
            q0.append(0)
            p0.append(0)
            
            n += 1

    else:
        
        n = 0

        while n < N:
            
            q0.append(fn(n, N, a, d0))
            p0.append(0)
    
            n += 1

    asol = LA.solve(A, q0)

    B = [[w[i]*A[j,i] for i in range(N)] for j in range(N)]

    bsol = LA.solve(B, p0)

    q = []
    
    n = 0
    
    while n < N:
        q.append(qn(N, n, t, w, A, asol, bsol))
        
        n += 1
        
    return q

def distance(N, t):
    dist = []
    
    n = 0
    
    sol = qvec(t, N)
    
    while n < N:
        
        dist.append(sol[(n+1) % N] - sol[n])
        
        n += 1
        
    return dist

def draw_figs(N, t, label):

    sols = qvec(t, N)

    xpos = []

    ypos = []
    
    dist = distance(N, t)

    n = 0

    while n < N:
        xpos.append(np.cos(2*np.pi*n/N + sols[n]))
        ypos.append(np.sin(2*np.pi*n/N + sols[n]))
    
        n += 1
    
    xpos.append(np.cos(sols[0]))
    ypos.append(np.sin(sols[0]))
    
    rcol = []
    bcol = []
    
    n = 0
    
    while n < N:
        
        if dist[n] > 0:
            
            bcol.append(dist[n])
            rcol.append(0)
            
        else:
        
            rcol.append(-dist[n])
            bcol.append(0)
        
        n += 1
        
    bnorm = max(bcol)
    rnorm = max(rcol)
    
    # print(rcol, bcol)
    
    fig, ax = plt.subplots()
    
    plt.text(0,0,'N={}'.format(N),ha='center')
    
    # ax = fig.add_axes([0.0,0.0,1.0,1.0])
    
    ax.axis('off')
    ax.set_xlim([-1.25,1.25])
    ax.set_ylim([-1.25,1.25])
    ax.set_aspect('equal')
    
    n = 0
    
    while n < N:
        
        if N < 100:
        
            plt.plot([xpos[n], xpos[n+1]], [ypos[n], ypos[n+1]], 
                  color = (1/2*(1+rcol[n]/rnorm-bcol[n]/bnorm),
                           1/2*(1-rcol[n]/rnorm-bcol[n]/bnorm),
                           1/2*(1-rcol[n]/rnorm+bcol[n]/bnorm)))
            
        else:
            
            plt.plot([xpos[n], xpos[n+1]], [ypos[n], ypos[n+1]], 
                  color = (1-bcol[n]/bnorm, 1-rcol[n]/rnorm-bcol[n]/bnorm, 1-rcol[n]/rnorm))
        
        n += 1
        
    rad = 1/(2*N)
    
    if N < 100:
        
        i = 0
        
        while i < N:
    
            circle = plt.Circle((xpos[i],ypos[i]),rad,
                          color = 'black')
            ax.add_artist(circle)
            
            i += 1
        
    filename = 'img{}'.format(label)
    
    if os.path.isdir('clftN{0}'.format(N)) == False:
        os.mkdir('clftN{0}'.format(N))
        
    if os.path.isfile('clftN{0}/{1}'.format(N,filename)) == True:
        os.remove('clftN{0}/{1}'.format(N,filename))
    
    fig.savefig('clftN{0}/{1}'.format(N,filename), dpi = 500)
    
    plt.close(fig)
    
plt.ioff()

N = 11 #Number of oscillators

tmax = 20 #max time to run simulation

t = 0 #starting time

label = 1 #starting label for output images

while t < tmax:
    
    draw_figs(N, t, label)
    t += 0.1
    label += 1