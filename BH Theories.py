#!/usr/bin/env python
# coding: utf-8

# # Black Hole Spin Comparisons

# In[1]:


from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


# ### Constants

# In[2]:


G=6.67*10**-11
c=299792458
Msun=1.989*10**30


# ### G-Mode

# In[3]:


def FunGMode():
    r=np.arange(2.5, 11.83, 0.1,  dtype=np.complex)
    #r=2.5
    a2=(r**.5)*(5*r-4)/2
    a1=r*(1-6*r)
    a0=r**(2.5)*(8-r)/2
    Q=1/3*a1-1/9*a2**2
    R=1/6*(a1*a2-3*a0)-(1/27)*a2**3
    check=Q**3+R**2
    ##(Q**3+r**2)<0 so all roots are real
    S1=(R+(check)**(1/2))**(1/3)
    S2=(R-(check)**(1/2))**(1/3)
    Z1=(S1+S2)-(a2/3)
    Z2=(-1/2)*(S1+S2)-(a2/3)+(3**(1/2)*1j/2)*(S1-S2)
    a=(-1/2)*(S1+S2)-(a2/3)-(3**(1/2)*1j/2)*(S1-S2)
    ## Z3 is correct root

    alphar=1-6*r**-1+8*a*r**(-3/2)-3*a**2*r**(-2)

    gmode=(alphar**(1/2)*(c**3/(G*Msun))*(r**(3/2)+a)**-1*(1/(20*np.pi)))

    #0<L/Ledd<.5
    L=.5
    E=.1*L
    lowerbound=gmode*(1-E)
    
    return gmode, lowerbound, a


# In[4]:


def ZeroGMode():
    r=11.83
    a2=(r**.5)*(5*r-4)/2
    a1=r*(1-6*r)
    a0=r**(2.5)*(8-r)/2
    Q=1/3*a1-1/9*a2**2
    R=1/6*(a1*a2-3*a0)-(1/27)*a2**3
    check=Q**3+R**2
    ##(Q**3+r**2)<0 so all roots are real
    S1=(R+(check)**(1/2))**(1/3)
    S2=(R-(check)**(1/2))**(1/3)
    Z1=(S1+S2)-(a2/3)
    Z2=(-1/2)*(S1+S2)-(a2/3)+(3**(1/2)*1j/2)*(S1-S2)
    a=(-1/2)*(S1+S2)-(a2/3)-(3**(1/2)*1j/2)*(S1-S2)
    ## Z3 is correct root

    alphar=1-6*r**-1+8*a*r**(-3/2)-3*a**2*r**(-2)

    gmode=(alphar**(1/2)*(c**3/(G*Msun))*(r**(3/2)+a)**-1*(1/(20*np.pi)))

    #0<L/Ledd<.5
    L=.5
    E=.1*L
    lowerbound=gmode*(1-E)
    return gmode, lowerbound, a


# In[15]:


# plot gmode
gmode, lowerbound, a = FunGMode()
plt.fill_between(a, lowerbound, gmode, facecolor='turquoise')
plt.plot(a, gmode, 'turquoise')
plt.plot(a, lowerbound, 'turquoise')

g, low, aa = ZeroGMode()
plt.plot(aa, g, 'bo')
plt.plot(aa, low, 'bo')

print("for a = ", aa, ", upperbound =  ", g, "and lowerbound = ", low)


# ### 3:2 Resonance Theory

# In[6]:


def Fun32():
    r=np.arange(49/5-4*(11/5)**(1/2), (49/5) + 4*(11/5)**(1/2), .1)
    a=(1/39)*(44*(r)**(1/2)-(5)**(.5)*(r*(39*r-34))**(1/2))

    alphar=1-6*r**-1+8*a*r**(-3/2)-3*a**2*r**(-2)
    alphath=1-4*a*r**(-3/2)+3*a**2*r**(-2)

    radial=(alphar**(1/2)*(c**3/(G*Msun))*(r**(3/2)+a)**-1*(1/(20*np.pi)))

    vertical=((alphath)**(1/2)*(c**3/(G*Msun))*(r**(3/2)+a)**-1*(1/(20*np.pi)))
    
    return radial, vertical, a


# In[7]:


def zero32():
    r=(49/5) + 4*(11/5)**(1/2)
    a=(1/39)*(44*(r)**(1/2)-(5)**(.5)*(r*(39*r-34))**(1/2))

    alphar=1-6*r**-1+8*a*r**(-3/2)-3*a**2*r**(-2)
    alphath=1-4*a*r**(-3/2)+3*a**2*r**(-2)

    radial=(alphar**(1/2)*(c**3/(G*Msun))*(r**(3/2)+a)**-1*(1/(20*np.pi)))

    vertical=((alphath)**(1/2)*(c**3/(G*Msun))*(r**(3/2)+a)**-1*(1/(20*np.pi)))
    
    return radial, vertical, a


# In[8]:


rad, vert, aa = zero32()
print("for a = ", aa, ", vertical = ", vert, "and rad = ", rad)

radial, vertical, a = Fun32()
plt.plot(a, radial, 'gold')
plt.plot(a, vertical, 'gold', ls='--')
plt.plot(aa, rad, 'bo')
plt.plot(aa, vert, 'bo')


# ### ISCO

# In[9]:


def FunISCO():
    r=np.arange(1, 9, 0.1)
    a=((-1/6)*r**2*(-8/r**(3/2)+(64/(r**3)+(12*(1-6/r))/(r**2))**(.5)))
    ISCO=((r**(3/2)+a)**(-1)*(1/(2*np.pi))*(c**3/(10*Msun*G)))
    return ISCO, a


# In[10]:


def zeroISCO():
    r=9
    a=((-1/6)*r**2*(-8/r**(3/2)+(64/(r**3)+(12*(1-6/r))/(r**2))**(.5)))
    ISCO=((r**(3/2)+a)**(-1)*(1/(2*np.pi))*(c**3/(10*Msun*G)))
    return ISCO, a


# In[11]:


zeroISCO, aa = zeroISCO()
print("for a = ", aa, ", ISCO = ", zeroISCO)

ISCO, a = FunISCO()
plt.plot(a, ISCO, 'k')
plt.plot(aa, zeroISCO, 'bo')


# ### BH Angular Frequency

# In[12]:


def FunAngFreq():
    a=np.arange(-1, .998, 0.01)
    angfreq=(c*a)/((Msun)*(G/c**2)*(160*np.pi*(1+(1-a**2)**.5)))
    return angfreq, a


# In[13]:


def ZeroAngFreq():
    a=-1
    angfreq=(c*a)/((Msun)*(G/c**2)*(160*np.pi*(1+(1-a**2)**.5)))
    return angfreq, a


# In[14]:


zeroang, aa = ZeroAngFreq()
print("for a = ", aa, ", BH angular frequency = ", zeroang)

angfreq, a = FunAngFreq()
plt.plot(a, angfreq, 'k--')
plt.plot(aa, zeroang, 'bo')

