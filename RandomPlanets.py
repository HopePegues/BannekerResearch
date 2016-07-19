import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import interpolate

StarName, Distance, Mass, Radius, Period, Ewha, Lhalbol = np.loadtxt('BigData.cat',delimiter=',', unpack=True)

def RandomPlanet(num):    
    #INPUT: num = number of stars in our list
    #Randomly assign 0 or 1 to each star of my list, 1 representing a planet's prescence
    #OUTPUT: stwpl, PLMass, PLPeriod
    #stwpl = array with the indexes of stars that have planets = SelectedStars
    #PLMass = randomly assigned mass of of each planet
    #PLTPeriod = Period of each planets
 
    plts = np.random.randint(0,high=2,size = len(num))
    stwpl, = np.where(plts==1)
    PLMass=(np.random.random(size=len(stwpl)))*1.5 + 0.5
    PLPeriod = 15.
    return stwpl, PLMass,PLPeriod

SelectedStars,PLMass,PLPeriod = RandomPlanet(StarName)
Mass_SS = Mass[SelectedStars]

def Rossby(rotational_period,mass):
    #Calulates the Rossby Numbers for all of the stars in our list
    Mo = 1.
    logtau = 1.16 - 1.49 * np.log10(mass/Mo) - 0.54 * (np.log10(mass/Mo))**2.
    tau = 10.**logtau
    Ro = rotational_period/tau
    return Ro

Ro = Rossby(Period,Mass)

def RVSignal(star_mass,planet_mass,planet_period):
    #Calculates the radial velocities of our stars with planets
    G = 6.67*10.**-11.
    M_e = 5.972*10.**24.
    M_s = 1.989*10.**30.
    sec = 86400.
    a = (2.*math.pi)**(1./3.)*((G**(1./3.))*M_e)/((sec**(1./3.))*M_s**(2./3.))
    vs = a*((planet_mass**3.)/(planet_period*star_mass**2.))**(1./3.)
    return vs

starv = RVSignal(Mass_SS,PLMass,PLPeriod)

def Likelihood(P):
	#Generates our likelihood function for our planets' periods
	
	if P > 10.:
		return 1.
	if 1. < P < 10.:
		return P/10. #For now
	if P < 1.:
		return 0.

def PlanetPeriod(num):
	#Randomly generate a planet within the bounds of its likelihood from Likelihood function
	#Return the period of that planet if it is under the curve of its likelihood

	period = 10.**np.random.uniform(-1,2,num)
	logP = np.array([Likelihood(i) for i in period])
	y = np.random.uniform(0,1,num)
	#print logP,y
	return y[y < logP],period[y<logP]

def InverseTransform(CProb):
	c = 2./3.
	if CProb < c/2.:
		return np.sqrt(2./c*CProb)
	else:
		return (CProb+c/2.)/c
		
def PlanetPeriod2(num):
	y = np.random.uniform(0,1,num)
	logperiod = np.array([InverseTransform(i) for i in y])
	return 10.**logperiod

#Using Morton&Swift's data to get the cdf of planet radii
xdata,ydata,x3,x4 = np.loadtxt('FitThisData.txt', unpack=True)
cml = np.zeros_like(xdata)
for i,ix in enumerate(xdata):
    hi = np.trapz(ydata[:i+1],xdata[:i+1])
    cml[i] = hi
cml = cml/np.trapz(ydata,xdata)

spline = interpolate.splrep(cml,xdata)
rad = interpolate.splev(np.random.uniform(0,1,10000),spline)
	
periods=PlanetPeriod2(20000)
plt.hist(np.log10(periods),50)
value, periods = PlanetPeriod(20000)
plt.hist(np.log10(periods),50, alpha=0.5)
plt.show()
plt.hist(rad)
plt.show()

plt.hist(starv)
plt.show()