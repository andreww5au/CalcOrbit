# CalcOrbit
3D model of the Solar System, written using VPython

This solar system simulator calculates an initial position and velocity vector for every planet in 3D cartesian 
coordinates using an approximate ephemeris and algorithm from JPL. It also calculates positions and velocities 
for an arbitrary number of additional bodies defined by ranges of orbital parameters (semi major axis, eccentricity, 
etc), with all coordinates calculated for J2000.0 (January 1st, 2000, 12:00 UTC). The library that calculates 
the initial positions and velocities is 'planephem.py'.

Once the initial state is defined, an N-body gravitational model is started, with a default timestep of one day. 
For each timestep, the gravitational force exerted on every body, by every other body, is calculated, and the 
momenta and positions of each body is updated. To save time, the asteroids are modelled as massless test particles,
so while the planets extert forces on the asteroids, the reverse forces are ignored. All of the N-body processing 
code is in bruteforce2.py.

This N-body model is just a toy, as far as real N-body simulations are concerned - real N-body code is highly 
optimised, and in a low level language, not Python. This model is designed to teach the concepts, and keep the 
model visible, in real time.

The calc-orbit.py program recognises the following key presses:

0-9: Changes the viewpoint centre to either the Sun (0) or Mercury-Pluto (1-9), and follow that object, if 
it's moving. 

t: Toggle 'tracking mode', which automatically rotates the camera to keep the viewing angle constant when 
following a moving planet.

u: Toggle the 'up' direction used to manage camera rotation in 'tracking mode' when following a moving planet.

p: Pause the animation, or resume if it's already paused.

c: Clear the current dot-trail.

e: Display a set of arrows showing calculation error vectors - the difference between the true (ephemeris) 
position of a planet on that date, and the calculated position in the N-body simulation.

b/s: Make the displayed size of the planets and asteroids bigger or smaller, respectively. The 
displayed sizes are purely for display purposes, and are unrelated to the real physical sizes, which are 
only used to calculate they capture radius for collisions. NOTE - only the planet sizes
are changed - I think it's because vpython simple_sphere() objects can't have their radius
changed on the fly. 

Mouseclick: Click on an object to make it start leaving a trail. Click on empty space
to turn off the trail. It's easier to click on an object if you pause the animation 
(press 'p') first, click, then unpause (press 'p' again).


Changing the model parameters:
All parameters are defined in the code itself, in calc-orbit2.py. Configurable lines of code include:

################################# Configurable Parameters for the model ###############################

AU = 1.49597870700e11 #Gravitational constant
t = 2451545.0         #Start modelling at this date (2451545.0 = J2000.0)
dt = 1.0              #Calculation time step, in days (default 1.0)
Substeps = 5          #Number of calculated position updates between each display update
PreCalc = 0           #Number of days to pre-calculate before starting visual display
PrInterval = 1000     #Number of days between displaying timing info
DotInterval = 10      #Number of days between displaying a dot for asteroid trails
dispmult = 0.1        #The radius of each sphere displayed is this constant 
                      # multiplied by the display radius in the object parameters
capmult = 5           #Multiply the physical body radius by this amount to get the capture radius
ShowPath = False      #Display dots/paths as orbits. Can be toggled during runtime with 'p' key.
PlanetsToModel = ['Sun', 'Earth', 'Mars', 'Jupiter']
################ Configure the properties of the 'Asteroids' (massless test particles) ############
Nasteroids = 500      #Number of random asteroids
RandomA = False       #If false, the semi-major axes for all asteroids are evenly distributed, instead of random
rpos = None            #Start angle for orbital position, in degrees ('None' for random positions)
amin, amax = 2.45, 2.55  #bounds for semi-major axis distribution
emin, emax = -0.0, 0.0   #Bounds for eccentricity variation (zero for circular orbits)
Imin, Imax = 0.0, 0.0    #Bounds for inclination angle variation, in degrees (zero for orbits in the plane of the ecliptic)

The important parameters for asteroid orbits are amin, amax, which define the minimum and maximum value for 
semi-major axis, emin,emax which define the minimum and maximum eccentricity, and Imin,Imax which define the
min and max orbital inclination angle. All values are randomly chosen for each of the Nasteroids generated, 
within the ranges specified. The only exception is that if 'RandomA' is False, instead of a random choice of 
semi-major axis, a non-random even spread of semi-major axes are chosen distributed evenly over the specified range.

The 'rpos' parameter can either be a value in degrees (0-360) defining the starting point in it's orbit for all 
asteroids, or 'None' for random starting positions. Using a fixed value of 'rpos' and fixed 'e' and 'I' values 
(max and min values equal) gives a population of asteroids that vary only in semi-major axis, and allow very 
strong visible evidence of minor perturbations due to planets.

CAVEATS:
This code was written a long time ago - if I was writing it now, I'd just use skyfield (https://rhodesmill.org/skyfield/) 
or AstroPy (https://www.astropy.org/) to calculate the initial solar system positions instead of writing code to use the
JPL ephemeris. The code was initially written using 'classic' VPython (version 6 or lower), so some of the features
may not work well (or at all) using VPython 7.
