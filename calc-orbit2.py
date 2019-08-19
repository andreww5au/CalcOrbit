#!/usr/bin/env python

"""Solar system simulator, written by Dr Andrew Williams <andrew@longtable.org>.
   
   Uses:
     planephem.py    calculates initial 3D coordinates and velocities for planets and random asteroids, and
     bruteforce2.py  performs the actual N-body modelling via a call to 'Process'
     jupsurf.tga     A texture map for Jupiter to show surface features.

   Released under the GPL, feel free to use it for non-commercial purposes. I'd appreciate an email with
   what you're using it for, and any nifty improvements you've made.
"""

from copy import copy
import time
import traceback

import vpython
from vpython import vector as v

import planephem
import bruteforce2 as nbody

################################# Configurable Parameters for the model ###############################

AU = 1.49597870700e11  # Astronomical Unit (mean Earth-Sun distance) in metres.
t = 2451545.0  # Start modelling at this date (2451545.0 = J2000.0)
dt = 1.0  # Calculation time step, in days (default 1.0)
Substeps = 10  # Number of calculated position updates between each display update
PreCalc = 0  # Number of years to pre-calculate before starting visual display
PrInterval = 10  # Number of days between displaying timing info

dispmult = 0.1  # The radius of each sphere displayed is this constant
#    multiplied by the display radius in the object parameters
capmult = 5  # Multiply the physical body radius by this amount to get the capture radius

# PlanetsToModel = ['Sun', 'Earth', 'Mars', 'Jupiter']

PlanetsToModel = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']

################ Configure the properties of the 'Asteroids' (massless test particles) ############

Nasteroids = 1000  # Number of random asteroids

RandomA = False  # If false, the semi-major axes for all asteroids are evenly distributed, instead of random
rpos = 12  # Start angle for orbital position, in degrees ('None' for random positions)

amin, amax = 3.0, 5.5  # bounds for semi-major axis distribution
emin, emax = 0.0, 0.00  # Bounds for eccentricity variation (zero for circular orbits)
Imin, Imax = -0.0, 0.0  # Bounds for inclination angle variation, in degrees (zero for orbits in the plane of the ecliptic)

############################ Set up display window and lighting ####################################


win = 800  # Display window size, in pixels
scene = vpython.canvas(title="Orbit", width=win, height=win)  # Create window
scene.lights = []  # Dump default scene lighting
scene.ambient = vpython.color.gray(0.3)  # Low level ambient lighting
sunlight = vpython.local_light(pos=v(0, 0, 0), color=vpython.color.white)  # Sunlight from centre


########################### Class and function definitions #############################################
class Body(planephem.Planet):
    """Subclass the 'Planet' object in the ephemeris library to add methods for setting up the
       3D Visual Python objects to display it.
    """

    def __init__(self, *args, **keywords):
        planephem.Planet.__init__(self, *args, **keywords)
        self.color = vpython.color.white
        self.dradius = 1.0
        self.texture = None
        self.rings = False
        self.up = (0, 0, 1)
        self.body = None

    def setbody(self, color=None, dradius=None, rings=None, simple=False, up=None, texture=None):
        """Define the sphere representing the object itself. The radius is given in units of
           'dispmult', a global parameter allowing human-friendly sizes. Store initial
           position, mass, and momentum from the ephemeris library.

           if simple is True, use a simple_sphere object with fewer vertices, to render faster.
        """
        if color is None:
            color = self.color
        else:
            self.color = color
        if dradius is None:
            dradius = self.dradius
        else:
            self.dradius = dradius
        if texture is None:
            texture = self.texture
        else:
            self.texture = texture
        if rings is None:
            rings = self.rings
        else:
            self.rings = rings
        if up is None:
            up = self.up
        else:
            self.up = up

        if rings:
            # self.body = vpython.frame(pos=v(*self.pos))
            rings = vpython.cylinder(pos=v(0, 0, 0),
                                     radius=dradius * dispmult * 3,
                                     axis=v(0, 1, 0),             # (0, 1, 2)
                                     length=0.1 * dispmult,
                                     color=color,
                                     texture=texture,
                                     up=v(*up))
            if simple:
                planet = vpython.simple_sphere(pos=v(0, 0, 0),
                                               radius=dradius * dispmult,
                                               color=color,
                                               up=v(*up))
            else:
                planet = vpython.sphere(pos=v(0, 0, 0),
                                        radius=dradius * dispmult,
                                        color=color,
                                        texture=texture,
                                        up=v(*up))
            self.body = vpython.compound([rings, planet])
            self.body.pos = v(*self.pos)
            self.body.axis = v(0, 1, 2)
        else:
            if simple:
                self.body = vpython.simple_sphere(pos=v(*self.pos),
                                                  radius=dradius * dispmult,
                                                  color=color,
                                                  up=v(*up))
            else:
                self.body = vpython.sphere(pos=v(*self.pos),
                                           radius=dradius * dispmult,
                                           color=color,
                                           texture=texture,
                                           up=v(*up))
        self.body.cradius = self.rad * capmult
        self.body.mass = self.M


def ShowErrors():
    """For each of Sun, Mercury,...,Ceres, recalculate the position for the current time using the
       ephemeris library, and display an arrow indicating the difference between the current modelled
       position and the calculated ephemeris position. The errors will gradually increase, because of
       the non-real effects of any random asteroids, because many real bodies are not in the model
       (real asteroids/comets/moons, nearby stars, etc), and because the initial velocity calculation
       is fairly crude, much less accurate than the position.
    """
    global arrows
    sunpos = copy(bodies['Sun'].body.pos)
    for ar in arrows:
        ar.visible = False
    arrows = []
    for name in Planets:
        if name in bodies.keys():
            bod = bodies[name]
            bod.calc(Teph=t)
            curpos = v(copy(bod.body.pos))
            errvec = v(*bod.pos) + sunpos - curpos
            if bod.body:
                arrows.append(vpython.arrow(pos=curpos, axis=errvec, shaftwidth=0.5 * dispmult))


def scalesizes(factor):
    """Scale the displayed sizes of all bodies up or down by the provided factor
    """
    # TODO - simple_sphere() objects don't change size when their radius attribute is changed. Need
    #        to deleted every asteroid simple_sphere() and and create new ones, with the new size.
    for b in bodies.values():
        b.body.radius *= factor


def ProcessKeys(event):
    """Handle any keystrokes during the model.
    """
    global paused, Tracking, scenecenter, TrackAngle, TrackMag, scene, tlabel, Trail
    try:  # Key pressed:
        k = event.key
        if (k == '0') and 'Sun' in PlanetsToModel:  # Keys 0-9 keep the view centered on Sun-Pluto respectively
            scenecenter = bodies['Sun']
        elif (k == '1') and 'Mercury' in PlanetsToModel:
            scenecenter = bodies['Mercury']
        elif (k == '2') and 'Venus' in PlanetsToModel:
            scenecenter = bodies['Venus']
        elif (k == '3') and 'Earth' in PlanetsToModel:
            scenecenter = bodies['Earth']
        elif (k == '4') and 'Mars' in PlanetsToModel:
            scenecenter = bodies['Mars']
        elif (k == '5') and 'Jupiter' in PlanetsToModel:
            scenecenter = bodies['Jupiter']
        elif (k == '6') and 'Saturn' in PlanetsToModel:
            scenecenter = bodies['Saturn']
        elif (k == '7') and 'Uranus' in PlanetsToModel:
            scenecenter = bodies['Uranus']
        elif (k == '8') and 'Neptune' in PlanetsToModel:
            scenecenter = bodies['Neptune']
        elif (k == '9') and 'Pluto' in PlanetsToModel:
            scenecenter = bodies['Pluto']
        elif (k == 't'):  # If tracking is on, keep the same viewpoint when following an object
            print("'t' pressed, Tracking=%s, TrackAngle=%s" % (Tracking, TrackAngle))
            if not Tracking:
                bvec = v(*nbody.hp[scenecenter.i]).norm()
                TrackAngle = scene.camera.axis.norm() - bvec
                TrackMag = scene.camera.axis.mag
                Tracking = True
                print("Tracking")
            else:
                Tracking = False
                print("Not Tracking")
        elif k == 'u':
            if scene.up == v(0, 0, 1):
                scene.up = v(0, 1, 0)
            else:
                scene.up = v(0, 0, 1)
        elif k == 'b':  # Up and down arrows scale the displayed object sizes up and down
            scalesizes(1.2)
        elif k == 's':
            scalesizes(0.8)
        elif (k == 'p'):  # Toggle the display of orbital paths and dots behind objects
            if paused:
                paused = False
                print('Resumed')
            else:
                paused = True
                print('Paused')
        elif (k == 'c'):
            if Trail is not None:
                Trail.clear()
                Trail.stop()
                Trail = None
        elif (k == 'e'):  # Show differences between modelled and ephemeris positions
            ShowErrors()
        if k in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:  # New scene center object, update display window
            scene.center = scenecenter.body.pos
            tlabel.pos = scene.center
    except:   # Mouse click
        print("Exception in ProcessKeys: %s" % traceback.format_exc())
        pass


def ProcessClicks(evt):
    global Trail
    b = scene.mouse.pick
    if b is not None:  # If we've clicked on an actual object, and not the background
        try:
            if Trail is not None:
                Trail.clear()
                Trail.stop()
            print('Trailing: %s' % b)
            Trail = vpython.attach_trail(b, type='points', radius=0.02, pps=10, retain=600)
        except:
            pass  # Can't leave a trail for Saturn, or a breadcrumb point
    else:
        Trail.clear()
        Trail.stop()
        Trail = None


############################### Define and set up the bodies in the model ###############################

# juptex = vpython.materials.texture(data=vpython.materials.loadTGA('jupsurf.tga'), mapping='spherical')

allplanets = planephem.GetPlanets(oclass=Body, Teph=t)  # Get Sun, planets, Pluto and Ceres in a python dict

# Create Visual Python display objects for each of the planets we want to model
planets = {}
for name in PlanetsToModel:
    planets[name] = allplanets[name]
    if name == 'Sun':
        planets[name].setbody(color=v(1.0, 1.0, 0.5), dradius=3)  # , material=vpython.materials.emissive)
    elif name == 'Mercury':
        planets[name].setbody(color=vpython.color.gray(0.5), dradius=0.6)
    elif name == 'Venus':
        planets[name].setbody(color=vpython.color.white, dradius=1)
    elif name == 'Earth':
        planets[name].setbody(color=v(0, 0.5, 1.0), dradius=1)
    elif name == 'Mars':
        planets[name].setbody(color=v(1.0, 0.2, 0.2), dradius=1)
    elif name == 'Jupiter':
        planets[name].setbody(color=vpython.color.orange, dradius=3, texture={'file':'jupsurf.png', 'mapping':'spherical'})
    elif name == 'Saturn':
        planets[name].setbody(color=vpython.color.orange, dradius=1, rings=True)
    elif name == 'Uranus':
        planets[name].setbody(color=vpython.color.white, dradius=2)
    elif name == 'Neptune':
        planets[name].setbody(color=vpython.color.blue, dradius=2)
    elif name == 'Pluto':
        planets[name].setbody(color=vpython.color.white, dradius=0.8)
    elif name == 'Ceres':
        planets[name].setbody(color=vpython.color.green, dradius=0.8)

# Add planets to the nbody simulation engine as heavy bodies (the masses of these bodies affect all objects)
for name, obj in planets.items():
    obj.i = nbody.AddHeavyBody(name=name, body=obj.body, pos=obj.pos, vel=obj.vel, M=obj.M, r=obj.rad * capmult)

# Get asteroids

if RandomA:
    asteroids = planephem.GetRandomAsteroids(oclass=Body, Teph=t,
                                             amin=amin, amax=amax,
                                             emin=emin, emax=emax,
                                             Imin=Imin, Imax=Imax,
                                             rpos=rpos,
                                             N=Nasteroids)
else:
    asteroids = planephem.GetEvenAsteroids(oclass=Body, Teph=t,
                                           amin=amin, amax=amax,
                                           emin=emin, emax=emax,
                                           Imin=Imin, Imax=Imax,
                                           rpos=rpos,
                                           N=Nasteroids)

# Define visual python display objects and orbits for them
for mp in asteroids.values():
    mp.setbody(color=vpython.color.white, dradius=0.5, simple=True)

# Add asteroids to the nbody simulation engine as light bodies (modelled as massless particles)
for name, obj in asteroids.items():
    obj.i = nbody.AddLightBody(name=name, body=obj.body, pos=obj.pos, vel=obj.vel)

# create one dict called 'bodies' containing all objects - planets and asteroids
bodies = planets
bodies.update(asteroids)

############################### Initialise the Nbody code ready for processing ########################
nbody.Init(sun=bodies['Sun'].i)

############################## Set up initial variables before starting model #########################
startt = t
paused = False
scenecenter = bodies['Sun']
scene.autoscale = False
scene.autocenter = False
scene.show_rendertime = True
scene.bind('keydown', ProcessKeys)
scene.bind('click', ProcessClicks)
arrows = []
Tracking = False
TrackAngle = None
TrackMag = None
Trail = None
lasttime = time.time()
timingstring = ""
tlabel = vpython.label(pos=v(0, 0, 0), text='', xoffset=-150, yoffset=450, box=False, line=False)
Planets = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Ceres']

################## Optionally precalculate model for some time before starting display ################
i = PreCalc * 365.25 / (dt * Substeps)
while i > 0:
    nbody.Process(t, steps=Substeps, dt=dt)
    if i % 1000 == 0:
        print("Steps remaining: %d" % i)
    t += dt * Substeps
    i -= 1
nbody.UpdatePositions(vclass=v)


########################################### Main loop #################################################


def mainloop():
    global t, tlabel, lasttime, timingstring
    if int(startt - t) % PrInterval == 0:  # Occasionally update timing information
        timingstring = "%4.3f seconds per %d days" % (time.time() - lasttime, PrInterval)
        lasttime = time.time()
        # print(timingstring)

    # Use the Nbody engine to run the model for 'Substeps' steps, each of time 'dt'
    nbody.Process(t, steps=Substeps, dt=dt)
    t = t + dt * Substeps
    tlabel.text = "Year=%6.4f (%d JD). Timing: %s" % ((t - 2451545.0) / 365.246 + 2000.0, t, timingstring)

    # Update positions of display objects
    nbody.UpdatePositions(vclass=v)

    # Update display window viewing center and viewing angle
    if scenecenter != bodies['Sun']:
        scene.center = scenecenter.body.pos
        tlabel.pos = scene.center
        if Tracking:
            bvec = v(*nbody.hp[scenecenter.i]).norm()
            scene.camera.axis = (TrackAngle + bvec) * TrackMag
            scene.camera.pos = scenecenter.body.pos - scene.camera.axis


while 1:
    vpython.rate(50)
    if not paused:
        mainloop()
