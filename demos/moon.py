# Demonstration of tidal locking in the Moon's orbit
#
# There are three complete rotations of the planet for every two orbits of
# the sun.
#
# Hopelessly out of proper scale, of course, it's meant to illustrate the motion
#
# Andrew Williams, 2002

import vpython
from vpython import vector as v
import math
import time

win = 700
scene = vpython.canvas(title="Mercury's Orbit", width=win, height=win)

earth = vpython.sphere()
earth.pos = v(0, 0, 0)
earth.radius = 1.5
earth.color = vpython.color.blue

moon = vpython.sphere()
moon.radius = 0.5
moon.color = vpython.color.red
moon.a = 5  # Semi-major axis
moon.e = 0.0  # eccentricity of orbit (0=circle)
moon.pos = v(moon.a * (1 - moon.e), 0, 0)
startpos = v(moon.a * (1 - moon.e), 0, 0)

moon.sidperiod = 27.321661  # Sidereal Period (moon's 'year', in Earth days)
moon.rotperiod = 27.321661  # Length of the planets sidereal day, in Earth days
marrow = vpython.arrow(pos=moon.pos, axis=v(-1, 0, 0))

# tlabel=label(pos=v(7,7,0), text='', xoffset=0, yoffset=0, box=0)

stheta = 0.0
rtheta = 0.0
dt = 0.01
t = 0.0

time.sleep(5)  # Pause for a bit before starting

while 1:  # Do the animation
    vpython.rate(100)
    t = t + dt
    moon.pos = vpython.rotate(startpos, angle=stheta)
    moon.pos.mag = moon.a * (1 - moon.e * moon.e) / (1 + moon.e * math.cos(stheta))
    marrow.pos = moon.pos
    marrow.axis = vpython.rotate(v(-1, 0, 0), angle=rtheta)
    #  tlabel.text="Time (days) %06.2f" % t

    stheta = stheta + (dt / moon.sidperiod) * 2 * math.pi
    rtheta = rtheta + (dt / moon.rotperiod) * 2 * math.pi
    if stheta > 2 * math.pi:
        stheta = stheta - 2 * math.pi
    if rtheta > 2 * math.pi:
        rtheta = rtheta - 2 * math.pi
