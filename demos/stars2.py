import math
import random

import numpy

import vpython
from vpython import vector as v

# Stars interacting gravitationally

win = 1000

Nstars = 300  # change this to have more or fewer stars

G = 6.7e-11  # Universal gravitational constant

# Typical values
Msun = 2E30  # 1 solar mass, in kg
Rsun = 1E10  # approx 7 solar radii
Rtrail = 2e8
L = 4e11
vsun = 0.2 * math.sqrt(G * Msun / Rsun)

scene = vpython.canvas(title="Stars", width=win, height=win,
                       range=2 * L, forward=v(-1, -1, -1))

xaxis = vpython.curve(pos=[v(0, 0, 0), v(L, 0, 0)], color=v(0.5, 0.5, 0.5))
yaxis = vpython.curve(pos=[v(0, 0, 0), v(0, L, 0)], color=v(0.5, 0.5, 0.5))
zaxis = vpython.curve(pos=[v(0, 0, 0), v(0, 0, L)], color=v(0.5, 0.5, 0.5))

Stars = []
colors = [vpython.color.red, vpython.color.green, vpython.color.blue,
          vpython.color.yellow, vpython.color.cyan, vpython.color.magenta]
poslist = []
plist = []
mlist = []
rlist = []

for i in range(Nstars):
    x = L * random.gauss(0, 1)
    y = L * random.gauss(0, 1)
    z = L * random.gauss(0, 1)
    d = math.sqrt(x ** 2 + y ** 2 + z ** 2)
    r = Rsun
    star = vpython.sphere(pos=v(x, y, z), radius=r, color=colors[i % 6])
    star.trail = vpython.curve(pos=[star.pos], color=colors[i % 6], radius=Rtrail)
    Stars.append(star)
    mass = Msun
    px = mass * vsun * (-1 + 2 * random.random()) * (1 - d / (1.7 * L))
    py = mass * vsun * (-1 + 2 * random.random()) * (1 - d / (1.7 * L))
    pz = mass * vsun * (-1 + 2 * random.random()) * (1 - d / (1.7 * L))
    poslist.append([x, y, z])
    plist.append([px, py, pz])
    mlist.append(mass)
    rlist.append(r)

pos = numpy.array(poslist)
p = numpy.array(plist)
m = numpy.array(mlist)
m.shape = (Nstars, 1)
radius = numpy.array(rlist)

vcm = sum(p) / sum(m)  # velocity of center of mass
p = p - m * vcm  # make total initial momentum equal zero

t = 0.0
dt = 10000.0
Nsteps = 0
pos = pos + (p / m) * (dt / 2.)  # initial half-step
Nhits = 0

autoscale = 0
autocenter = 0

while 1:
    # Compute all forces on all stars
    r = pos - pos[:, numpy.newaxis]  # all pairs of star-to-star vectors
    for n in range(Nstars):
        r[n, n] = 1e6  # otherwise the self-forces are infinite
    rmag = numpy.sqrt(numpy.add.reduce(r * r, -1))  # star-to-star scalar distances
    #    hit = less_equal(rmag,radius+radius[:,newaxis])-identity(Nstars)
    #    hitlist = sort(nonzero(hit.flat)) # 1,2 encoded as 1*Nstars+2
    hitlist = []

    F = -G * m * m[:, numpy.newaxis] * r / rmag[:, :, numpy.newaxis] ** 3  # all force pairs
    for n in range(Nstars):
        F[n, n] = 0  # no self-forces
    p = p + sum(F, 1) * dt

    # Having updated all momenta, now update all positions         
    pos = pos + (p / m) * dt

    # Update positions of display objects; add trail
    for i in range(Nstars):
        Stars[i].pos = v(*pos[i])

    # If any collisions took place, merge those stars
    for ij in hitlist:
        i, j = divmod(ij, Nstars)  # decode star pair
        if not Stars[i].visible: continue
        if not Stars[j].visible: continue
        # m[i] is a one-element list, e.g. [6e30]
        # m[i,0] is an ordinary number, e.g. 6e30
        newpos = (pos[i] * m[i, 0] + pos[j] * m[j, 0]) / (m[i, 0] + m[j, 0])
        newmass = m[i, 0] + m[j, 0]
        newp = p[i] + p[j]
        newradius = Rsun * ((newmass / Msun) ** (1. / 3.))
        iset, jset = i, j
        if radius[j] > radius[i]:
            iset, jset = j, i
        Stars[iset].radius = newradius
        m[iset, 0] = newmass
        pos[iset] = newpos
        p[iset] = newp
        Stars[jset].trail.visible = 0
        Stars[jset].visible = 0
        p[jset] = [0, 0, 0]
        m[jset, 0] = Msun * 1E-30  # give it a tiny mass
        Nhits = Nhits + 1
        pos[jset] = [10. * L * Nhits, 0, 0]  # put it far away
    t = t + dt
    vpython.rate(100)
