
import math
from random import random

import numpy

import vpython
from vpython import vector as v

# Stars interacting gravitationally

win = 600

Nstars = 20  # change this to have more or fewer stars

G = 6.7e-11  # Universal gravitational constant

# Typical values
Msun = 2E30
Rsun = 2E9
Rtrail = 2e8
L = 4e10
vsun = 0.8 * math.sqrt(G * Msun / Rsun)

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
    x = -L + 2 * L * random()
    y = -L + 2 * L * random()
    z = -L + 2 * L * random()
    r = Rsun / 2 + Rsun * random()
    star = vpython.sphere(pos=v(x, y, z), radius=r, color=colors[i % 6])
    star.trail = vpython.curve(pos=[star.pos], color=colors[i % 6], radius=Rtrail)
    Stars.append(star)
    mass = Msun * r ** 3 / Rsun ** 3
    px = mass * (-vsun + 2 * vsun * random())
    py = mass * (-vsun + 2 * vsun * random())
    pz = mass * (-vsun + 2 * vsun * random())
    poslist.append([x, y, z])
    plist.append([px, py, pz])
    mlist.append(mass)
    rlist.append(r)

pos = numpy.array(poslist)
p = numpy.array(plist)
m = numpy.array(mlist)
m.shape = (Nstars, 1)
radius = numpy.array(rlist)

vcm = numpy.sum(p) / numpy.sum(m)  # velocity of center of mass
p = p - m * vcm  # make total initial momentum equal zero

dt = 1000.0
Nsteps = 0
pos = pos + (p / m) * (dt / 2.)  # initial half-step
Nhits = 0

while True:
    vpython.rate(100)

    # Compute all forces on all stars
    r = pos - pos[:, numpy.newaxis]  # all pairs of star-to-star vectors
    for n in range(Nstars):
        r[n, n] = 1e6  # otherwise the self-forces are infinite
    rmag = numpy.sqrt(numpy.add.reduce(r * r, -1))  # star-to-star scalar distances
    hit = numpy.less_equal(rmag, radius + radius[:, numpy.newaxis]) - numpy.identity(Nstars)
    hitlist = numpy.sort(numpy.nonzero(hit.flat)[0]).tolist()  # 1,2 encoded as 1*Nstars+2
    F = G * m * m[:, numpy.newaxis] * r / rmag[:, :, numpy.newaxis] ** 3  # all force pairs

    for n in range(Nstars):
        F[n, n] = 0  # no self-forces
    p = p + numpy.sum(F, 1) * dt

    # Having updated all momenta, now update all positions         
    pos = pos + (p / m) * dt

    # Update positions of display objects; add trail
    for i in range(Nstars):
        Stars[i].pos = v(*pos[i])
        if Nsteps % 10 == 0:
            Stars[i].trail.append(pos=Stars[i].pos)

    # If any collisions took place, merge those stars
    for ij in hitlist:
        i, j = divmod(ij, Nstars)  # decode star pair
        if (not Stars[i].visible) or (not Stars[j].visible):
            continue
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

    Nsteps += 1
