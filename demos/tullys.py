import vpython
from vpython import vector as v
import random

win = 1000
lab = None

print("""
Visualise the galaxies in the Tully Nearby Galaxy Catalog.

Right button drag to rotate view.
Left button drag up or down to move in or out.

Galaxy orientations are random, and the sizes are scaled up by 10
""")


def mouseclick(m):
    global lab
    pickobj = scene.mouse.pick
    if pickobj:
        if lab:
            lab.visible = 0
        lab = vpython.label(pos=pickobj.pos, xoffset=20, yoffset=20, text=pickobj.name)
    else:
        if lab:
            lab.visible = 0


scene = vpython.canvas(title="2367 Local Galaxies", width=win, height=win, range=30, forward=v(-1, -1, -1))

A = 3.0
xaxis = vpython.curve(pos=[v(0, 0, 0), (A, 0, 0)], color=v(0.5, 0.5, 0.5))
yaxis = vpython.curve(pos=[v(0, 0, 0), (0, A, 0)], color=v(0.5, 0.5, 0.5))
zaxis = vpython.curve(pos=[v(0, 0, 0), (0, 0, A)], color=v(0.5, 0.5, 0.5))

Stars = [vpython.cylinder(pos=v(0, 0, 0), radius=30 / 100.0, color=v(1.0, 0.0, 0.0),
                          length=30 / 1000.0,
                          axis=v(0.0, 1.0, 0.0))]
Stars[0].name = "Milky Way"

scene.bind('click', mouseclick)

for l in open('C:/Users/Andrew//Documents/PyCharm/andrew/OU and OSS code/tully.csv', 'r').readlines():
    name, gtype, radvel, MB, diam, x, y, z, lum = tuple(l.split(','))
    if name == 'Name':
        continue
    gtype, radvel, MB, diam = int(gtype), int(radvel), float(MB), float(diam)
    x, y, z = float(x), float(y), float(z)
    if x == 0 and y == 0 and z == 0:  # Discard LMC and SMC, too close to MW
        continue
    if diam < 2.0:
        diam = 2.0
    if gtype < -2 or gtype > 9:
        Stars.append(vpython.sphere(pos=v(x, y, z), radius=diam / 100.0, color=v(1.0, 1.0, 1.0)))
    else:
        Stars.append(vpython.cylinder(pos=v(x, y, z), radius=diam / 100.0, color=v(1.0, 1.0, 1.0),
                                      length=diam / 1000.0,
                                      axis=v(-1 + 2 * random.random(),
                                             -1 + 2 * random.random(),
                                             -1 + 2 * random.random())))
    Stars[-1].name = name  # Specify the name of the last galaxy added

while 1:
    vpython.rate(10)
