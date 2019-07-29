"""Click once in the display window with the mouse, and a line will be drawn
   giant to dwarf star. Half a second later, a second line will be drawn. Use
   this to demonstrate that equal areas are swept out in equal times.
"""

import vpython
from vpython import vector as v

win = 1000
scene = vpython.canvas(title="Orbit", width=win, height=win, range=4e11, forward=v(0, -1, 0), up=v(-1, 0, 0))

giant = vpython.sphere()
giant.pos = v(-1e11, 0, 0)
giant.radius = 2e10
giant.color = vpython.color.red
giant.mass = 1e30
giant.p = v(0, 0, -1e2) * giant.mass

dwarf = vpython.sphere()
dwarf.pos = v(1.5e11, 0, 0)
dwarf.radius = 1e10
dwarf.color = vpython.color.yellow
dwarf.mass = 1e28
dwarf.p = -giant.p

giant.orbit = vpython.curve(color=giant.color, radius=2e9)
dwarf.orbit = vpython.curve(color=dwarf.color, radius=2e9)

dt = 86400

autoscale = 0
autocenter = 0

startbar = None
endbar = None
starttick = -1


def mouseclick(evt):
    global startbar, starttick
    if startbar:
        startbar.visible = 0
    startbar = vpython.cylinder(pos=giant.pos, axis=(dwarf.pos - giant.pos), radius=3e9)
    if endbar:
        endbar.visible = 0
    starttick = 50


scene.bind('click', mouseclick)

while 1:
    vpython.rate(100)

    dist = dwarf.pos - giant.pos
    force = 6.7e-11 * giant.mass * dwarf.mass * dist / vpython.mag(dist) ** 3
    giant.p = giant.p + force * dt
    dwarf.p = dwarf.p - force * dt

    giant.pos += giant.p / giant.mass * dt
    giant.orbit.append(pos=giant.pos)
    dwarf.pos += dwarf.p / dwarf.mass * dt
    dwarf.orbit.append(pos=dwarf.pos)

    starttick = starttick - 1
    if starttick == 0:
        endbar = vpython.cylinder(pos=giant.pos, axis=(dwarf.pos - giant.pos), radius=3e9)
