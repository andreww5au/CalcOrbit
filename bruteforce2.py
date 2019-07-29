"""N-body modelling engine, written by Dr Andrew Williams <andrew@longtable.org>.
   
   Not intended for serious use, as the title suggests it's a dumb brute-force algorithm in pure Python,
   intended for teaching purposes.

   Released under the GPL, feel free to use it for non-commercial purposes. I'd appreciate an email with
   what you're using it for, and any nifty improvements you've made.
"""

import numpy

# try:
#  import psyco
#  psyco.full()
#  from psyco.classes import *
# except:
#  print "Psycho not installed."


# Simulation radius in AU
EscapeDistance = 20

# Physical constants:
G = 6.67259e-11
AU = 1.49597870700e11

# data type used for floating point arrays
dtp = numpy.dtype('float64')

# Lists used to store bodies added to the simulation. The data in these lists is used by the 'Init'
# function to initialise the numpy arrays used for the actual processing, and not updated or used
# after that. For position and momentum vectors during the simulation (after the first call to 'Process'),
# use the numpy arrays (hpos, hp, lpos, lp).
hnamelist = []
lnamelist = []

hposlist = []
hplist = []
hmlist = []
hrlist = []
hbodies = {}

lposlist = []
lplist = []
lbodies = {}


################ Functions to set up the initial simulation ###############################
def AddHeavyBody(name='', body=None, pos=None, vel=None, M=None, r=None):
    """Add an object to the full Nbody simulation.

       Arguments are:
         name - a string containing the name of the body (optional)
         body - an arbitrary object with a 'pos' attribute. This 'pos' attribute
                is updated regularly by the Nbody code. Typically used to pass
                a Visual Python object that should be updated with the position
                from the running simulation
         pos -  The initial position of the object (NOT taken from the 'body.pos' value)
         vel -  The velocity vector (in AU/day)
         M -    The mass of the body, in kilograms
         r -    The capture radius of the body, in metres, used for collision detection
    """
    hnamelist.append(name)
    hposlist.append(numpy.array(pos))
    hplist.append(numpy.array(vel) * M)
    hmlist.append(M)
    hrlist.append(r)
    i = len(hnamelist) - 1
    hbodies[i] = body
    return i


def AddLightBody(name='', body=None, pos=None, vel=None, M=None, r=None):
    """Add an object to the Nbody simulation, but modelled as a massless particle.

       Arguments are:
         name - a string containing the name of the body (optional)
         body - an arbitrary object with a 'pos' attribute. This 'pos' attribute
                is updated regularly by the Nbody code. Typically used to pass
                a Visual Python object that should be updated with the position
                from the running simulation
         pos -  The initial position of the object (NOT taken from the 'body.pos' value)
         vel -  The velocity vector (in AU/day)
         M -    ignored
         r -    ignored
    """
    lnamelist.append(name)
    lposlist.append(numpy.array(pos))
    lplist.append(numpy.array(vel))
    i = len(lnamelist) - 1
    lbodies[i] = body
    return i


####################### Functions to handle collisions and escapes ###############################
def CollideLight(t, ip, imp):
    """Handle collision of a 'Light' (massless) body with a 'Heavy' body.
       Arguments are:
         t -   The time of the collision event, in JD
         imp - The index of the massless body
         ip -  The index of the heavy body
    """
    if imp in lbodies.keys():
        lremaining[imp] = False
        print("%s: Minor planet %s colided with %s." % ((t - 2451545.0) / 365.246 + 2000.0, lnamelist[imp], hnamelist[ip]))
        lbodies[imp].visible = False
        try:
            lbodies[imp].trail_object.pos = []
        except:
            pass
        del lbodies[imp]


def EscapeLight(t, imp):
    """Handle the case where a 'Light' (massless) body exceeds the boundary threshold
       for the Nbody simulation.

       Arguments are:
         t -   The time of the collision event, in JD
         imp - The index of the massless body
    """
    if imp in lbodies.keys():
        lremaining[imp] = False
        print("%s: Minor planet %s left the Solar System." % ((t - 2451545.0) / 365.246 + 2000.0, lnamelist[imp]))
        lbodies[imp].visible = False
        try:
            lbodies[imp].trail_object.pos = []
        except:
            pass
        del lbodies[imp]


######################### Called externally during the simulation ################################
def UpdatePositions(vclass=None):
    """For each of the 'body' objects provided by the AddHeavyBody or AddLightBody calls,
       update the 'pos' attribute in that object with the current position of the body in
       the simulation.

       The vclass parameter should be a vpython vector.

       This function is not called by the 'Process' code, it must be called externally, from
       the same main loop that calls the 'Process' function.
    """
    for i, body in hbodies.items():
        body.pos = vclass(*tuple(hpos[i]))
    for i, body in lbodies.items():
        body.pos = vclass(*tuple(lpos[i]))


########################## Initialise the simulation before processing ########################
def Init(sun=0):
    """Initialise the numpy arrays containing the positions, momenta, masses and radii of both light and
       heavy objects added to the simulation. Also precalculate some arrays that stay depend on parameters
       that stay constant during the simulation (object masses and radii). The only argument, 'sun', is the
       index of the largest mass (in the 'heavy body' list) in the simulation, nominally the Sun. The
       initial momentum of this mass is discarded, and a new momentum is assigned to make the total
       momentum of the simulation sum to zero, for no net motion of the system.

       Arguments are:
         sun - index of the SUn in the Heavy Body list.
    """
    global Nheavy, hpos, hp, hm, hradius2, Nlight, lpos, lp, lremaining, hid, hdumn, hmassprod, lmassprod

    # Heavy bodies
    Nheavy = len(hposlist)  # Number of 'heavy' bodies
    hpos = numpy.array(hposlist, dtp)  # Nheavy*3 array of heavy body positions, in AU
    hp = numpy.array(hplist, dtp)  # Nheavy*3 array of heavy body momenta, in kg*AU/day
    hm = numpy.array(hmlist, dtp)  # Nheavy*1 array of  masses, in kg
    hm.shape = (Nheavy, 1)  # Reshape to (1 by Nheavy) vs. (Nheavy by 1)
    hradius2 = numpy.array(hrlist) / AU
    hradius2 *= hradius2  # Nheavy*1 array of capture radii squared, in (AU^2)
    hradius2 = hradius2[:, numpy.newaxis]  # Add an extra axis, so we don't need to do this during processing

    # Light bodies
    Nlight = len(lposlist)  # Number of 'light' bodies
    lpos = numpy.array(lposlist, dtp)  # Nlight*3 array of light body positions, in AU
    lp = numpy.array(lplist, dtp)  # Nlight*3 array of light body momenta, in kg*AU/day
    lremaining = numpy.ones(Nlight, bool)  # Boolean array with 'True' indicating a light body that's
    #  still in the simulation

    # Constant arrays, created here to save time during processing
    hid = numpy.arange(Nheavy)
    hdumn = numpy.ones(Nheavy, dtp)

    # Make the total system momentum sum to zero
    hp[sun] = numpy.zeros(3, dtp)
    memtot = sum(hp, 0)
    hp[sun] = -memtot

    # Calculate G*M1*M2 for every heavy-heavy and heavy-light pair of bodies, using the
    # units of AU for distance, and days for time, instead of metres and seconds
    hmassprod = G * hm * hm[:, numpy.newaxis] * 86400 * 86400 / (AU * AU * AU)
    lmassprod = G * hm[:, numpy.newaxis] * 86400 * 86400 / (AU * AU * AU)


###################################### Main simulation ######################################
def Process(t, steps=1, dt=1.0, numcores=None):  # Euler-Cromer method
    """Calculate 'steps' iterations in the Nbody simulation, starting at time 't' with a time
       delta of 'dt' (in days) between each step.
       If any light bodies approach to within the capture radius of any heavy bodies, call
       'CollideLight' to handle the collision. If any light bodies exceed 'EscapeDistance'
       (in AU) from (0,0,0), call 'EscapeLight' to handle the escape.

       Arguments are:
         t -      Initial time in days (arbitrary zero point)
         steps -  Number of iterations to process
         dt -     Time step for each iteration, in days
    """
    global hpos, hp, lpos, lp
    for j in range(steps):
        # First handle the motion of the 'Heavy' bodies
        r = hpos[:, numpy.newaxis] - hpos  # all pairs of heavy-heavy position vectors
        rmag3 = numpy.add.reduce(r * r, -1)  # square of scalar distances between heavy objects
        rmag3 *= numpy.sqrt(rmag3)  # Scalar distances cubed
        rmag3[hid, hid] = hdumn  # Set all self-self distances to 1, not 0, to avoid zero division
        #  When multiplied by position vectors, self-self vectors will be
        #  zero, so self-self forces will also be zero
        F = hmassprod * r / rmag3[:, :, numpy.newaxis]  # all heavy-heavy force vector pairs

        hp += sum(F, 1) * dt  # For each heavy object, sum the forces from all the other heavy objects
        #  and add it to that objects current momentum
        hpos += hp * dt / hm  # Position = position + momentum*dt/mass

        if Nlight:
            # Now handle the motion of the 'Light' bodies
            r = hpos[:, numpy.newaxis] - lpos  # all pairs of heavy-light position vectors
            rmag3 = numpy.add.reduce(r * r, -1)  # square of scalar distances

            # While we have the squares of the heavy-light scaler distances, detect any collisions
            colliding = numpy.logical_and(rmag3 < hradius2, lremaining).nonzero()
            if len(colliding[0]):
                for ip, imp in zip(colliding[0], colliding[1]):
                    CollideLight(t, ip, imp)

            rmag3 *= numpy.sqrt(rmag3)  # rmag3 is now the heavy-light scalar distances cubed
            F = lmassprod * r / rmag3[:, :, numpy.newaxis]  # all heavy-light force vector pairs

            lp += sum(F, 0) * dt  # For each light object, sum the forces from all the heavy objects
            #  and add it to that objects current momentum
            lpos += lp * dt  # position = position + momentum*dt (assume mass is 1.0)

    if Nlight:
        # Use the new position vectors to detect any light bodies that have left the modelling volume
        escaped = numpy.logical_and(numpy.logical_or.reduce(abs(lpos) > EscapeDistance, -1), lremaining).nonzero()
        if len(escaped[0]):
            for imp in escaped[0]:
                EscapeLight(t, imp)
