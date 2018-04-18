
from math import sin, cos, tan, sqrt, pi
import random

AU = 1.49597870700e11   # Astronomical unit in metres
Msun = 19890850e23
G = 6.67259e-11


"""
   Ephemeris library, written by Dr Andrew Williams <andrew@physics.uwa.edu.au>.
   
   Released under the GPL, feel free to use it for non-commercial purposes. I'd appreciate an email with
   what you're using it for, and any nifty improvements you've made.

   Calculate ephemeris data from the JPL planet ephemeris and from standard elliptical orbital elements for
   other objects. 

   It uses the JPL mean ephemeris data from http://ssd.jpl.nasa.gov/txt/p_elem_t2.txt
   Documentation on using this data is from http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf

"""

################ JPL ephemeris data copied into Python structs, with extra rows for the Sun. ############
sun = (   0.00000000,     0.00000000,     0.00000000,       0.00000000,     0.00000000,     0.00000000)
sund = (  0.00000000,     0.00000000,     0.00000000,       0.00000000,     0.00000000,     0.00000000)
mer = (   0.38709843,     0.20563661,     7.00559432,     252.25166724,    77.45771895,    48.33961819)
merd = (  0.00000000,     0.00002123,    -0.00590158,  149472.67486623,     0.15940013,    -0.12214182)
ven = (   0.72332102,     0.00676399,     3.39777545,     181.97970850,   131.76755713,    76.67261496)
vend = ( -0.00000026,    -0.00005107,     0.00043494,   58517.81560260,     0.05679648,    -0.27274174)
ear = (   1.00000018,     0.01673163,    -0.00054346,     100.46691572,   102.93005885,    -5.11260389)
eard = ( -0.00000003,    -0.00003661,    -0.01337178,   35999.37306329,     0.31795260,    -0.24123856)
mar = (   1.52371243,     0.09336511,     1.85181869,      -4.56813164,   -23.91744784,    49.71320984)
mard = (  0.00000097,     0.00009149,    -0.00724757,   19140.29934243,     0.45223625,    -0.26852431)
jup = (   5.20248019,     0.04853590,     1.29861416,      34.33479152,    14.27495244,   100.29282654)
jupd = ( -0.00002864,     0.00018026,    -0.00322699,    3034.90371757,     0.18199196,     0.13024619)
sat = (   9.54149883,     0.05550825,     2.49424102,      50.07571329,    92.86136063,   113.63998702)
satd = ( -0.00003065,    -0.00032044,     0.00451969,    1222.11494724,     0.54179478,    -0.25015002)
ura = (  19.18797948,     0.04685740,     0.77298127,     314.20276625,   172.43404441,    73.96250215)
urad = ( -0.00020455,    -0.00001550,    -0.00180155,     428.49512595,     0.09266985,     0.05739699)
nep = (  30.06952752,     0.00895439,     1.77005520,     304.22289287,    46.68158724,   131.78635853)
nepd = (  0.00006447,     0.00000818,     0.00022400,     218.46515314,     0.01009938,    -0.00606302)
plu = (  39.48686035,     0.24885238,    17.14104260,     238.96535011,   224.09702598,   110.30167986)
plud = (  0.00449751,     0.00006016,     0.00000501,     145.18042903,    -0.00968827,    -0.00809981)

sunx = (  0.00000000,   0.00000000,   0.00000000,   0.00000000)
merx = (  0.00000000,   0.00000000,   0.00000000,   0.00000000)
venx = (  0.00000000,   0.00000000,   0.00000000,   0.00000000)
earx = (  0.00000000,   0.00000000,   0.00000000,   0.00000000)
marx = (  0.00000000,   0.00000000,   0.00000000,   0.00000000)
jupx = ( -0.00012452,   0.06064060,  -0.35635438,  38.35125000)
satx = (  0.00025899,  -0.13434469,   0.87320147,  38.35125000)
urax = (  0.00058331,  -0.97731848,   0.17689245,   7.67025000)
nepx = ( -0.00041348,   0.68346318,  -0.10162547,   7.67025000)
plux = ( -0.01262724,   0.00000000,   0.00000000,   0.00000000)


###################### Extra data (mass and radius) from Perth Observatory astronomy almanac ############
sunp = (  19890850e23,      696000e3   )
merp = (  3.302153e23,     2440.00e3   )
venp = (  48.68959e23,     6051.84e3   )
earp = (  60.47679e23,     6371.01e3   )   # Mass of Earth _and_ Moon
marp = (  6.419079e23,     3389.92e3   )
jupp = (  18991.62e23,       69911e3   )
satp = (  5686.515e23,       58232e3   )
urap = (  868.4831e23,       25362e3   )
nepp = (  1024.655e23,       24624e3   )
plup = (  0.147e23,           1151e3   )


################# Trig functions in degrees to suit values in JPL ephemeris #############################
def sind(a):
  return sin(a*pi/180.0)


def cosd(a):
  return cos(a*pi/180.0)


def tand(a):
  return tan(a*pi/180.0)


################## Collect JPL ephemeris and other data into tuples in a dictionary ####################
planetelements = {'Sun':(sun,sund,sunx,sunp),
                  'Mercury':(mer,merd,merx,merp),
                  'Venus':(ven,vend,venx,venp),
                  'Earth':(ear,eard,earx,earp),
                  'Mars':(mar,mard,marx,marp),
                  'Jupiter':(jup,jupd,jupx,jupp),
                  'Saturn':(sat,satd,satx,satp),
                  'Uranus':(ura,urad,urax,urap),
                  'Neptune':(nep,nepd,nepx,nepp),
                  'Pluto':(plu,plud,plux,plup)    }


########################### Class and Function definitions #############################################
def vfit(plist):
  """Carry out a linear regression fit on the (t,p) tuples in X, Y, and Z, to find
     a mean velocity vector. This is done to estimate the initial velocity of a body
     by calculating positions in 3D over a few days of orbit before and after the nominal
     time, and fitting the calculated positions for the velocity vector.
  """
  N = len(plist)
  xb, Xyb,Yyb,Zyb = 0.0, 0.0, 0.0, 0.0
  Xsxy,Ysxy,Zsxy = 0.0, 0.0, 0.0
  sxx = 0.0
  
  for x,p in plist:
    xb += x
    Xyb += p[0]
    Yyb += p[1]
    Zyb += p[2]
    sxx += x*x
    Xsxy += x*p[0] 
    Ysxy += x*p[1] 
    Zsxy += x*p[2] 
  xb = xb/N
  sxx -= N*xb*xb
  Xsxy -= N*xb*Xyb
  Ysxy -= N*xb*Yyb
  Zsxy -= N*xb*Zyb
  vel = (Xsxy/sxx, Ysxy/sxx, Zsxy/sxx)
#  print name, plist, vel
  return vel


class Planet(object):
  """Define a class that represents a body in orbit, with methods to calculate 
     position and velocity in 3D Cartesian coordinates.
  """
  def __init__(self, name, etype='elliptical', elements=None, 
               a=1.0, e=0.0, I=0.0, M0=2451545.0, epoch=2451545.0, w=0.0, omg=0.0,
               mass=1e23, rad=1e3):
    """If etype is 'planet', save the JPL ephemeris data in 'elements', otherwise use
       'a', 'e', etc directly, where the arguments are:
          a: Semi major axis in Astronomical Units
          e: Eccentricity
          I: Inclination in degrees
          M0: Mean anomaly at epoch
          Epoch: Epoch of elements
          w: Argument of perihelion
          omg: Longitude of ascending node
    """
    self.pos = (0.0, 0.0, 0.0)
    self.vel = (0.0, 0.0, 0.0)
    if etype == 'planet':
      elm, eld, elx, elp = elements
      self.a0, self.e0, self.I0, self.L0, self.longperi0, self.longnode0 = elm
      self.ad, self.ed, self.Id, self.Ld, self.longperid, self.longnoded = eld
      self.b, self.c, self.s, self.f = elx
      self.M, self.rad = elp
      self.degperday = self.Ld/36525
    elif etype == 'elliptical':
      self.a = a
      self.e = e
      self.I = I
      self.M0 = M0
      self.M = mass
      self.rad = rad
      self.epoch = epoch
      self.w = w
      self.omg = omg
      P = sqrt(4*pi*pi*a*a*a*AU*AU*AU / (G*Msun))
      self.Md = 86400*360.0/P     # Rate of change of mean anomaly, in degrees per day.
      self.degperday = self.Md
    self.name = name
    self.body = None
    self.etype = etype

  def calcpos(self, Teph=2451545.0):
    """Calculate the 3D Cartesian coordinates at time Teph
    """
    if self.etype == 'planet':
      T = (Teph-2451545.0)/36525
      a = self.a0 + T*self.ad
      e = self.e0 + T*self.ed
      I = self.I0 + T*self.Id
      L = self.L0 + T*self.Ld
      b, c, s, f = self.b, self.c, self.s, self.f
      wb = self.longperi0 + T*self.longperid
      omg = self.longnode0 + T*self.longnoded

      w = wb - omg
      M = L - wb + b*T*T + c*cosd(f*T) + s*sind(f*T)
    elif self.etype == 'elliptical':
      a = self.a
      e = self.e
      I = self.I
      w = self.w
      omg = self.omg
      M = self.M0 + (Teph-self.epoch)*self.Md
    else:
      return

    while M > 180.0:
      M -= 360.0
    while M < -180.0:
      M += 360.0             # Take to module -180 -> +180
    estar = e*180.0 / pi
    dE = 1000
    E = M + estar*sind(M)
    while dE > 1e-6:
      dM = M - (E - estar*sind(E))
      dE = dM / (1 - e*cosd(E))
      E = E + dE

    xp = a*(cosd(E)-e)
    yp = a*sqrt(1-e*e)*sind(E)
    zp = 0.0

    xecl = (cosd(w)*cosd(omg)-sind(w)*sind(omg)*cosd(I))*xp + (-sind(w)*cosd(omg)-cosd(w)*sind(omg)*cosd(I))*yp
    yecl = (cosd(w)*sind(omg)+sind(w)*cosd(omg)*cosd(I))*xp + (-sind(w)*sind(omg)+cosd(w)*cosd(omg)*cosd(I))*yp
    zecl = (sind(w)*sind(I))*xp + (cosd(w)*sind(I))*yp

    self.pos = (xecl, yecl, zecl)

  def calc(self, Teph=2451545.0):
    """Calculate the 3D Cartesian coordinates and velocity vector at time Teph
    """
    if self.name != 'Sun':
      plist = []
      for i in [-0.5,0,0.5]:
        self.calcpos(Teph=Teph+i/self.degperday)
        plist.append((i/self.degperday,self.pos))
      self.vel = vfit(plist)
    else:
      self.vel = (0.0, 0.0, 0.0)
    self.calcpos(Teph=Teph)


####################### Helper functions to collect and return planets and asteroids #########################
      
def GetPlanets(oclass=Planet, Teph=2451545.0):
  """Return a dict of objects of type 'oclass' for the Sun, Mercury, ... Pluto
     using the JPL ephemeris, plus an object for Ceres using the real orbital elements.
  """
  planets = {}
  for name,e in planetelements.items():
    p = oclass(name, etype='planet', elements=e)
    p.calc(Teph=Teph)
    planets[name] = p

  ceres = oclass('Ceres',etype='elliptical',
                 epoch=53600+2400000.5, a=2.7659794, e=0.08001400,
                 I=10.58604, w=73.39245, omg=80.40973, M0=86.9543947,
                 mass=0.0095e23, rad=948e3)
  ceres.calc(Teph=Teph)
  planets['Ceres'] = ceres
  return planets


def GetRandomAsteroids(oclass=Planet, Teph=2451545.0,
                       amin=1.5, amax=5.2,
                       emin=-0.1, emax=0.1,
                       Imin=-0.5, Imax=0.5,
                       rpos=None,
                       N=100):
  """Return a dict of objects of type 'oclass' for N random asteroids. You can fine tune the
     distributions of orbital elements and masses however you want...
  """
  asteroids = {}
  for mp in range(N):
    a = random.random()*(amax-amin) + amin
    e = random.random()*(emax-emin) + emin
    I = random.random()*(Imax-Imin) + Imin
    if rpos is None:
      omg = random.random()*360.0 - 180.0
      M0 = random.random()*360.0 - 180.0
      w = random.random()*360.0 - 180.0
    else:
      omg = rpos
      M0 = rpos
      w = rpos
    mass = 1e9
    rad = 50e3

    p = oclass(name='MP%d' % mp, etype='elliptical', epoch=2451545.0, a=a, e=e, I=I, w=w, omg=omg, M0=M0, mass=mass, rad=rad)
    p.calc(Teph=Teph)
    asteroids['MP%d' % mp] = p
  return asteroids

  
def GetEvenAsteroids(oclass=Planet, Teph=2451545.0,
                     amin=1.5, amax=5.2,
                     emin=-0.1, emax=0.1,
                     Imin=-0.5, Imax=0.5,
                     rpos=None,
                     N=100):
  """Return a dict of objects of type 'oclass' for N random asteroids. You can fine tune the
     distributions of orbital elements and masses however you want...
  """
  asteroids = {}
  for mp in range(N):
    a = amin + (float(mp)/N)*(amax-amin)
    e = random.random()*(emax-emin) + emin
    I = random.random()*(Imax-Imin) + Imin
    if rpos is None:
      omg = random.random()*360.0 - 180.0
      M0 = random.random()*360.0 - 180.0
      w = random.random()*360.0 - 180.0
    else:
      omg = rpos
      M0 = rpos
      w = rpos
    mass = 1e9
    rad = 50e3

    p = oclass(name='MP%d' % mp, etype='elliptical', epoch=2451545.0, a=a, e=e, I=I, w=w, omg=omg, M0=M0, mass=mass, rad=rad)
    p.calc(Teph=Teph)
    asteroids['MP%d' % mp] = p
  return asteroids
