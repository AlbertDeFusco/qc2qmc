#!/usr/bin/env python
from math import pi,sqrt

# The GAMESS primitive normalization
#
# found in inputa.src::ATOMS
#
def primNormGAMESS(prim):
  ang=prim[0]
  alpha=prim[1]
  coeff=prim[2]

  facs =  pi**(3./2.) / ( 2.*alpha * sqrt(2.*alpha))

  # s functions
  if(ang == 0):
    norm = facs

  # p functions
  if(ang == 1):
    norm = 0.5 * facs / (2.*alpha)

  # d functions
  if(ang == 2):
    norm = 0.75 * facs / (4.*alpha**2)

  # f functions
  if(ang == 3):
    norm = 1.875 * facs / (4.*alpha**3)

  return sqrt(norm)

# The GAMESS contraction normalization
#
# In inputa.src::ATOMS, the primitives are normalized first
#
def conNormGAMESS(contraction,normFunction,primNorm=False):
  norm=0.
  for i in contraction:
    ang = i[0]
    ci = i[2] / normFunction(i)
    for j in contraction:
      cj = j[2] / normFunction(j)

      ee=i[1] + j[1]
      fac = ee*sqrt(ee)

      # s primitive
      if(ang == 0):
	prim = ci*cj/fac

      # p primitive
      if(ang == 1):
	prim = 0.5*ci*cj/(ee*fac)

      # d primitive
      if(ang == 2):
	prim = 0.75*ci*cj/(ee*ee*fac)

      # f primitive
      if(ang == 3):
	prim = 1.875*ci*cj/(ee*ee*ee*fac)

      norm = norm + prim

  norm = sqrt(norm*pi**(3./2.))


  return [(i[0],i[1],i[2]/norm/normFunction(i)) if primNorm
      else (i[0],i[1],i[2]/norm)
      for i in contraction]


# CASINO/Crystal contraction normalization
#
# From Mike Towler in CASINO/examples/generic/gauss_dfg
#
def conNormCASINO(contraction): #,normFunction,primNorm=False)
  norm=0.
  for i in contraction:
    ang = i[0]
    ai = i[1]
    di = i[2] #/ normFunction(i)
    for j in contraction:
      aj = j[1]
      dj = j[2] #/ normFunction(j)

      norm = norm + \
          di*dj * ( 2.*sqrt(aj*aj) / (ai + aj) )**(ang + 3./2.)

  norm = sqrt(norm)

  return [(i[0],i[1],i[2]/norm) for i in contraction]

def primNormCASINO(prim):
  ang=prim[0]
  alpha=prim[1]
  coeff=prim[2]

  norm = sqrt(2.**(ang+3./2.)*alpha**(ang+3./2.)) / pi**(3./4.) \
         * sqrt(2.**ang / dFac(2*ang-1))

  return norm


# The double factorial
def dFac(n):
  try:
    return reduce(lambda x,y: y*x, range(n,1,-2))
  except:
    return 1


def checkContraction(contraction):
  ang=contraction[0][0]
  for prim in contraction:
    thisAng=prim[0]
    if (thisAng != ang):
      print "All primitives must have the same angular momentum"
      sys.exit(1)



#Provide the contraction as a list of tuples
# [ (l, alpha, d ), etc. ]
contraction = [ (2,46.1353741080831021977,0.06678829454430918743),
                (2,20.2682182253994397729,0.23122499388298942708),
                (2,6.09459166525985575420,5.07995919900226523237) ]


print "Angular Momentum %3d" % contraction[0][0]
print "  Input Contraction"
for prim in contraction:
  print "    Exp %15.8f; Coeff %15.8f" % \
      (prim[1],prim[2])

print
print
print "  Normalized primitives"
for prim in contraction:
  gmsFac = primNormGAMESS(prim)
  qmcFac  = primNormCASINO(prim)

  print "    Exp %15.8f; GAMESS Coeff %15.8f; CASINO Coeff %15.8f" % \
      (prim[1],prim[2]/gmsFac,prim[2]*qmcFac)

print
print
print "  GAMESS Normalized contractions"
for prim in conNormGAMESS(contraction,primNormGAMESS):
  print "    Exp %15.8f; Coeff %15.8f" % \
      (prim[1],prim[2])
print
print
print "  GAMESS primitive-normalized contractions"
for prim in conNormGAMESS(contraction,primNormGAMESS,True):
  print "    Exp %15.8f; Coeff %15.8f" % \
      (prim[1],prim[2])
print
print
print "  CASINO Normalized contractions"
for prim in conNormCASINO(contraction):
  print "    Exp %15.8f; Coeff %15.8f" % \
      (prim[1],prim[2])
print
print
print "  CASINO primitive-normalized contractions"
for prim in conNormCASINO(contraction):
  print "    Exp %15.8f; Coeff %15.8f" % \
      (prim[1],prim[2]*primNormCASINO(prim))
