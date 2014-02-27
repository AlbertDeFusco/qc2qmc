#!/usr/bin/env python
import sys
from math import pi,sqrt

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

def conNorm(contraction,normFunction):
  out=list()
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

  for i in contraction:
    out.append((i[0],
                i[1],
                i[2] / norm))

  return out



def dFac(n):
  try:
    return reduce(lambda x,y: y*x, range(n,1,-2))
  except:
    return 1

def primNormBe(prim):
  ang=prim[0]
  alpha=prim[1]
  coeff=prim[2]

  piAlpha = (2.*alpha/pi)**(3./2.)

  return sqrt(piAlpha * (4.*alpha)**(ang) / dFac(2.*ang-1))
  

def checkContraction(contraction):
  ang=contraction[0][0]
  for prim in contraction:
    thisAng=prim[0]
    if (thisAng != ang):
      print "All primitives must have the same angular momentum"
      sys.exit(1)



#Provide the contraction as a list of tuples
# [ (l, alpha, c ), etc. ]
contraction = list()

try:
  contraction = sys.argv[1]
except:
  contraction = [ (2,46.1353741080831021977,0.06678829454430918743),
                  (2,20.2682182253994397729,0.23122499388298942708),
		  (2,6.09459166525985575420,5.07995919900226523237) ]

checkContraction(contraction)


print "Angular Momentum %3d" % contraction[0][0]
print "  Input Contraction"
for prim in contraction:
  print "    Exponent %15.8f; %15.8f" % \
      (prim[1],prim[2])

print
print
print "  Normalized primitives"
for prim in contraction:
  gmsFac = primNormGAMESS(prim)
  BeFac  = primNormBe(prim)

  print "    Exponent %15.8f; GMS %15.8f; Be %15.8f" % \
      (prim[1],prim[2]/gmsFac,prim[2]/BeFac)

print
print
print "  GAMESS Normalized contractions"
for prim in conNorm(contraction,primNormGAMESS):
  print "    Exponent %15.8f; %15.8f" % \
      (prim[1],prim[2])
print
print
print "  GAMESS Normalized contractions"
for prim in conNorm(contraction,primNormBe):
  print "    Exponent %15.8f; %15.8f" % \
      (prim[1],prim[2])
