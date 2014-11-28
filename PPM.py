#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, sys, time
from math import *

########### import main module named Maxwell from ./build/ directory ################
sys.path.append("build")
from Maxwell import *

########### Set the size of computational grid ###########
########### It's size will be nx*ny*nz cells ##########
nbx,nby,nbz = 24,1,48 #8, 1, 24
nx,ny,nz = 2**MaxRank*nbx, 2**MaxRank*nby, 2**MaxRank*nbz   # MaxRank is predifined constant parameter in file src/params.hpp
Bricks = getBricks()
Bricks.set(nbx,nby,nbz);

########### Components to drop 3D field components at constant time moments #######################################
########### (available electric field components:  Ex,Ey,Ez,
###########            magnetic field components:  Hx,Hy,Hz,
###########            Pointing vector components: Sx,Sy,Sz,
###########                             Material:  Mat       )   #######
DropP = getDropP()
#DropP.Fields4Drop = Components(["Ex","Ey","Ez","Ef","Mat"])
DropP.Fields4Drop = Components(["Mat"])

########### Components to save in sensors ##############
#DropP.Fields4Sensors = Components(["Ex","Ey","Ez", "Hx","Hy","Hz"])
DropP.Fields4Sensors = Components(["Ex","Ey","Ez","Hx","Hy","Hz","Ef"])

########### Directory for log files (LogPath) and for drop files (DropPath, sensors files are in DropPath too) ######################
calc_path="./PPM/la"+(sys.argv[1].replace(",","."))+"_a"+(sys.argv[2].replace(",","."))
if not os.path.exists(calc_path): os.mkdir(calc_path)
DropP.LogPath  = calc_path
DropP.DropPath = calc_path

########### structures for control main parameters, Signal parameters, and material parameters ########
########### structures for Signal and material parameters could be created or modified in files Signal.hpp and materials.hpp #######
pars = getPars()
parsSignal, parsMat = getParsSrc(), getParsMat()

########### set spatial and temporal steps (don't forget about Courant condition, light speed = 1) ##########
pars.setGrid(dx = 2*5./512, dy = 4*5./512, dz = 0.25*5./512, dt = 0.7/512)
pars.deep1d = 8
pars.deep3d = 5

########### number of programm threads will be created ##########
pars.Nthreads = 16
pars.timeSep = 3

XLength = pars.dx*pars.Nx*(2**MaxRank)
YLength = pars.dy*pars.Ny*(2**MaxRank)
ZLength = pars.dz*pars.Nz*(2**MaxRank)
########### set some signal parameters (look in file src/Signal.hpp) ##############
parsSignal.setPars()
parsMat.setPars()

la = float(sys.argv[1].replace(",","."));			# 750 --- 850 nm
alpha = float(sys.argv[2].replace(",","."));
nAir = 1.00028
parsSignal.Lambda = la/nAir;
parsSignal.Xcnt = 0.5*XLength;
parsSignal.Ycnt = 0.5*YLength;
parsSignal.Th = 30.0;					# 100.3 fs
parsSignal.F0 = 1.0/parsSignal.Th; 			
parsSignal.Omega = 2.0*pi*parsSignal.F0;		# 2*PI*mu  
parsSignal.omega = 2.0*pi/parsSignal.Lambda		# 2*PI*c/la  
parsSignal.alpha = alpha*pi/180.0; 		
parsSignal.k0 = nAir*parsSignal.omega/1.;		# omega/c
parsSignal.kx = parsSignal.k0*sin(parsSignal.alpha);		
parsSignal.ky = 0.0;		
parsSignal.kz = parsSignal.k0*cos(parsSignal.alpha);		
parsSignal.Px = parsSignal.A0*cos(parsSignal.alpha); 
parsSignal.Py = 0.0; 			
parsSignal.Pz = parsSignal.A0*sin(parsSignal.alpha);

########### set material parameters (look in file src/materials.hpp) ###############
parsMat.X0 = 0.5*parsMat.Xlength;
parsMat.Y0 = 0.5*parsMat.Ylength;

parsMat.zGGG = 0.0;
parsMat.thinGGG =  0.2;

parsMat.thinBIG = 2.2;
parsMat.zBIG = parsMat.zGGG + 0.5*parsMat.thinGGG + 0.5*parsMat.thinBIG;

parsMat.thinGold = 0.04;
parsMat.zGold = parsMat.zGGG + 0.5*parsMat.thinGGG + parsMat.thinBIG + 0.5*parsMat.thinGold;
parsMat.period = 0.73;
parsMat.cell = 0.63;
parsMat.shift = 0.5*parsMat.cell;

########### Set sensors coordinates (in computational grid) ############
print "Xstep\tYstep"
print int(1./16.*XLength/pars.dx), int(XLength/pars.dx),"\t", int(1./16.*YLength/pars.dy), int(YLength/pars.dy)
zSensor = 20;
#for ix in range(0, int(XLength/pars.dx), int(1./16.*XLength/pars.dx)):
#	for iy in range( 0, int(YLength/pars.dy), int(1./16.*YLength/pars.dy) ):
#		DropP.setSensor(ix, iy, int(zSensor) );
DropP.setSensor(int(0.5*XLength/pars.dx), int(0.5*YLength/pars.dy), int(zSensor) );
InitializeCoffs();

print "Xl, Yl, Zl"
print pars.dx*pars.Nx*(2**MaxRank)
print pars.dy*pars.Ny*(2**MaxRank)
print pars.dz*pars.Nz*(2**MaxRank)
########### Set material dielectic and magnetic coefficients ##########
epsVac,muVac = Array(3),Array(3)
for i in 0,1,2: epsVac[i] = 1.; muVac[i] = 1.;
#setCoffsIso(IndVacum, eps=1., mu=1.);

epsDiel = Array(3);
for i in 0,1,2: epsDiel[i] = 6.+i;
setCoffsIso   (IndAir, eps = (nAir)**2,		mu = 1.);
setCoffsIso   (IndBIG, eps = 6.066-0.0011*la,	mu = 1.);
setCoffsIso   (IndGGG, eps = (1.97)**2, 	mu = 1.);

#setCoffsIso   (IndDiel, eps = 1.2*1.2, mu = 1.);
#setCoffsAniso (IndDiel, eps = epsDiel, mu = muVac);
#setCoffsDisp  (IndDiel, eps_inf = 1,   mu = 1, sigma=0, Deps1=8.93, omega1=3.42*2*pi, gamma1=0.425*2*pi, gamma1_=0.087*2*pi, Deps2=1.855    , omega2=2.72*2*pi, gamma2=0.123*2*pi, gamma2_=2.678*2*pi)

ld1 = LorenzDisp()
ld2 = LorenzDisp()
ld1.setDrude(Deps=8.5, wp=4.9*pi, gp=0.21, gp_=1e-3)
ld2.set(Deps=1.0, wp=0.5*2*pi, gp=30*2*pi, gp_=340*pi)
setCoffsDisp(IndGold, eps_inf = 1, mu = 1, LorenzD = DispParams([ld1,ld2]), sigma=0)

print "Center", parsSignal.Xcnt, parsSignal.Ycnt
print pars.dx, pars.dy
########### Initialize Data with zero ##############
InitializeData()

########### Begin count processing ##########
########### Here will be counted about 300*(2^MaxRank) time steps #########
def run():
  print "Starting..."
  for iT in range( int(30*parsSignal.Th) ):
    if iT> 0: UpdateOneStep(iT)
    if iT> 30: DropP.Fields4Drop = Components([]) #Components(["Ex","Ey","Ez","Hx","Hy","Hz"])
run()
