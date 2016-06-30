#=============================================================================
# Author: Benjamin Buckman
# E=mail: buckman.12@osu.edu
# Description: Integrate galprop 3-D emissivity maps for arbitrary inclination through line of sight.
# Output: fits file of image
#=============================================================================
from scipy import *
import sys
sys.path.append('./pyfits/lib')

#=============================================================================
# synch_line_of_sight(1,2,3,4,5)
# Parameters:
#		1) emissivity filename 
#		2) fits output filename 
#		3) inclination in degrees [phi,theta,psi]
#		4) physical distance to object (kpc)
#		5) instrument half power beam width (hpbw) in arcsec
#==============================================================================
def synch_line_of_sight(emissfilename, fileout, inclination, objectdist, hpbw):
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	inc_rad = inclination*math.pi/180. #inclination in radians
	
	#READING FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	ymin = head['CRVAL2']
	zmin = head['CRVAL3']
	log_n_min = head['CRVAL4'] #log10(nu_min)
	
	xdel = head['CDELT1']
	ydel = head['CDELT2']
	zdel = head['CDELT3']
	nfactor = head['CDELT4'] #factor
	
	xnum = head['NAXIS1']
	ynum = head['NAXIS2']
	znum = head['NAXIS3']
	nnum = head['NAXIS4']
	
	nmin = 10.**log_n_min  #Hz
	
	xmax = xmin+xdel*xnum
	ymax = ymin+ydel*ynum
	zmax = zmin+zdel*znum
	
	xsize = xmax-xmin
	ysize = ymax-ymin
	zsize = zmax-zmin

	orPos = np.zeros(3) #original position (x,y,z) kpc	
	orxbins = np.zeros(xnum)
	orybins = np.zeros(ynum)
	orzbins = np.zeros(znum)
	for n in range(0,xnum):
		orxbins[n] = xmin+xdel*n
	for n in range(0,ynum):
		orybins[n] = ymin+ydel*n
	for n in range(0,znum):
		orzbins[n] = zmin+zdel*n
	
	original_mapsize = math.sqrt(xsize**2.+ysize**2.+zsize**2.) #kpc
	rotate_mapsize = 1.01*original_mapsize #kpc

	#OUTPUT IMAGE------------------------------------------------------
	outRes = 41 #number of pixels = outRes**2
	int_steps = 100 #number of integration steps
	
	kpcPix = rotate_mapsize/outRes #kpc/pixel_dim
	image = np.zeros([nnum,outRes,outRes]) #output map (pixels)
	imPos = np.zeros(3) #image position [x,y,z] kpc
	
	xminIm = -rotate_mapsize/2. #kpc
	yminIm = -rotate_mapsize/2. #kpc
	zminIm = -rotate_mapsize/2. #kpc
	
	dxIm = rotate_mapsize/(outRes-1.) #kpc
	dyIm = rotate_mapsize/(outRes-1.) #kpc
	dzIm = rotate_mapsize/(int_steps-1.) #kpc
	
	cmxPix = dxIm*kpc2cm #cm/pixel
	cmyPix = dyIm*kpc2cm #cm/pixel
	cmzPix = dzIm*kpc2cm #cm/pixel
	
	cmDist = objectdist*kpc2cm
	solidangle = 1./((objectdist*kpc2cm)**2.) #4pi*ratio of 1cm^2 to entire sphere
	volPix = cmxPix*cmyPix*cmzPix #cm^3
	
	srPix = dxIm*dyIm/objectdist**2. #steradian/pixel
	Pix2sac = srPix*Sr2sqArcSec #arcsec^2/pixel
	sac2Pix = 1./Pix2sac #pixel/arsec^2
	
	beamConversion = Pix2sac*math.pi*(hpbw/2.)**2./math.log(2.)
	
	print('ArcSec^2/pixel: '+str(Pix2sac))
	print('Beam Conversion: '+str(beamConversion))
	
	#INTEGRATION---------------------------------------------------------
	print('=====================================================')
	print('Beginning Integration')
	print('-----------------------------------------------------')
	
	for ip in range(0,nnum):  # p
		frequency=nmin*nfactor**ip
		
		for ix in range(0,outRes):  # x
			for iy in range(0,outRes):  # y
				flux = 0.											# total flux
				for iz in range(0,int_steps):  # z
					imPos[0] = xminIm+ix*dxIm
					imPos[1] = yminIm+iy*dyIm
					imPos[2] = zminIm+iz*dzIm
				
					orPos = rotate_coordinates(imPos,inc_rad[0],inc_rad[1],inc_rad[2])
				
					if (orPos[0] >= xmin-xdel/2. and orPos[0] <= xmax+xdel/2. and
							orPos[1] >= ymin-ydel/2. and orPos[1] <= ymax+ydel/2. and
							orPos[2] >= zmin-zdel/2. and orPos[2] <= zmax+zdel/2.):
						orxindex = min(range(len(orxbins)), key=lambda i: abs(orxbins[i]-orPos[0]))
						oryindex = min(range(len(orybins)), key=lambda i: abs(orybins[i]-orPos[1]))
						orzindex = min(range(len(orzbins)), key=lambda i: abs(orzbins[i]-orPos[2]))
						data=scidata[ip][orzindex][oryindex][orxindex]
					else: 
						data=0
					
					flux += volPix*data/(cmDist+cmzPix*iz)**2.
					#end z
				image[ip][ix][iy]=flux
				#end y
			#end x
		print('Frequency (MHz): '+str(frequency/1.e+6))
		print('Total Luminosity (erg/s): '+str(np.sum(image[ip])*4.*math.pi*cmDist**2.*frequency))
		print('Total Flux (Jy): '+str(np.sum(image[ip])*1.0e+23))
		print('-----------------------------------------------------')
	#end p
			
	print('Integration complete...')
	print('=====================================================')
	
	#Write array to FITS image--------------------------------------------------  
	hduout = py.PrimaryHDU(beamConversion*image) # in ?/beam
	hduout.header['UNIT'  ] = 'erg cm^-2 sr^-1 s^-1 Hz^-1 hpbw'
	hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
	hduout.header['INCLN1'] = (inclination[0], 'Inclination angle 1')
	hduout.header['INCLN2'] = (inclination[1], 'Inclination angle 2')
	hduout.header['INCLN3'] = (inclination[2], 'Inclination angle 3')
	hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
	hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CRVAL3'] = (nmin, 'Minimum Frequency')
	hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CDELT3'] = (nfactor, 'frequency=CRVAL3*CDELT3^index')
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 
	print('synch_line_of_sight COMPLETED')

	############################################################################
	#############   END synch_line_of_sight  ###################################
	############################################################################
	



#=============================================================================
# gamma_line_of_sight(1,2,3,4)
# Parameters:
#		1) emissivity filename 
#		2) fits output filename 
#		3) inclination in degrees [phi,theta,psi]
#		4) physical distance to object (kpc)
#==============================================================================	
def gamma_line_of_sight(emissfilename, fileout, inclination, objectdist):
	import pyfits as py
	import numpy as np
	import math
	import os
	import sys
	from sys import stdout
	
	#constants
	pc2cm = 3.08568025e18 #cm/pc
	kpc2cm = 3.08568025e21 #cm/kpc
	radian2arcsec = 648000./math.pi #acsec/radian
	Sr2sqArcSec  = radian2arcsec**2.  # arcseconds^2/sr
	
	inc_rad = inclination*math.pi/180. #inclination in radians
	
	#READING FITS-------------------------------------------------
	hdulist=py.open(emissfilename)
	head=hdulist[0].header
	scidata=hdulist[0].data
	
	#Distances in kpc
	xmin = head['CRVAL1']
	ymin = head['CRVAL2']
	zmin = head['CRVAL3']
	log_e_min = head['CRVAL4'] #log10(e_min)
	
	xdel = head['CDELT1']
	ydel = head['CDELT2']
	zdel = head['CDELT3']
	log_efactor = head['CDELT4'] #factor
	
	xnum = head['NAXIS1']
	ynum = head['NAXIS2']
	znum = head['NAXIS3']
	enum = head['NAXIS4']
	
	emin = 10.**log_e_min  #MeV
	efactor = 10.**log_efactor
	
	xmax = xmin+xdel*xnum
	ymax = ymin+ydel*ynum
	zmax = zmin+zdel*znum
	
	xsize = xmax-xmin
	ysize = ymax-ymin
	zsize = zmax-zmin

	orPos = np.zeros(3) #original position (x,y,z) kpc	
	orxbins = np.zeros(xnum)
	orybins = np.zeros(ynum)
	orzbins = np.zeros(znum)
	for n in range(0,xnum):
		orxbins[n] = xmin+xdel*n
	for n in range(0,ynum):
		orybins[n] = ymin+ydel*n
	for n in range(0,znum):
		orzbins[n] = zmin+zdel*n
	
	original_mapsize = math.sqrt(xsize**2.+ysize**2.+zsize**2.) #kpc
	rotate_mapsize = 1.01*original_mapsize #kpc

	#OUTPUT IMAGE------------------------------------------------------
	outRes = 41 #number of pixels = outRes**2
	int_steps = 100 #number of integration steps
	
	kpcPix = rotate_mapsize/outRes #kpc/pixel_dim
	image = np.zeros([enum,outRes,outRes]) #output map (pixels)
	imPos = np.zeros(3) #image position [x,y,z] kpc
	
	xminIm = -rotate_mapsize/2. #kpc
	yminIm = -rotate_mapsize/2. #kpc
	zminIm = -rotate_mapsize/2. #kpc
	
	dxIm = rotate_mapsize/(outRes-1.) #kpc
	dyIm = rotate_mapsize/(outRes-1.) #kpc
	dzIm = rotate_mapsize/(int_steps-1.) #kpc
	
	cmxPix = dxIm*kpc2cm #cm/pixel
	cmyPix = dyIm*kpc2cm #cm/pixel
	cmzPix = dzIm*kpc2cm #cm/pixel
	
	cmDist = objectdist*kpc2cm
	solidangle = 1./((objectdist*kpc2cm)**2.) #4pi*ratio of 1cm^2 to entire sphere
	volPix = cmxPix*cmyPix*cmzPix #cm^3
	
	srPix = dxIm*dyIm/objectdist**2. #steradian/pixel
	Pix2sac = srPix*Sr2sqArcSec #arcsec^2/pixel
	sac2Pix = 1./Pix2sac #pixel/arsec^2
	
	#beamConversion = Pixsas*math.pi*(hpbw/2.)**2./math.log(2.)
	
	print('ArcSec^2/pixel: '+str(Pix2sac))
	#print('Beam Conversion: '+str(beamConversion))
	
	#INTEGRATION---------------------------------------------------------
	print('=====================================================')
	print('Beginning Integration')
	print('-----------------------------------------------------')
	
	for ip in range(0,enum):  # p
		energy=emin*efactor**ip
		
		for ix in range(0,outRes):  # x
			for iy in range(0,outRes):  # y
				flux = 0.											# total flux
				for iz in range(0,int_steps):  # z
					imPos[0] = xminIm+ix*dxIm
					imPos[1] = yminIm+iy*dyIm
					imPos[2] = zminIm+iz*dzIm
				
					orPos = rotate_coordinates(imPos,inc_rad[0],inc_rad[1],inc_rad[2])
				
					if (orPos[0] >= xmin-xdel/2. and orPos[0] <= xmax+xdel/2. and
							orPos[1] >= ymin-ydel/2. and orPos[1] <= ymax+ydel/2. and
							orPos[2] >= zmin-zdel/2. and orPos[2] <= zmax+zdel/2.):
						orxindex = min(range(len(orxbins)), key=lambda i: abs(orxbins[i]-orPos[0]))
						oryindex = min(range(len(orybins)), key=lambda i: abs(orybins[i]-orPos[1]))
						orzindex = min(range(len(orzbins)), key=lambda i: abs(orzbins[i]-orPos[2]))
						data=scidata[ip][orzindex][oryindex][orxindex]
					else: 
						data=0
					
					flux += volPix*data/(cmDist+cmzPix*iz)**2.
					#end z
				image[ip][ix][iy]=flux
				#end y
			#end x
		print('Energy (MeV): '+str(energy))
		print('Total Luminosity (erg/s): '+str(np.sum(image[ip])*4.*math.pi*cmDist**2.))
		print('Total Flux (MeV^2 cm^-2 sr^-1 s^-1 MeV^-1): '+str(np.sum(image[ip])))
		print('-----------------------------------------------------')
	#end p
			
	print('Integration complete...')
	print('=====================================================')
	
	#Write array to FITS image--------------------------------------------------  
	hduout = py.PrimaryHDU(image) # in ?/beam
	hduout.header['UNIT'  ] = 'MeV^2 cm^-2 sr^-1 s^-1 MeV^-1'
	hduout.header['OBJDIS'] = (objectdist, 'Distance to Object (kpc)')
	hduout.header['INCLN1'] = (inclination[0], 'Inclination angle 1')
	hduout.header['INCLN2'] = (inclination[1], 'Inclination angle 2')
	hduout.header['INCLN3'] = (inclination[2], 'Inclination angle 3')
	#hduout.header['HPBW'  ] = (hpbw, 'Half power beam width')
	hduout.header['CRVAL1'] = (xminIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CRVAL2'] = (yminIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CRVAL3'] = (emin, 'Minimum Energy (MeV)')
	hduout.header['CDELT1'] = (dxIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CDELT2'] = (dyIm/objectdist*radian2arcsec, 'arcseconds')
	hduout.header['CDELT3'] = (efactor, 'energy=CRVAL3*CDELT3^index')
	
	hduoutlist = py.HDUList([hduout])
	if os.access(fileout, os.F_OK ):  os.remove(fileout)
	hduoutlist.writeto(fileout)   
	print('FITS image output to '+str(fileout)) 
	print('gamma_line_of_sight COMPLETED')

	############################################################################
	#############   END gamma_line_of_sight  ###################################
	############################################################################



	
# Rotates a vector(x,y,z) by Euler angles
def rotate_coordinates(vector,phi,theta,psi):
    import math
    import numpy as np
    output = np.zeros(3)
    output[0] = vector[0]*math.cos(theta)*math.cos(psi) + vector[1]*math.cos(theta)*math.sin(psi) - vector[2]*math.sin(theta)
    output[1] = vector[0]*(-math.cos(phi)*math.sin(psi)+math.sin(phi)*math.sin(theta)*math.cos(psi)) + vector[1]*(math.cos(phi)*math.cos(psi)+math.sin(phi)*math.sin(theta)*math.sin(psi)) + vector[2]*math.sin(phi)*math.cos(theta)
    output[2] = vector[0]*(math.sin(phi)*math.sin(psi)+math.cos(phi)*math.sin(theta)*math.cos(psi))+vector[1]*(-math.sin(phi)*math.cos(psi)+math.cos(phi)*math.sin(theta)*math.sin(psi))+vector[2]*math.cos(phi)*math.cos(theta)
    return output

    
