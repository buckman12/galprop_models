import pyfits as py
import numpy as np
import math

r_0 = 0.200 #kpc
h_0 = 0.050 #kpc

kpc2cm = 3.086e+18*1.e+3 #cm/kpc
m_H1 = 1.67e-24 #g

sigma = 0.25 #g*cm^-2
rho = sigma/(2*h_0*kpc2cm) #g*cm^-3
H1 = rho/m_H1 #cm^-3

xmin = -2.
xmax = +2.
ymin = -2.
ymax = +2.
zmin = -1.0
zmax = +1.0

xnum = 401
ynum = 401
znum = 401

xdel = (xmax-xmin)/(xnum-1)
ydel = (ymax-ymin)/(ynum-1)
zdel = (zmax-zmin)/(znum-1)

gas_map=np.zeros((znum,ynum,xnum), dtype=np.float32)

#x=xmin
#for i in range(0,xnum):
#	y=ymin
#	for j in range(0,ynum):
#		z=zmin
#		for k in range(0,znum):
#			if x**2+y**2<r_0**2 and abs(z)<h_0:
#				gas_map[k][j][i]=H1
#			else:
#				gas_map[k][j][i]=H1*math.exp(-math.sqrt(x**2+y**2)/r_0)*math.exp(-abs(z)/h_0)
#			z+=zdel
#		y+=ydel
#	x+=xdel
	
hdu=py.PrimaryHDU(gas_map)

hdu.header['CRVAL1']=xmin
hdu.header['CRVAL2']=ymin
hdu.header['CRVAL3']=zmin

hdu.header['CDELT1']=xdel
hdu.header['CDELT2']=ydel
hdu.header['CDELT3']=zdel

hdu.header['CTYPE1']='x in kpc'
hdu.header['CTYPE2']='y in kpc'
hdu.header['CTYPE3']='z in kpc'

hdulist=py.HDUList(hdu)
hdulist.writeto('CO_m82_torus_v1.fits')


