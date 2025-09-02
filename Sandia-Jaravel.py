#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['BD','Conv'],0,'Arg : ')
(                             [ BD , CONV ])=Arg

from numpy import *
import sys
import time
import pygmsh as pg
import gmsh   as gm
import Meshing as me
import meshio as mo

import matplotlib.pyplot as plt
import matplotlib

t0=time.time()
#===================================================================================
#                     Parameters
#===================================================================================

gdim = 2  # Geometric dimension of the mesh

#=====> files
d_mesh='MESH/'
d_plot='PLOT/'
d_data='DATA/'
name='Sandia-Jaravel'

#=====> Dimensions
D0= 7.2e-3
D1=18.2e-3
D2=50.0e-3
Lt= 2
Lc=2e-2
Ls=1e-2
ep= 1.0e-3

#=====> Mesh Size
h0=1e-3
h1=1e-3
h2=1e-3
hf=1e-2

hfc=1e-3
hfs=1e-3

#=====> Flame refinement params
Lr=10e-2
Dr= 5e-3

#=====> Boundary layer params
rbd=1.1
rb0=0.5
rbz=0.75

if False :
	#===================================================================================
	util.Section( 'Conversion : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	d_fluent='/mnt/d/FLUENT/Sandia-Jaravel/MESH/'
	# f_in=d_fluent+'burner-axisymmetric-2d.msh'
	# f_ou=d_fluent+'burner-axisymmetric-2d-ASCII.msh'

	mesh=mo.gmsh.read(d_mesh+name+'-MESH.msh')
	# mesh=mo.ansys.read(f_in)
	# mo.ansys.write(f_ou,binary=False)
	mo.ansys.write(d_mesh+name+'-ANSYS.msh',mesh,binary=False)
	# mo.cgns.write(d_mesh+name+'-MESH.cgns',mesh)
	sys.exit('=> File converted')

#===================================================================================
util.Section( 'Geometrie : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

#=====> X dim
rp=0.5*ep
r0=0.5*D0+rp
r1=0.5*D1+rp
r2=0.5*D2

#=====> Boundary layer
Lbm=     r0- rp        
Lbp=0.5*(r1-(rp+r0+ep))
Lbc=     r2-(rp+r1   ) 
Nmr=int(me.Nbl(rbd,Lbm,rb0*h0))+1 ; him=me.h0(Lbm,rbd,Nmr) ; hnm=him*rbd**(Nmr-1) ; print('=> Main  jet : h {:.2e} , h0 {:.2e} , N : {} , hn {:.2e}'.format(h0,him,Nmr,hnm))
Npr=int(me.Nbl(rbd,Lbp,rb0*h1))+1 ; hip=me.h0(Lbp,rbd,Npr) ; hnp=hip*rbd**(Npr-1) ; print('=> Pilot jet : h {:.2e} , h0 {:.2e} , N : {} , hn {:.2e}'.format(h1,hip,Npr,hnp)) ; rb=hip/hnp
Ncr=int(me.Nbl(rbd,Lbc,rb0*h2))+1 ; hic=me.h0(Lbc,rbd,Ncr) ; hnc=hic*rbd**(Ncr-1) ; print('=> Coflow    : h {:.2e} , h0 {:.2e} , N : {} , hn {:.2e}'.format(h2,hic,Ncr,hnc))
Nmz=int(round(Lc/(hnm*rbz),0))
Npz=int(round(Lc/(hnp*rbz),0))
Ncz=int(round(Ls/(hnc*rbz),0))

if BD : sys.exit('=> Boundary layer')

#=====> Y dim
y1=Lc-Ls
y2=Lc
y3=Lt+Lc

gm.initialize()
#=====> Points
Points=[
gm.model.geo.add_point(0    ,y2,0,meshSize=hnm), #===> Main Jet
gm.model.geo.add_point(0    , 0,0,meshSize=h0),
gm.model.geo.add_point(r0-rp, 0,0,meshSize=h0),
gm.model.geo.add_point(r0-rp,y2,0,meshSize=him),
gm.model.geo.add_point(r0+rp,y2,0,meshSize=hip), #===> Pilot Jet
gm.model.geo.add_point(r0+rp, 0,0,meshSize=h1),
gm.model.geo.add_point(r1-rp, 0,0,meshSize=h1),
gm.model.geo.add_point(r1-rp,y2,0,meshSize=hip),
gm.model.geo.add_point(r1+rp,y2,0,meshSize=hic), #===> Coflow
gm.model.geo.add_point(r1+rp,y1,0,meshSize=h2),
gm.model.geo.add_point(r2   ,y1,0,meshSize=h2),
gm.model.geo.add_point(r2   ,y2,0,meshSize=hnc),
gm.model.geo.add_point(r2   ,y3,0,meshSize=hf), #===> Outlet
gm.model.geo.add_point(0    ,y3,0,meshSize=hf) ]
Np=len(Points)
gm.model.geo.synchronize()

#=====> Lines
lml=gm.model.geo.add_line(Points[ 0],Points[ 1]) #===> Main left
lmc=gm.model.geo.add_line(Points[ 1],Points[ 2]) #===> Main Inlet
lmr=gm.model.geo.add_line(Points[ 2],Points[ 3]) #===> Main right
lme=gm.model.geo.add_line(Points[ 3],Points[ 4]) #===> Main leap
lmo=gm.model.geo.add_line(Points[ 3],Points[ 0]) #===> Main outlet
lpl=gm.model.geo.add_line(Points[ 4],Points[ 5]) #===> Pilot left
lpc=gm.model.geo.add_line(Points[ 5],Points[ 6]) #===> Pilot Inlet
lpr=gm.model.geo.add_line(Points[ 6],Points[ 7]) #===> Pilot right
lpe=gm.model.geo.add_line(Points[ 7],Points[ 8]) #===> Pilot leap
lpo=gm.model.geo.add_line(Points[ 7],Points[ 4]) #===> Pilot outlet
lcl=gm.model.geo.add_line(Points[ 8],Points[ 9]) #===> Coflow left
lcc=gm.model.geo.add_line(Points[ 9],Points[10]) #===> Coflow Inlet
lcr=gm.model.geo.add_line(Points[10],Points[11]) #===> Coflow right
lco=gm.model.geo.add_line(Points[11],Points[ 8]) #===> Coflow outlet
lor=gm.model.geo.add_line(Points[11],Points[12]) #===> Outlet right
lou=gm.model.geo.add_line(Points[12],Points[13]) #===> Outlet
lol=gm.model.geo.add_line(Points[13],Points[0 ]) #===> Outlet left
gm.model.geo.synchronize()

#=====> transfinite
# gm.model.geo.mesh.setTransfiniteCurve(lmc,Nmr+1,meshType='Progression',coef=1/rbd) #===> Main Jet
# gm.model.geo.mesh.setTransfiniteCurve(lmo,Nmr+1,meshType='Progression',coef=rbd  )
# gm.model.geo.mesh.setTransfiniteCurve(lmr,Nmz+1)
# gm.model.geo.mesh.setTransfiniteCurve(lml,Nmz+1)
# gm.model.geo.mesh.setTransfiniteCurve(lcc,Ncr+1,meshType='Progression',coef=rbd  ) #===> Coflow
# gm.model.geo.mesh.setTransfiniteCurve(lco,Ncr+1,meshType='Progression',coef=1/rbd)
# gm.model.geo.mesh.setTransfiniteCurve(lcr,Ncz+1)
# gm.model.geo.mesh.setTransfiniteCurve(lcl,Ncz+1)
# gm.model.geo.mesh.setTransfiniteCurve(lpc,2*Npr+1,meshType='Bump',coef=rb) #===> Pilot
# gm.model.geo.mesh.setTransfiniteCurve(lpo,2*Npr+1,meshType='Bump',coef=rb)
# gm.model.geo.mesh.setTransfiniteCurve(lpr,  Npz+1)
# gm.model.geo.mesh.setTransfiniteCurve(lpl,  Npz+1)
# gm.model.geo.synchronize()

#=====> Surfaces
lm=gm.model.geo.add_curve_loop([lml,lmc,lmr,lmo]) #===> Main Jet
sm=gm.model.geo.add_plane_surface([lm])
lp=gm.model.geo.add_curve_loop([lpl,lpc,lpr,lpo]) #===> Pilot Jet
sp=gm.model.geo.add_plane_surface([lp])
lc=gm.model.geo.add_curve_loop([lcl,lcc,lcr,lco]) #===> Coflow
sc=gm.model.geo.add_plane_surface([lc])
lo=gm.model.geo.add_curve_loop([lol,-lmo,lme,-lpo,lpe,-lco,lor,lou]) #===> Outlet
so=gm.model.geo.add_plane_surface([lo])
gm.model.geo.synchronize()

#=====> transfinite
# gm.model.geo.mesh.setTransfiniteSurface(sm)
# gm.model.geo.mesh.setTransfiniteSurface(sp)
# gm.model.geo.mesh.setTransfiniteSurface(sc)
# gm.model.geo.synchronize()

#=====> Recombine
# gm.model.geo.mesh.setRecombine(gdim,sm)
# gm.model.geo.mesh.setRecombine(gdim,sp)
# gm.model.geo.mesh.setRecombine(gdim,sc)
# gm.model.geo.synchronize()

#=====> Field
gm.model.mesh.field.add("Box", 1) # Flame center
gm.model.mesh.field.setNumber(1, "VIn" , hfc )
gm.model.mesh.field.setNumber(1, "VOut", hf  )
gm.model.mesh.field.setNumber(1, "XMin", 0  )
gm.model.mesh.field.setNumber(1, "XMax", r0 )
gm.model.mesh.field.setNumber(1, "YMin", y2   )
gm.model.mesh.field.setNumber(1, "YMax", y2+Lr)
gm.model.mesh.field.setNumber(1, "Thickness",  1e-2)
gm.model.mesh.field.add("Box", 2) # Flame side
gm.model.mesh.field.setNumber(2, "VIn" , hfs )
gm.model.mesh.field.setNumber(2, "VOut", hf  )
gm.model.mesh.field.setNumber(2, "XMin", r0   )
gm.model.mesh.field.setNumber(2, "XMax", r0+Dr)
gm.model.mesh.field.setNumber(2, "YMin", y2   )
gm.model.mesh.field.setNumber(2, "YMax", y2+Lr)
gm.model.mesh.field.setNumber(2, "Thickness",  1e-2)
# gm.model.mesh.field.add("BoundaryLayer", 3)
# gm.model.mesh.field.setNumber(3, "AnisoMax"  , 1e3 )
# gm.model.mesh.field.setNumber(3, "Quads"     , 1   )
# gm.model.mesh.field.setNumber(3, "Thickness" , 1   )
# gm.model.mesh.field.setNumber(3, "CurvesList", lmr )
# gm.model.mesh.field.setNumber(3, "NbLayers"  , 5   )
# gm.model.mesh.field.setNumber(3, "Ratio"     , 1.1 )
# gm.model.mesh.field.setNumber(3, "Size"      , h0*rb0 )
# gm.model.mesh.field.setNumber(3, "SizeFar"   , h0  )
gm.model.mesh.field.add("Min", 4)
gm.model.mesh.field.setNumbers(4, "FieldsList", [1,2])
gm.model.mesh.field.setAsBackgroundMesh(4)

#===================================================================================
util.Section( 'Physical groups : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.model.addPhysicalGroup(1,[lml,lol],tag=1000,name='Axis'        )
gm.model.addPhysicalGroup(1,[lmc]    ,tag=2001,name='Inlet-Main'  )
gm.model.addPhysicalGroup(1,[lpc]    ,tag=2002,name='Inlet-Pilot' )
gm.model.addPhysicalGroup(1,[lcc]    ,tag=2003,name='Inlet-Coflow')
gm.model.addPhysicalGroup(1,[lmr]    ,tag=3001,name='Wall-Main'   )
gm.model.addPhysicalGroup(1,[lpl,lpr],tag=3002,name='Wall-Pilot'  )
gm.model.addPhysicalGroup(1,[lcl]    ,tag=3003,name='Wall-Coflow' )
gm.model.addPhysicalGroup(1,[lme]    ,tag=4001,name='Leap-Main'   )
gm.model.addPhysicalGroup(1,[lpe]    ,tag=4002,name='Leap-Pilot'  )
gm.model.addPhysicalGroup(1,[lcr,lor],tag=5001,name='External'    )
gm.model.addPhysicalGroup(1,[lou]    ,tag=5002,name='Outlet'      )
gm.model.addPhysicalGroup(2,[sm]     ,tag=6001,name='Main'        )
gm.model.addPhysicalGroup(2,[sp]     ,tag=6002,name='Pilot'       )
gm.model.addPhysicalGroup(2,[sc]     ,tag=6003,name='Coflow'      )
gm.model.addPhysicalGroup(2,[so]     ,tag=6004,name='Flame'       )

#===================================================================================
util.Section( 'Meshing : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.option.setNumber("Mesh.Algorithm", 6)
gm.option.setNumber("Mesh.SmoothRatio",1)
gm.option.setNumber("Mesh.SaveElementTagType",2)
gm.model.mesh.setSmoothing(gdim,so,100)
gm.model.mesh.generate(gdim)
gm.model.mesh.optimize("Netgen")

#===================================================================================
util.Section( 'Writing : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.model.geo.synchronize()
gm.write(d_mesh+name+'-Geo.geo_unrolled')
gm.write(d_mesh+name+'-MESH.msh')
gm.write(d_mesh+name+'-MESH.bdf')
# gm.write(d_mesh+name+'-MESH.unv')
# gm.write(d_mesh+name+'-MESH.cgns')

gm.finalize()

if CONV :
	#===================================================================================
	util.Section( 'Conversion : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	mesh=mo.gmsh.read(d_mesh+name+'-MESH.msh')
	mo.ansys.write(d_mesh+name+'-ANSYS.msh',mesh)
	# mo.cgns.write(d_mesh+name+'-MESH.cgns',mesh)