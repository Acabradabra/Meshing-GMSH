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
Lt= 0.5
Lc=2e-2
Ls=1e-2
ep= 1.0e-3

#=====> Mesh Size
h0=1e-4
h1=1e-4
h2=2e-4
hf=2e-3

hfc=1e-3
hfs=1e-3

#=====> Flame refinement params
Lr=5e-2
Lf=2*ep
Dr= 5e-3

#=====> Boundary layer params
rb=1.1
rf=1
ra=0.7
N0=10
N1=10
N2=10

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
(h00,Lt0)=me.h_smooth(rb,N0,rf*h0) ; print('=> Main  jet    : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h00,Lt0,ra*h0/h00))
(h10,Lt1)=me.h_smooth(rb,N1,rf*h1) ; print('=> Main  Pilot  : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h10,Lt1,ra*h1/h10))
(h20,Lt2)=me.h_smooth(rb,N2,rf*h2) ; print('=> Main  Coflow : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h20,Lt2,ra*h2/h20))
N0z=int(round(Lc/(h0*ra),0))+1
N1z=int(round(Lc/(h1*ra),0))+1
N2z=int(round(Ls/(h2*ra),0))+1

N0r=int(round( (0.5*D0        -  Lt0)/h0 ,0))+1
N1r=int(round( (0.5*D1-(r0+rp)-2*Lt1)/h1 ,0))+1
N2r=int(round( (r2    -(r1+rp)-  Lt2)/h2 ,0))+1

if BD : sys.exit('=> Boundary layer')

#=====> Y dim
y1=Lc-Ls
y2=Lc
y3=Lt+Lc

gm.initialize()
#=====> Points
Points=[
gm.model.geo.add_point(0        ,y2,0,meshSize=h0 ), #===> Main Jet
gm.model.geo.add_point(0        , 0,0,meshSize=h0 ),
gm.model.geo.add_point(r0-rp-Lt0, 0,0,meshSize=h0 ),
gm.model.geo.add_point(r0-rp    , 0,0,meshSize=h00),
gm.model.geo.add_point(r0-rp    ,y2,0,meshSize=h00),
gm.model.geo.add_point(r0-rp-Lt0,y2,0,meshSize=h0 ),
gm.model.geo.add_point(r0+rp    ,y2,0,meshSize=h10), #===> Pilot Jet
gm.model.geo.add_point(r0+rp+Lt1,y2,0,meshSize=h1 ), 
gm.model.geo.add_point(r0+rp    , 0,0,meshSize=h10),
gm.model.geo.add_point(r0+rp+Lt1, 0,0,meshSize=h1 ),
gm.model.geo.add_point(r1-rp-Lt1, 0,0,meshSize=h1 ),
gm.model.geo.add_point(r1-rp    , 0,0,meshSize=h10),
gm.model.geo.add_point(r1-rp    ,y2,0,meshSize=h10),
gm.model.geo.add_point(r1-rp-Lt1,y2,0,meshSize=h1 ),
gm.model.geo.add_point(r1+rp    ,y2,0,meshSize=h20), #===> Coflow
gm.model.geo.add_point(r1+rp+Lt2,y2,0,meshSize=h2 ),
gm.model.geo.add_point(r1+rp    ,y1,0,meshSize=h20),
gm.model.geo.add_point(r1+rp+Lt2,y1,0,meshSize=h2 ),
gm.model.geo.add_point(r2       ,y1,0,meshSize=h2 ),
gm.model.geo.add_point(r2       ,y2,0,meshSize=h2 ),
gm.model.geo.add_point(r2       ,y3,0,meshSize=hf ), #===> Outlet
gm.model.geo.add_point(0        ,y3,0,meshSize=hf ) ]
Np=len(Points)
gm.model.geo.synchronize()

#=====> Lines
lml=gm.model.geo.add_line(Points[ 0],Points[ 1]) #===> Main left
lmc=gm.model.geo.add_line(Points[ 1],Points[ 2]) #===> Main Inlet flow
lm1=gm.model.geo.add_line(Points[ 2],Points[ 3]) #===> Main Inlet BL
lmr=gm.model.geo.add_line(Points[ 3],Points[ 4]) #===> Main right
lm2=gm.model.geo.add_line(Points[ 4],Points[ 5]) #===> Main Outlet BL
lmo=gm.model.geo.add_line(Points[ 5],Points[ 0]) #===> Main Outlet flow
lme=gm.model.geo.add_line(Points[ 4],Points[ 6]) #===> Main leap
lpl=gm.model.geo.add_line(Points[ 6],Points[ 8]) #===> Pilot left
lp2=gm.model.geo.add_line(Points[ 6],Points[ 7]) #===> Pilot Outlet BL L
lp1=gm.model.geo.add_line(Points[ 8],Points[ 9]) #===> Pilot inlet BL L
lpc=gm.model.geo.add_line(Points[ 9],Points[10]) #===> Pilot inlet flow
lp3=gm.model.geo.add_line(Points[10],Points[11]) #===> Pilot inlet BL R
lpr=gm.model.geo.add_line(Points[11],Points[12]) #===> Pilot right
lp4=gm.model.geo.add_line(Points[12],Points[13]) #===> Pilot Outlet BL R
lpo=gm.model.geo.add_line(Points[13],Points[ 7]) #===> Pilot Outlet
lpe=gm.model.geo.add_line(Points[12],Points[14]) #===> Pilot leap
lc2=gm.model.geo.add_line(Points[14],Points[15]) #===> Coflow Outlet BL 
lcl=gm.model.geo.add_line(Points[14],Points[16]) #===> Coflow Left
lc1=gm.model.geo.add_line(Points[16],Points[17]) #===> Coflow Inlet BL
lcc=gm.model.geo.add_line(Points[17],Points[18]) #===> Coflow Inlet
lcr=gm.model.geo.add_line(Points[18],Points[19]) #===> Coflow right
lco=gm.model.geo.add_line(Points[19],Points[15]) #===> Coflow outlet
lor=gm.model.geo.add_line(Points[19],Points[20]) #===> Outlet right
loo=gm.model.geo.add_line(Points[20],Points[21]) #===> Outlet
lol=gm.model.geo.add_line(Points[21],Points[ 0]) #===> Outlet left
gm.model.geo.synchronize()

#=====> transfinite
gm.model.geo.mesh.setTransfiniteCurve(lm1,N0+1,meshType='Progression',coef=1/rb) #===> Main Jet
gm.model.geo.mesh.setTransfiniteCurve(lm2,N0+1,meshType='Progression',coef=  rb)
gm.model.geo.mesh.setTransfiniteCurve(lmc,N0r)
gm.model.geo.mesh.setTransfiniteCurve(lmo,N0r)
gm.model.geo.mesh.setTransfiniteCurve(lmr,N0z)
gm.model.geo.mesh.setTransfiniteCurve(lml,N0z)
gm.model.geo.mesh.setTransfiniteCurve(lp2,N1+1,meshType='Progression',coef=rb) #===> Pilot
gm.model.geo.mesh.setTransfiniteCurve(lp1,N1+1,meshType='Progression',coef=rb)
gm.model.geo.mesh.setTransfiniteCurve(lpr,N1z)
gm.model.geo.mesh.setTransfiniteCurve(lpl,N1z)
gm.model.geo.mesh.setTransfiniteCurve(lp4,N1+1,meshType='Progression',coef=  rb)
gm.model.geo.mesh.setTransfiniteCurve(lp3,N1+1,meshType='Progression',coef=1/rb)
gm.model.geo.mesh.setTransfiniteCurve(lpc,N1r)
gm.model.geo.mesh.setTransfiniteCurve(lpo,N1r)
gm.model.geo.mesh.setTransfiniteCurve(lc1,N2+1,meshType='Progression',coef=rb) #===> Coflow
gm.model.geo.mesh.setTransfiniteCurve(lc2,N2+1,meshType='Progression',coef=rb)
gm.model.geo.mesh.setTransfiniteCurve(lcc,N2r)
gm.model.geo.mesh.setTransfiniteCurve(lco,N2r)
gm.model.geo.mesh.setTransfiniteCurve(lcr,N2z)
gm.model.geo.mesh.setTransfiniteCurve(lcl,N2z)
gm.model.geo.synchronize()

#=====> Surfaces
lm=gm.model.geo.add_curve_loop([lml,lmc,lm1,lmr,lm2,lmo]) #===> Main Jet
sm=gm.model.geo.add_plane_surface([lm])
lp=gm.model.geo.add_curve_loop([lpl,lp1,lpc,lp3,lpr,lp4,lpo,-lp2]) #===> Pilot Jet
sp=gm.model.geo.add_plane_surface([lp])
lc=gm.model.geo.add_curve_loop([lcl,lc1,lcc,lcr,lco,-lc2]) #===> Coflow
sc=gm.model.geo.add_plane_surface([lc])
lo=gm.model.geo.add_curve_loop([lol,-lmo,-lm2,lme,lp2,-lpo,-lp4,lpe,lc2,-lco,lor,loo]) #===> Outlet
so=gm.model.geo.add_plane_surface([lo])
gm.model.geo.synchronize()

#=====> transfinite
gm.model.geo.mesh.setTransfiniteSurface(sm,cornerTags=[Points[ 0],Points[ 1],Points[ 3],Points[ 4]])
gm.model.geo.mesh.setTransfiniteSurface(sp,cornerTags=[Points[ 6],Points[ 8],Points[11],Points[12]])
gm.model.geo.mesh.setTransfiniteSurface(sc,cornerTags=[Points[14],Points[16],Points[18],Points[19]])
gm.model.geo.synchronize()

#=====> Recombine
gm.model.geo.mesh.setRecombine(gdim,sm)
gm.model.geo.mesh.setRecombine(gdim,sp)
gm.model.geo.mesh.setRecombine(gdim,sc)
gm.model.geo.synchronize()

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
gm.model.mesh.field.add("Box", 3) # Flame side
gm.model.mesh.field.setNumber(3, "VIn" , h00 )
gm.model.mesh.field.setNumber(3, "VOut", hf  )
gm.model.mesh.field.setNumber(3, "XMin", r0-rp)
gm.model.mesh.field.setNumber(3, "XMax", r0+rp)
gm.model.mesh.field.setNumber(3, "YMin", y2   )
gm.model.mesh.field.setNumber(3, "YMax", y2+Lf)
gm.model.mesh.field.setNumber(3, "Thickness",0.2)
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
gm.model.mesh.field.setNumbers(4, "FieldsList", [1,2,3])
gm.model.mesh.field.setAsBackgroundMesh(4)

#===================================================================================
util.Section( 'Physical groups : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.model.addPhysicalGroup(1,[lml,lol]    ,tag=1000,name='Axis'        )
gm.model.addPhysicalGroup(1,[lmc,lm1]    ,tag=2001,name='Inlet-Main'  )
gm.model.addPhysicalGroup(1,[lpc,lp1,lp3],tag=2002,name='Inlet-Pilot' )
gm.model.addPhysicalGroup(1,[lcc,lc1]    ,tag=2003,name='Inlet-Coflow')
gm.model.addPhysicalGroup(1,[lmr]        ,tag=3001,name='Wall-Main'   )
gm.model.addPhysicalGroup(1,[lpl,lpr]    ,tag=3002,name='Wall-Pilot'  )
gm.model.addPhysicalGroup(1,[lcl]        ,tag=3003,name='Wall-Coflow' )
gm.model.addPhysicalGroup(1,[lme]        ,tag=4001,name='Leap-Main'   )
gm.model.addPhysicalGroup(1,[lpe]        ,tag=4002,name='Leap-Pilot'  )
gm.model.addPhysicalGroup(1,[lcr,lor]    ,tag=5001,name='External'    )
gm.model.addPhysicalGroup(1,[loo]        ,tag=5002,name='Outlet'      )
gm.model.addPhysicalGroup(2,[sm,sp,sc,so],tag=6000,name='Fluid'       )
# gm.model.addPhysicalGroup(2,[sm]         ,tag=6001,name='Main'        )
# gm.model.addPhysicalGroup(2,[sp]         ,tag=6002,name='Pilot'       )
# gm.model.addPhysicalGroup(2,[sc]         ,tag=6003,name='Coflow'      )
# gm.model.addPhysicalGroup(2,[so]         ,tag=6004,name='Flame'       )

#===================================================================================
util.Section( 'Meshing : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.option.setNumber("Mesh.Algorithm", 7)
gm.option.setNumber("Mesh.SmoothRatio",1)
gm.option.setNumber('Mesh.OptimizeThreshold',0.6)
gm.option.setNumber('Mesh.OptimizeNetgen',1)
gm.option.setNumber("Mesh.SaveElementTagType",2)
gm.option.setNumber("Mesh.MeshSizeExtendFromBoundary",0)
gm.option.setNumber("Mesh.MeshSizeFromPoints",0)
gm.option.setNumber("Mesh.MeshSizeFromCurvature",0)
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

gm.finalize()

if CONV :
	#===================================================================================
	util.Section( 'Conversion : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	mesh=mo.gmsh.read(d_mesh+name+'-MESH.msh')
	mo.ansys.write(d_mesh+name+'-ANSYS.msh',mesh)
	# mo.cgns.write(d_mesh+name+'-MESH.cgns',mesh)