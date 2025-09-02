#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import os
import sys
############################################################################################################################################
#========================>     Parseur
############################################################################################################################################
def Parseur(Varg,n0,strEXIT) :
	Sysa=sys.argv[1:]
	Narg=len(Varg)
	Arg,Arg0=[],''
	for arg in Varg :
		arg='--'+arg ; Arg0+=arg+' '
		if arg in Sysa : Arg.append(True) ; Sysa.pop(Sysa.index(arg))
		else           : Arg.append(False)
	NSysa=len(Sysa)
	exit=('--help' in Sysa or NSysa<n0)
	# if   exit and strEXIT : sys.exit(strEXIT)
	if   exit : sys.exit('\n\n '+strEXIT+'\nOpt : '+Arg0+'\n\n')
	# if   exit and strEXIT : sys.exit('\n\n '+strEXIT+'  ,  ARG : '+Arg0+'\n\n')
	# elif exit             : sys.exit('\n\n ARG : '+Arg0+'\n\n')
	Section('Arg (Sysa,Varg,Arg) :',1,5,'b')
	print(Sysa)
	print(Varg)
	print(Arg)
	print('\n')
	return(Sysa,NSysa,Arg)
#---------------------------------------------------------------------
def Section(txt,N0,N1,col) : 
	print(
		N0*'\n'
		+N1*'='+'> '
		+'\033[1;31m'*(col=='r')
		+'\033[1;33m'*(col=='y')
		+'\033[1;32m'*(col=='g')
		+'\033[1;34m'*(col=='b')
		+txt+'\033[0m'
		+N0*'\n'
	)
############################################################################################################################################
txtEXIT='arg : I \n\n options : \n\
--Plot  : Exporte les graphiques \n\
--Data  : Exporte le champ magnétique dans la plaque \n\
--Mesh  : Force remeshing (otherwise check for mesh file) \n\
--Clear : Erase all plots and data \n\
--All   : Plot + Data '
(Sysa,NSysa,Arg)=Parseur(['Plot','Data','Mesh','Clear','All'],0,txtEXIT)
[                          PLOT , DATA , MESH , CLEAR , ALL ]=Arg
if ALL : 
	PLOT=True
	DATA=True
############################################################################################################################################
#========================>     Librairies
############################################################################################################################################
from dolfinx import default_scalar_type
from dolfinx.fem import (dirichletbc, Expression, Function, functionspace, locate_dofs_topological)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile
from dolfinx.io.gmshio import model_to_mesh
from dolfinx.mesh import compute_midpoints, locate_entities_boundary
from dolfinx.plot import vtk_mesh

from basix.ufl import element, mixed_element
from ufl import TestFunction, TrialFunction, as_vector, dot, dx, grad, inner, SpatialCoordinate
from mpi4py import MPI
from scipy.special import jv
from petsc4py import PETSc

import gmsh
import numpy as np
import pyvista
import time

real_type = PETSc.RealType  # type: ignore
t0=time.time()
############################################################################################################################################
#========================>     Parameters
############################################################################################################################################
gdim = 2  # Geometric dimension of the mesh

#=====> files
d_mesh='MESH/'
d_plot='PLOT/'
d_data='DATA/'
f_mesh=d_mesh+'Mesh.msh'

#=====> Geometry
r0=0 #5e-3
R=2e-2
rw=1e-3
Ly=5e-2
Lx=10e-2

#=====> Plate
Dz=2e-2
Rp=5e-2
Ep=1e-2

#=====> Electrycity
I=1
# I=float(Sysa[0])

#=====> tags
t_air=0
t_Itp=1
t_Itm=2
t_pla=3

#=====> Mesh Size
h0=1e-4
h1=1e-3

############################################################################################################################################
print('\033[1;31m=====> Start Mesh : {:.3f} s\033[0m'.format(time.time()-t0))
############################################################################################################################################

gmsh.initialize()
model_rank = 0
mesh_comm = MPI.COMM_WORLD
if mesh_comm.rank == model_rank and (not os.path.exists(f_mesh) or MESH ) :

	#=====> Background
	Back = gmsh.model.occ.addRectangle(r0,-Ly,0,Lx-r0,2*Ly)
	gmsh.model.occ.synchronize()

	#=====> Wires
	Wires=[]
	# Wire_L = gmsh.model.occ.addDisk(-R, 0, 0, rw,rw) ; Wires.append((2,Wire_L))
	Wire_R = gmsh.model.occ.addDisk( R, 0, 0, rw,rw) ; Wires.append((2,Wire_R))
	gmsh.model.occ.synchronize()

	#=====> Metal plate
	Plate  = gmsh.model.occ.addRectangle(r0,Dz,0,Rp-r0,Ep)
	gmsh.model.occ.synchronize()

	#=====> Domain
	whole_domain = gmsh.model.occ.fragment( [(2,Back)],Wires+[(2,Plate)])
	gmsh.model.occ.synchronize()

	#=====> Physical groups
	gmsh.model.addPhysicalGroup(2,[Plate+1], tag=t_air) # Generalized
	# gmsh.model.addPhysicalGroup(2,[Wire_L], tag=t_Itp)
	gmsh.model.addPhysicalGroup(2,[Wire_R], tag=t_Itm)
	gmsh.model.addPhysicalGroup(2,[Plate ], tag=t_pla)

	#=====> Field
	edges = gmsh.model.getBoundary(Wires, oriented=False)
	gmsh.model.mesh.field.add("Distance", 1)
	gmsh.model.mesh.field.setNumbers(1, "EdgesList", [e[1] for e in edges])
	gmsh.model.mesh.field.add("Threshold", 2)
	gmsh.model.mesh.field.setNumber(2, "IField", 1)
	gmsh.model.mesh.field.setNumber(2, "LcMin", h0)
	gmsh.model.mesh.field.setNumber(2, "LcMax", h1)
	gmsh.model.mesh.field.setNumber(2, "DistMin", 4 * rw)
	gmsh.model.mesh.field.setNumber(2, "DistMax", 10 * rw)
	gmsh.model.mesh.field.add("Box", 3) # Plate
	gmsh.model.mesh.field.setNumber(3, "VIn" , h0 )
	gmsh.model.mesh.field.setNumber(3, "VOut", h1 )
	gmsh.model.mesh.field.setNumber(3, "XMin",  r0)
	gmsh.model.mesh.field.setNumber(3, "XMax",  Rp)
	gmsh.model.mesh.field.setNumber(3, "YMin",  Dz)
	gmsh.model.mesh.field.setNumber(3, "YMax",  Dz+Ep)
	gmsh.model.mesh.field.setNumber(3, "Thickness",  0.01)
	gmsh.model.mesh.field.add("Box", 4) # Central axis
	gmsh.model.mesh.field.setNumber(4, "VIn" , h0 )
	gmsh.model.mesh.field.setNumber(4, "VOut", h1 )
	gmsh.model.mesh.field.setNumber(4, "XMin",  r0)
	gmsh.model.mesh.field.setNumber(4, "XMax",    R)
	gmsh.model.mesh.field.setNumber(4, "YMin",  -Ly)
	gmsh.model.mesh.field.setNumber(4, "YMax",   Ly)
	gmsh.model.mesh.field.setNumber(4, "Thickness",  0.01)
	gmsh.model.mesh.field.add("Min", 5)
	gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2])
	gmsh.model.mesh.field.setAsBackgroundMesh(5)

	#=====> Meshing
	gmsh.option.setNumber("Mesh.Algorithm", 7)
	gmsh.model.mesh.generate(gdim)
	gmsh.model.mesh.optimize("Netgen")

	#=====> Export
	gmsh.write(f_mesh)
	gmsh.write(d_mesh+'Geo.geo_unrolled')
	# sys.exit()
	mesh, ct, _ = model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)
	gmsh.finalize()
	#========================> Paraview
	# print('===> Paraview export')
	# with XDMFFile(MPI.COMM_WORLD, d_mesh+"mt.xdmf", "w") as xdmf:
	#     xdmf.write_mesh(mesh)
	#     xdmf.write_meshtags(ct, mesh.geometry)
elif mesh_comm.rank == model_rank :
	gmsh.model.add("Mesh from file")
	gmsh.merge(f_mesh)
	mesh, ct, _ = model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)
	# output = model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
	gmsh.finalize()
############################################################################################################################################
print('\033[1;31m=====> Start Resolution : {:.3f} s\033[0m'.format(time.time()-t0))
############################################################################################################################################

#========================> Parameter
Q  = functionspace(mesh, ("DG", 0))
mu = Function(Q)
J  = Function(Q)
J.x.array[:] = 0
material_tags = np.unique(ct.values)
for tag in material_tags:
    cells = ct.find(tag)
    # if   tag == t_air : mu_ = 4 * np.pi * 1e-7  # Vacuum
    # elif tag == t_Itm : mu_ = 6.3e-3            # Copper
    # elif tag == t_pla : mu_ = 1.26e-6           # Iron
    mu_ = 4 * np.pi * 1e-7
    mu.x.array[cells]                  = np.full_like(cells, mu_, dtype=default_scalar_type)
    if   tag==t_Itp : J.x.array[cells] = np.full_like(cells,   I, dtype=default_scalar_type)
    elif tag==t_Itm : J.x.array[cells] = np.full_like(cells,  -I, dtype=default_scalar_type)
#========================> weak problem
# Coor=mesh.geometry.x
V = functionspace(mesh, ("Lagrange", 1))
# V = functionspace(mesh, ("Lagrange", 2))
# V = functionspace(mesh, ("DG", 1))
Coor=V.mesh.geometry.x

# degree=1
# curl_el = element("N1curl",   mesh.basix_cell(), degree, dtype=real_type)
# lagr_el = element("Lagrange", mesh.basix_cell(), degree, dtype=real_type)
# V = functionspace(mesh, mixed_element([curl_el, lagr_el]))
# V = functionspace(mesh, curl_el)

tdim = mesh.topology.dim
unk=tdim-1
facets = locate_entities_boundary(mesh, unk, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(V, unk, facets)
bc = dirichletbc(default_scalar_type(0), dofs, V)

u = TrialFunction(V)
v = TestFunction(V)
r,z=SpatialCoordinate(mesh)
# r = Function(V)
# r.x.array[:]=Coor[:,0]
# r.x.array[:]=1
a = (1 / mu) * dot(grad(u), grad(v)) * r * dx
L = J * v * r * dx
#========================> Solving
# sys = PETSc.Sys()  # type: ignore
# use_superlu = PETSc.IntType == np.int64
# if sys.hasExternalPackage("mumps") and not use_superlu:  # type: ignore
# 	mat_factor_backend = "mumps"
# elif sys.hasExternalPackage("superlu_dist"):  # type: ignore
# 	mat_factor_backend = "superlu_dist"
# else:
# 	if mesh_data.mesh.comm.size > 1:
# 		raise RuntimeError("This demo requires a parallel linear algebra backend.")
# 	else:
# 		mat_factor_backend = "petsc"

A_t = Function(V)
problem = LinearProblem(a, L, u=A_t, bcs=[bc],
	# 	petsc_options={
	# 	"ksp_type": "preonly",
	# 	"pc_type": "lu",
	# 	"pc_factor_mat_solver_type": mat_factor_backend,
	# },
)
# problem = LinearProblem(
# 	a,
# 	L,
# 	u=A_t,
# 	bcs=[],
# 	petsc_options={
# 		"ksp_type": "preonly",
# 		"pc_type": "lu",
# 		"pc_factor_mat_solver_type": mat_factor_backend,
# 	},
# )
problem.solve()
#========================> Magnetic field
W = functionspace(mesh, ("DG", 0, (mesh.geometry.dim, )))
B = Function(W)
# B_expr = Expression( as_vector(( -A_t.dx(1) , (A_t*r).dx(0)/r )) , W.element.interpolation_points())
B_expr = Expression( as_vector(( -A_t.dx(1) , A_t/r+A_t.dx(0) )) , W.element.interpolation_points())
B.interpolate(B_expr)

# A2= Function(W)
# A_expr=Expression(A_t,W.element.interpolation_points())
# A2.interpolate(A_expr)

############################################################################################################################################
print('\033[1;31m=====> Start Visualisation : {:.3f} s\033[0m'.format(time.time()-t0))
############################################################################################################################################
if CLEAR :
	os.system('rm '+d_plot+'/*')

if PLOT or DATA :
	#========================> Librairies
	import matplotlib.pyplot as plt
	import matplotlib        as mtp
	from matplotlib.patches import Circle
	cmap=mtp.colormaps['inferno']
	#========================> Data reshape
	top_imap = mesh.topology.index_map(mesh.topology.dim)
	num_cells = top_imap.size_local + top_imap.num_ghosts
	mesh.topology.create_connectivity(mesh.topology.dim, mesh.topology.dim)
	midpoints = compute_midpoints(mesh, mesh.topology.dim, np.arange(num_cells, dtype=np.int32))
	
	num_dofs = W.dofmap.index_map.size_local + W.dofmap.index_map.num_ghosts
	assert (num_cells == num_dofs)
	Bv = np.zeros((num_dofs, 3), dtype=np.float64)
	Bv[:, :mesh.geometry.dim] = B.x.array.real.reshape(num_dofs, W.dofmap.index_map_bs)
	NB=np.hypot(Bv[:,0],Bv[:,1])
	NA=A_t.x.array
	#========================> Interpolation
	Trin=mtp.tri.Triangulation( Coor[:,0],Coor[:,1] ) 
	Tric=mtp.tri.Triangulation( midpoints[:,0],midpoints[:,1] ) 
	Br  =mtp.tri.LinearTriInterpolator( Tric,Bv[:,0] )
	Bz  =mtp.tri.LinearTriInterpolator( Tric,Bv[:,1] )
	At  =mtp.tri.LinearTriInterpolator( Trin,NA      )
	#========================> Discretization
	Nv=int(1e3)
	Hr=np.linspace( 0,R,Nv)
	Hz=Nv*[0]
	Vr=Nv*[0]
	Vz=np.linspace(-R,R,Nv)
	#========================> plot lim
	xl0,xl1= 5.0e-3,2.5e-2
	yl0,yl1=-1.5e-2,1.5e-2

if PLOT :
	#========================> Visualisation

	#=====> Potential 
	figA,axA=plt.subplots(figsize=(10,10))
	axA.set_xlim(xl0,xl1)
	axA.set_ylim(yl0,yl1)
	axA.set_aspect('equal')
	fv=axA.tricontourf( Trin,NA,levels=int(1e3),cmap=cmap)
	figA.savefig(d_plot+'NA-I{:.2f}.png'.format(I)) ; plt.close(figA)

	Xt=np.mean( Tric.x[Tric.triangles],axis=1 )
	Tric.set_mask( abs(Xt)<=r0 )
	#=====> Field
	figF,axF=plt.subplots(figsize=(10,10))
	axF.set_xlim(xl0,xl1)
	axF.set_ylim(yl0,yl1)
	axF.set_aspect('equal')
	fv=axF.tricontourf( Tric,NB,levels=int(1e3),cmap=cmap)
	# axF.quiver( midpoints[:,0],midpoints[:,1],Bv[:,0],Bv[:,1],color='w' )
	
	#=====> Quiver
	figV,axV=plt.subplots(figsize=(10,10))
	axV.set_xlim(xl0,xl1)
	axV.set_ylim(yl0,yl1)
	axV.set_aspect('equal')
	axV.quiver( midpoints[:,0],midpoints[:,1],Bv[:,0],Bv[:,1] )
	figV.savefig(d_plot+'VB-I{:.2f}.png'.format(I)) ; plt.close(figV)

	#========================> Validation
	Br_v=Br(Vr,Vz)
	Bz_v=Bz(Vr,Vz)
	Br_h=Br(Hr,Hz)
	Bz_h=Bz(Hr,Hz)
	At_h=At(Hr,Hz)
	At_t=At(Hr,Nv*[ 5e-3])
	At_b=At(Hr,Nv*[-5e-3])
	
	# W_l=Circle((-R,0),rw,facecolor='none',edgecolor='r')
	W_r=Circle(( R,0),rw,facecolor='none',edgecolor='r')
	# axF.add_patch(W_l)
	axF.add_patch(W_r)
	axF.quiver(Vr,Vz,Br_v,Bz_v,color='w')
	axF.quiver(Hr,Hz,Br_h,Bz_h,color='w')
	figF.savefig(d_plot+'NB-I{:.2f}.png'.format(I)) ; plt.close(figF)
	
	figH,axH=plt.subplots()
	axH.plot(Hr,Bz_h,'k.',Hr,Br_h,'r.')
	figH.savefig(d_plot+'BH-I{:.2f}.png'.format(I)) ; plt.close(figH)
	
	figA,axA=plt.subplots()
	axA.plot(Hr,At_h,'k')
	axA.plot(Hr,At_t,'r')
	axA.plot(Hr,At_b,'b')
	figA.savefig(d_plot+'AH-I{:.2f}.png'.format(I)) ; plt.close(figA)
	
	# figV,axV=plt.subplots()
	# axV.plot(Vz,Bz_v,'k.',Vz,Br_v,'r.')
	# figV.savefig(d_plot+'BV-I{:.2f}.png'.format(I)) ; plt.close(figV)

#========================> Plate
if PLOT or DATA :
	Np=int(1e3)
	Pr=np.linspace( r0,Rp , Np )
	#=====> Magnetic field
	Br_pt=Br(Pr,Np*[Dz+    Ep]) ; Bz_pt=Bz(Pr,Np*[Dz+    Ep]) ; NB_pt=np.hypot(Br_pt,Bz_pt)
	Br_pm=Br(Pr,Np*[Dz+0.5*Ep]) ; Bz_pm=Bz(Pr,Np*[Dz+0.5*Ep]) ; NB_pm=np.hypot(Br_pm,Bz_pm)
	Br_pb=Br(Pr,Np*[Dz       ]) ; Bz_pb=Bz(Pr,Np*[Dz       ]) ; NB_pb=np.hypot(Br_pb,Bz_pb)
	#=====> Magnetic potential
	At_pt=At(Pr,Np*[Dz+    Ep]) 
	At_pm=At(Pr,Np*[Dz+0.5*Ep]) 
	At_pb=At(Pr,Np*[Dz       ]) 

if PLOT :
	figP,axP=plt.subplots()
	axP.plot( Pr,NB_pb,'b' , Pr,NB_pt,'r' )
	figP.savefig(d_plot+'Bplate-I{:.2f}.png'.format(I)) ; plt.close(figP)

if DATA :
	dat_plate=open(d_data+'BPlate-I{:.2f}.csv'.format(I),'w') 
	dat_plate.write('Radius,Bbot,Bmid,Btop,Abot,Amid,Atop\n')
	for n in range(Np) : dat_plate.write( '{:.12e} , {:.12e},{:.12e},{:.12e} , {:.12e},{:.12e},{:.12e}\n'.format( Pr[n] , NB_pb[n],NB_pm[n],NB_pt[n] , At_pb[n],At_pm[n],At_pt[n] ) )
	dat_plate.closed

############################################################################################################################################
print('\033[1;31m=====> Stop : {:.3f} s\033[0m'.format(time.time()-t0))
############################################################################################################################################