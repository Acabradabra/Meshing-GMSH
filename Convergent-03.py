#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import *
import sys
import pygmsh as pg
import gmsh   as gm

import matplotlib.pyplot as plt
import matplotlib

#===================================================================================
#                     Dimension
#===================================================================================

u=1e-3 ;
# Config='Long'
Config='Short'

# Size='Gros'
# Size='Mead0'
Size='Mead1'
# Size='Fine'

#================================> Geometry
if   Config=='Long' :
	Lin=30*u  ; #100*u ; // L pipe
	Lco=10*u  ; # L convergent
	Lij=10*u  ; # L inject
	Lwa=100*u ; # L wall
	# Ns=5      ; # N slices
	# r0=0.25   ;
elif Config=='Short' :
	Lin=30*u  ; #100*u ; // L pipe
	Lco=10*u  ; # L convergent
	Lij=10*u  ; # L inject
	Lwa=50*u  ; # L wall
	# Ns=0      ; # N slices
	# Ns=5      ; # N slices
	# Ns=10     ; # N slices
	# r0=0.20   ;

L0=10*u ; # Raf inlet
L1=10*u ; # Raf convergent
L2=10*u ; # Raf jet

Linj=Lin+Lco   ;
Lpipe=Linj+Lij ;
Ltot=Lpipe+Lwa ;

Din=10*u   ;
Dij=0.5*u  ;
Dou=10*u   ;
Dto=10*u   ;

# D0=3*u  ; # Raf jet
D0=1*u  ; # Raf jet

if Size=='Gros' :
	lin=1e-3; #1e-3 ; #5e-5 ; #1e-4 ;
	lpi=1e-3; #1e-3 ; #2e-4 ; #1e-4 ;
	lco=5e-4; #5e-4 ; #2e-5 ; #1e-4 ;
	lij=1e-4; #1e-4 ; #1e-5 ; #1e-4 ;
	lje=3e-4; #1e-3 ; #2e-5 ; #1e-4 ;
	lwa=1e-3; #1e-3 ; #2e-4 ; #1e-4 ;
	lou=6e-4; #1e-3 ; #5e-4 ; #1e-4 ;
	Nc =5   # N cells in corners
	Nxy=1*Nc
	Ns =5   # N slices
elif Size=='Mead0' :
	lin=2.5e-4 #5e-5
	lpi=1.0e-3 #2e-4
	lco=1.0e-4 #2e-5
	lij=5.0e-5 #1e-5
	lje=1.0e-4 #2e-5
	lwa=5.0e-4 #1e-4
	lou=5.0e-4 #1e-4
	Nc =1      # N cells in corners
	Nxy=5*Nc
	Ns=5       # N slices
elif Size=='Mead1' :
	lin=1.0e-4 #5e-5
	lpi=4.0e-4 #2e-4
	lco=4.0e-5 #2e-5
	lij=2.0e-5 #1e-5
	lje=4.0e-5 #2e-5
	lwa=2.0e-4 #1e-4
	lou=2.0e-4 #1e-4
	Nc =2      # N cells in corners
	Nxy=3*Nc
	Ns =10     # N slices
elif Size=='Fine' :
	lin=5e-5
	lpi=2e-4
	lco=2e-5
	lij=1e-5
	lje=2e-5
	lwa=1e-4 #2e-4
	lou=1e-4 #5e-4
	Nc  =5   # N cells in corners
	Nxy =1*Nc
	Ns  =10  # N slices

Lgeo=[Lin,Linj,Lpipe,Ltot,L0,L1,L2]
Dgeo=[Din,Dij,Dou,Dto,D0]
lgeo=[lin,lpi,lco,lij,lje,lwa,lou]

#================================> Boundary layer

Nt=5   ; # N added points for tightness
r=1.1  ; # Boundary layer ratio
Np0=17+Ns ;

RS0=1
r0=sqrt(3)/(4*r**(Ns-1)) ; # First ratio
print("r0 : ",r0) ;

h0_in=lin*r0 ; hN_in=h0_in*r**(Ns-1) ; ht_in=h0_in*(r**Ns-1)/(r-1) #; SN_in=lin*hN_in ; St_in=0.25*sqrt(3)*lin**2 ; RS_in=SN_in/St_in ; print('RS_in : ',RS_in)
h0_pi=lpi*r0 ; hN_pi=h0_pi*r**(Ns-1) ; ht_pi=h0_pi*(r**Ns-1)/(r-1) #; SN_pi=lpi*hN_pi ; St_pi=0.25*sqrt(3)*lpi**2 ; RS_pi=SN_pi/St_pi ; print('RS_pi : ',RS_pi)
h0_co=lco*r0 ; hN_co=h0_co*r**(Ns-1) ; ht_co=h0_co*(r**Ns-1)/(r-1) #; SN_co=lco*hN_co ; St_co=0.25*sqrt(3)*lco**2 ; RS_co=SN_co/St_co ; print('RS_co : ',RS_co)
h0_ij=lij*r0 ; hN_ij=h0_ij*r**(Ns-1) ; ht_ij=h0_ij*(r**Ns-1)/(r-1) #; SN_ij=lij*hN_ij ; St_ij=0.25*sqrt(3)*lij**2 ; RS_ij=SN_ij/St_ij ; print('RS_ij : ',RS_ij)
h0_je=lje*r0 ; hN_je=h0_je*r**(Ns-1) ; ht_je=h0_je*(r**Ns-1)/(r-1) #; SN_je=lje*hN_je ; St_je=0.25*sqrt(3)*lje**2 ; RS_je=SN_je/St_je ; print('RS_je : ',RS_je)
h0_wa=lwa*r0 ; hN_wa=h0_wa*r**(Ns-1) ; ht_wa=h0_wa*(r**Ns-1)/(r-1) #; SN_wa=lwa*hN_wa ; St_wa=0.25*sqrt(3)*lwa**2 ; RS_wa=SN_wa/St_wa ; print('RS_wa : ',RS_wa)
h0_ou=lou*r0 ; hN_ou=h0_ou*r**(Ns-1) ; ht_ou=h0_ou*(r**Ns-1)/(r-1) #; SN_ou=lou*hN_ou ; St_ou=0.25*sqrt(3)*lou**2 ; RS_ou=SN_ou/St_ou ; print('RS_ou : ',RS_ou)
H0=[h0_in,h0_pi,h0_co,h0_ij,h0_je,h0_wa,h0_ou]

RS=sqrt(3)/(4*r0*r**(Ns-1)) ;
print("RS : ",RS) ;

#================================> Normals
Xab=Lco    
Yab=Dij-Din
NAB=sqrt( Xab**2+Yab**2 )
nab=array([Xab,Yab])/NAB

Nx= Yab/NAB ;
Ny=-Xab/NAB ;
Nm=array([Nx,Ny])
N0=array([0 ,-1])
NA=0.5*(Nm+N0) ; NA/=-NA[1]
NB=0.5*(Nm+N0) ; NB/=-NB[1] #; print(NB)
Vn=[NA,Nm,NB,nab]

#================================> N point progression
lm0 =0.5*(lco+lij)
Nco = round( (NAB/lm0)+1           ,0 )
rco = round( (lco/lij)**(1/(Nco-1)),14)
Nco2= round(log(1+(rco-1)*NAB/lij)/log(rco),0)
lm1=0.5*(lij+lou)
Lwa=Dou-Dij
Nwa = round( (Lwa/lm1)+1           ,0 )
rwa = round( (lou/lij)**(1/(Nwa-1)),14)
Nwa2= round(log(1+(rwa-1)*Lwa/lij)/log(rwa),0)

Nij = round(  Lij/lij+1            ,0 )
lco2=lij*rco**(Nco-1)
NAB2=lij*( (rco**Nco -1)/(rco-1) )
NAB3=lij*( (rco**Nco2-1)/(rco-1) )

print('Nco0 : {:.0f}  ,  Nco2 : {:.0f}  ,  rco : {:.14f}'.format(Nco,Nco2,rco))
print('Nwa0 : {:.0f}  ,  Nwa2 : {:.0f}  ,  rwa : {:.14f}'.format(Nwa,Nwa2,rwa))
print('Nij  : {:.0f}'.format(Nij))
print('Er : {:.5e}'.format(abs(lco-lco2)/lco))
print('Er : {:.5e}'.format(abs(NAB-NAB2)/NAB))
print('Er : {:.5e}'.format(abs(NAB-NAB3)/NAB))

#================================> Input
# lm_in=(lin+lpi+lco)/3
lm_in=lco
Nin=round(Lin/lm_in,0)+1
print('Nin : {}'.format(Nin))

if '--Norm' in sys.argv :
	P0=array([0    ,Din])
	PA=array([Lin  ,Din])
	PB=array([Linj ,Dij])
	P1=array([Lpipe,Dij])
	AB=PB-PA
	print(AB)
	print(Xab,Yab)
	print(Xab*Nx+Yab*Ny)
	print(sqrt(Nx**2+Ny**2))
	
	Pm=0.5*(PA+PB)
	Pnm=Pm+Nm*1e-3
	PnA=PA+NA*1e-3
	PnB=PB+NB*1e-3
	
	fig,ax=plt.subplots()
	ax.plot( [0,Lin,Linj,Lpipe],[Din,Din,Dij,Dij],'-ok' )
	ax.plot( [Pm[0],Pnm[0]],[Pm[1],Pnm[1]],'-r' )
	ax.plot( [PA[0],PnA[0]],[PA[1],PnA[1]],'-r' )
	ax.plot( [PB[0],PnB[0]],[PB[1],PnB[1]],'-r' )
	ax.set_aspect(1)
	ax.set_xlim(0.025,0.045)
	plt.show()
	
	sys.exit()

#===================================================================================
#                     Geometry
#===================================================================================
def ht(h0,r,n) : return(h0*(r**n-1)/(r-1))
#===================================================================================
def Corner(Pm,N0,N1,n0,n1,nm,l,h,Nc) :
	Vp=[]
	P0=array(Pm)+Nc*l*array(N0)
	P1=array(Pm)+Nc*l*array(N1)
	Vp.append(geom.add_point( P0+h*array(n0),l ))
	Vp.append(geom.add_point( Pm+h*array(nm),l ))
	Vp.append(geom.add_point( P1+h*array(n1),l ))
	return( Vp )
#===================================================================================
def Bezier(P0,P1,Vp,geo,Nt) :
	Np=len(Vp)
	Vp2=[ P0 ]
	for p in range(Np-1) :
		Xl=linspace( Vp[p][0][0],Vp[p+1][0][0],Nt+2 )
		Yl=linspace( Vp[p][0][1],Vp[p+1][0][1],Nt+2 )
		ll=linspace( Vp[p][1]   ,Vp[p+1][1]   ,Nt+2 )
		# if p<Np-1 : end=1
		if p<Np-2 : end=1
		else      : end=0
		for n in range(Nt+end) :
			Vp2.append( geom.add_point( [Xl[n+1],Yl[n+1]] , ll[n+1]   ) )
	Vp2.append(P1)
	return( Vp2 )
#===================================================================================
def Top(P0,P1,Lgeo,Dgeo,lgeo,H0,Vn,r,n) :
	[Lin,Linj,Lpipe,Ltot,L0,L1,L2]=Lgeo
	[Din,Dij,Dou,Dto,D0]          =Dgeo
	[lin,lpi,lco,lij,lje,lwa,lou] =lgeo
	[h0_in,h0_pi,h0_co,h0_ij,h0_je,h0_wa,h0_ou]=H0
	[NA,Nm,NB,nab]=Vn
	if n>=0 :
		ht_in=ht(h0_in,r,n)
		ht_pi=ht(h0_pi,r,n)
		ht_co=ht(h0_co,r,n)
		ht_ij=ht(h0_ij,r,n)
		ht_je=ht(h0_je,r,n)
		ht_wa=ht(h0_wa,r,n)
		ht_ou=ht(h0_ou,r,n)
	else :
		ht_in=0
		ht_pi=0
		ht_co=0
		ht_ij=0
		ht_je=0
		ht_wa=0
		ht_ou=0
	Cwallx=[
		[[Lpipe+ht_ou,Dou-Nxy*lou],lou],
		[[Lpipe+ht_je,Dij+D0     ],lje],
		[[Lpipe+ht_ij,Dij+ Nc*lij],lij]]
	Cpipe=[
		[[Lin-Nc*lco,Din-ht_co],lco],
		[[Lin-L1    ,Din-ht_pi],lpi],
		[[L0        ,Din-ht_pi],lpi],
		[[0         ,Din-ht_in],lin]]

	#=====================> Corner Points
	Pcorxy=Corner( [Lpipe,Dou] , [1,0],[ 0,-1] , [0,-1],[1, 0],[1,-1] , lou,ht_ou,Nxy )
	Pcopi0=Corner( [Lpipe,Dij] , [0,1],[-1, 0] , [1, 0],[0,-1],[1,-1] , lij,ht_ij,Nc )
	Pcopi1=Corner( [Linj ,Dij] , [1,0], -nab   , [0,-1], Nm   ,NB     , lij,ht_ij,Nc )
	Pcopi2=Corner( [Lin  ,Din] ,  nab ,[-1, 0] ,  Nm   ,[0,-1],NA     , lco,ht_co,Nc )
	Pin=geom.add_point( Cpipe[-1][0],Cpipe[-1][1] )

	#=====================> Line Points
	Pwall =[ P0       ,geom.add_point([Ltot,Dou-ht_wa],lwa)]
	Pwally=[       Pwall[-1],Pcorxy[0]]
	Pwallx=Bezier(Pcorxy[-1],Pcopi0[0],Cwallx   ,geom,Nt)
	Ppipe0=[      Pcopi0[-1],Pcopi1[0]]
	Ppipe1=[      Pcopi1[-1],Pcopi2[0]]
	Ppipe2=Bezier(Pcopi2[-1],Pin      ,Cpipe    ,geom,Nt)

	Points=[Pcorxy,Pcopi0,Pcopi1,Pcopi2,Pwall ,Pwally,Pwallx,Ppipe0,Ppipe1,Ppipe2]

	#=====================> Lines
	l_wa = geom.add_line( Pwall[0], Pwall[1])
	l_wy = geom.add_line(Pwally[0],Pwally[1]) # Wall y
	l_cy = geom.add_line(Pcorxy[0],Pcorxy[1]) # corner y
	l_cx = geom.add_line(Pcorxy[1],Pcorxy[2]) # corner x
	l_wx = geom.add_bspline(Pwallx) # wall x
	l_c0 = geom.add_bspline(Pcopi0) # Corner pipe 0
	l_c1 = geom.add_bspline(Pcopi1) # Corner pipe 1
	l_c2 = geom.add_bspline(Pcopi2) # Corner pipe 2
	l_p0 = geom.add_line(Ppipe0[0],Ppipe0[1]) # Line pipe 0
	l_p1 = geom.add_line(Ppipe1[0],Ppipe1[1]) # Line pipe 1
	l_p2 = geom.add_bspline(Ppipe2) # Line pipe 2
	l_in = geom.add_line(Ppipe2[-1],P1) # inlet

	Lines=[l_wa,l_wy,l_cy,l_cx,l_wx,l_c0,l_c1,l_c2,l_p0,l_p1,l_p2,l_in]

	return(Points,Lines)
#===================================================================================

with pg.geo.Geometry() as geom:
	
	#===================================> Internal surface
	Paxis0=[]
	Paxis0.append(geom.add_point([0        ,0     ],lin)) # Axis
	Paxis0.append(geom.add_point([L0       ,0     ],lpi)) # Axis
	Paxis0.append(geom.add_point([Lin-L1   ,0     ],lpi)) # Axis
	Paxis0.append(geom.add_point([Lin      ,0     ],lco)) # Axis
	Paxis0.append(geom.add_point([Linj     ,0     ],lij)) # Axis
	Paxis0.append(geom.add_point([Lpipe    ,0     ],lij)) # Axis
	Paxis0.append(geom.add_point([Lpipe+L2 ,0     ],lje)) # Axis
	Paxis0.append(geom.add_point([Ltot     ,0     ],lwa)) # Axis,outlet
	l_ax0 = geom.add_bspline(Paxis0) # Axis
	(Point0,Line0)=Top(Paxis0[-1],Paxis0[0],Lgeo,Dgeo,lgeo,H0,Vn,r,Ns)
	[l_wa0,l_wy0,l_cy0,l_cx0,l_wx0,l_c00,l_c10,l_c20,l_p00,l_p10,l_p20,l_in0]=Line0
	[Pcorxy0,Pcopi00,Pcopi10,Pcopi20,Pwall0 ,Pwally0,Pwallx0,Ppipe00,Ppipe10,Ppipe20]=Point0

	# geom.set_transfinite_curve(l_c00,2*Nc+1,'Progression',1)
	# geom.set_transfinite_curve(l_c10,2*Nc+1,'Progression',1)
	# geom.set_transfinite_curve(l_c20,2*Nc+1,'Progression',1)

	lli=geom.add_curve_loop([l_ax0,l_wa0,l_wy0,l_cy0,l_cx0,l_wx0,l_c00,l_p00,l_c10,l_p10,l_c20,l_p20,l_in0])
	psi=geom.add_plane_surface(lli)
	
	#===================================> Boundary layer
	# Boundary_Surfaces=[]
	for n in range(Ns-1,-1,-1) :
		(Point1,Line1)=Top(Pwall0[-1],Ppipe20[-1],Lgeo,Dgeo,lgeo,H0,Vn,r,n)
		
		[l_wa1,l_wy1,l_cy1,l_cx1,l_wx1,l_c01,l_c11,l_c21,l_p01,l_p11,l_p21,l_in1]=Line1
		[Pcorxy1,Pcopi01,Pcopi11,Pcopi21,Pwall1 ,Pwally1,Pwallx1,Ppipe01,Ppipe11,Ppipe21]=Point1
		
		lls=geom.add_curve_loop([l_wa1,l_wy1,l_cy1,l_cx1,l_wx1,l_c01,l_p01,l_c11,l_p11,l_c21,l_p21,l_in1,-l_p20,-l_c20,-l_p10,-l_c10,-l_p00,-l_c00,-l_wx0,-l_cx0,-l_cy0,-l_wy0])
		pss=geom.add_plane_surface(lls)
		# Boundary_Surfaces.append(pss)
		
		[l_wa0,l_wy0,l_cy0,l_cx0,l_wx0,l_c00,l_c10,l_c20,l_p00,l_p10,l_p20,l_in0]=[l_wa1,l_wy1,l_cy1,l_cx1,l_wx1,l_c01,l_c11,l_c21,l_p01,l_p11,l_p21,l_in1]
		[Pcorxy0,Pcopi00,Pcopi10,Pcopi20,Pwall0 ,Pwally0,Pwallx0,Ppipe00,Ppipe10,Ppipe20]=[Pcorxy1,Pcopi01,Pcopi11,Pcopi21,Pwall1 ,Pwally1,Pwallx1,Ppipe01,Ppipe11,Ppipe21]

	# geom.set_recombined_surfaces(Boundary_Surfaces)

	#===================================> Output
	geom.env.removeAllDuplicates()

	geom.save_geometry('Geo.geo_unrolled')