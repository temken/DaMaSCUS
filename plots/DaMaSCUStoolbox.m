(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["DaMaSCUStoolbox`"]

ReducedMass::usage="ReducedMass[m1,m2] calculates m1 m2/(m1+m2).";
InUnits::usage="InUnits[n,unit] converts a number n in natural units into given unit.";
SphericalCoordinates::usage="SphericalCoordinates[r,\[Theta],\[Phi]] returns the vector with spherical coordinates (r,\[Theta],\[Phi]).";
Spherical2Cartesian::usage="Spherical2Cartesian[r,\[Theta],\[Phi]] returns the vector with spherical coordinates (r,\[Theta],\[Phi]).";
Cartesian2Spherical::usage="Cartesian2Spherical[{x,y,z}] returns the spherical coordinates (r,\[Theta],\[Phi]) corresponding to the vector (x,y,z).";
Cartesian2Spherical2::usage="Cartesian2Spherical[{x,y,z}] returns the spherical coordinates (r,cos(\[Theta]),\[Phi]) corresponding to the vector (x,y,z).";
Equat2Gal::usage="Equat2Gal[vector,T] converts a vector given in the equatorial frame of the epoch T to the galactic frame.";
Gal2Equat::usage="Gal2Equat[vector,T] converts a vector given in the galactic frame to the equatorial frame of the epoch T.";
Lab2Equat::usage="Lab2Equat[vector,{lat,lon},nJ2000] converts a vector given in the local laboratory frame to the equatorial frame. The laborator's position is given in latitude and longitude, its time is given in fractional days since J2000.0.";
Equat2Lab::usage="Equat2Lab[vector,{lat,lon},nJ2000] converts a vector given in the equatorial frame at nJ2000 to the local laboratory frame. The laborator's position is given in latitude and longitude, the time nJ2000 is given in fractional days since J2000.0.";
Lab2Gal::usage="Lab2Gal[vector,{lat,lon},nJ2000] converts a vector given in the local laboratory frame to the galactic rest frame. The laborator's position is given in latitude and longitude, its time is given in fractional days since J2000.0.";
Gal2Lab::usage="Gal2Lab[vector,{lat,lon},nJ2000] converts a vector given in the galactic frame at nJ2000 to the local laboratory frame. The laborator's position is given in latitude and longitude, the time nJ2000 is given in fractional days since J2000.0.";
FractionalDays::usage="FractionalDays[{d,m,y},{h,min,sec}] converts a given date and time into the number of fractional days since J2000.0 .";
GMST::usage="GMST[nJ2000] calculates the Greenwich Mean Sidereal Time for a given number of fractional days since J2000.0 .";
GAST::usage="GAST[nJ2000] calculates the Greenwich Apparent Sidereal Time for a given number of fractional days since J2000.0 .";
LAST::usage="LAST[nJ2000,\[Lambda]] calculates the Local Apparent Sidereal Time for a given number of fractional days since J2000.0 at the longitude \[Lambda].";
SecondsToClock::usage="SecondsToClock[seconds] converts a number between 0 and 86400 into a clock format.";
EarthVelocity::usage="EarthVelocity[nJ2000] calculates the earth velocity vector in the galactic rest frame for a given fractional day since J2000.0 .";
RotationVelocity::usage="RotationVelocity[n,{latitude,longitude}] calculates the earth rotation component of a laboratory at (lat,lon) for a given number of fractional days since J2000.0 .";
LabVelocity::usage="LabVelocity[nJ2000,{latitude,longitude}] returns the laboratory's velocity in the galactic rest frame, for a time given in fractiona days since J2000.0. The laboratory is located at (lat,lon).";
DetectorPosition::usage="DetectorPosition[nJ2000,{lat,lon},depth] returns the vector (in galactic rest frame) pointing towards a detector at (lat, lon) at a time given by the fractional days since J2000.0. ";
DetectorAngle::usage="DetectorAngle[last,date,{lat,lon}] returns the \[Theta] angle of a detector."
DetectorAngle2::"DetectorAngle2[nJ2000.0,{lat,lon}]";
EtaFunction::usage="EtaFunction[vmin,vEarth] returns the analytic \[Eta] function.";
vMinimum::usage="vMinimum[ER_,m\[Chi],A] calculates the minimal velocity of a DM particle of mass m\[Chi] to cause a recoil energy ER on a nucleus A"
VelocityDistribution::usage="VelocityDistribution[v,vE], where v is a velocity vector and vE is the Earth's velocity vector. Standard halo model velocity distribution function."
VelocityDistribution2::usage="VelocityDistribution2[v,cos\[Theta],\[Phi],vE]. Does the same as VelocityDistribution[], but takes spherical coordinates as argument."
SpeedDistribution::usage="SpeedDistribution[v,vE] returns the marginal speed distribution by integrating f over the direction angles."
AverageDMSpeed::usage="AverageDMSpeed[{dd,mm,yyyy},{h,m,s}] returns the average DM velocity."


(* ::Section::Initialization:: *)
(*(*(*(*(*Global Constants*)*)*)*)*)


(* ::Subsection::Initialization:: *)
(*(*(*(*(*Physical Parameters and Units*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*Units*)*)*)*)*)


(* ::Input::Initialization:: *)
GeV=1;
eV=10^-9 GeV;
keV = 10^-6 GeV;
MeV=10^-3 GeV;
TeV=10^3 GeV;
gram=5.617977528089887 10^23 GeV;
kg=10^3 gram;
cm=5.067 10^13 GeV^-1;
meter=100 cm;
km=10^3 meter;
fm=10^-15 meter;
pb=10^-36 cm^2;
sec=299792458 meter;
hr=3600 sec;
day=24hr;
yr=365.24 day;
erg=gram (cm/sec)^2;
Newton=kg meter/sec^2;
Joule=Newton meter;
Watt=Joule/sec;
Farad=Newton/Coulomb;
Tesla=(Newton sec)/(Coulomb meter);
Kelvin=8.62 10^-14 GeV;
pc=3.08567758 10^16 meter;
kpc=10^3 pc;
Mpc=10^6 pc;


(* ::Text::Initialization:: *)
(*(*(*(*(*Specific Parameters*)*)*)*)*)


(* ::Input::Initialization:: *)
mPl=1.2209 10^19 GeV;
\[Alpha]EM = 1/137.035999139;
GNewton=1/mPl^2;
GFermi=1.16637 10^-5 GeV^-2;
mProton=0.938GeV;
mElectron=0.511MeV;
mNucleon=0.932GeV;
(*ElementaryCharge=1.602 10^-19Coulomb;*)
ElementaryCharge=0.30282212;


(* ::Text::Initialization:: *)
(*(*(*(*(*Dark Matter Halo Parameters*)*)*)*)*)


(* ::Input::Initialization:: *)
v0=220 km/sec;
vesc=544 km/sec;
\[Rho]\[Chi]=0.3 GeV/cm^3;
Nesc=\[Pi] v0^2 (Sqrt[\[Pi]]v0 Erf[vesc/v0]-2vesc Exp[-vesc^2/v0^2]);


(* ::Text::Initialization:: *)
(*(*(*(*(*Astrophysical Parameters*)*)*)*)*)


(* ::Input::Initialization:: *)
m\[Earth]=5.972 10^24 kg;
r\[Earth]=6371km;
\[Rho]\[Earth]=5.51 gram/cm^3;
m\[Sun]=1.989 10^30  kg;
r\[Sun]=6.957 10^8  meter;
\[Rho]\[Sun]=151 gram/cm^3;


(* ::Section::Initialization:: *)
(*(*(*(*(*Functions*)*)*)*)*)


(* ::Input::Initialization:: *)
Begin["`Private`"]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*General physical functions*)*)*)*)*)


(* ::Input::Initialization:: *)
ReducedMass[m1_,m2_]:=(m1 m2)/(m1+m2)
InUnits[number_,units_]:=number/units


(* ::Subsection::Initialization:: *)
(*(*(*(*(*Coordinate Systems*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*Spherical Coordinates*)*)*)*)*)


(* ::Input::Initialization:: *)
SphericalCoordinates[r_,\[Theta]_,\[Phi]_]:={r Sin[\[Theta]]Cos[\[Phi]],r Sin[\[Theta]] Sin[\[Phi]],r Cos[\[Theta]]}
Spherical2Cartesian[r_,\[Theta]_,\[Phi]_]:={r Sin[\[Theta]]Cos[\[Phi]],r Sin[\[Theta]] Sin[\[Phi]],r Cos[\[Theta]]}
Cartesian2Spherical[{x_,y_,z_}]:={Sqrt[x^2+y^2+z^2],ArcCos[z/Sqrt[x^2+y^2+z^2]],Mod[ArcTan[x,y],2\[Pi]]}
Cartesian2Spherical2[{x_,y_,z_}]:={Sqrt[x^2+y^2+z^2],z/Sqrt[x^2+y^2+z^2],Mod[ArcTan[x,y],2\[Pi]]}


(* ::Text::Initialization:: *)
(*(*(*(*(*Going from equatorial rectangular coordinates to galactic coordinates and back (arxiv: 1312:1355v2)*)*)*)*)*)


(* ::Input::Initialization:: *)
Equat2Gal[vec_,T_:0]:=Module[{arcsecond,\[Zeta]A,zA,\[Theta]A,P,lCP,\[Alpha]GP,\[Delta]GP,M},
arcsecond=1/60^2 \[Degree];
\[Zeta]A:=2306.083227 arcsecond T + 0.298850 arcsecond T^2;
zA:=2306.077181 arcsecond T + 1.092735 arcsecond T^2;
\[Theta]A:=2004.191903 arcsecond T -0.429493 arcsecond T^2;
P:={{Cos[\[Zeta]A]Cos[\[Theta]A]Cos[zA]-Sin[\[Zeta]A]Sin[zA],-Sin[\[Zeta]A]Cos[\[Theta]A]Cos[zA]-Cos[\[Zeta]A]Sin[zA],-Sin[\[Theta]A]Cos[zA]},{Cos[\[Zeta]A]Cos[\[Theta]A]Sin[zA]+Sin[\[Zeta]A]Cos[zA],-Sin[\[Zeta]A]Cos[\[Theta]A]Sin[zA]+Cos[\[Zeta]A]Cos[zA],-Sin[\[Theta]A]Sin[zA]},{Cos[\[Zeta]A]Sin[\[Theta]A],-Sin[\[Zeta]A]Sin[\[Theta]A],Cos[\[Theta]A]}};
lCP=122.932 \[Degree];
\[Alpha]GP=192.85948 \[Degree];
\[Delta]GP=27.12825 \[Degree];
M={{-Sin[lCP] Sin[\[Alpha]GP]-Cos[lCP] Cos[\[Alpha]GP] Sin[\[Delta]GP],Sin[lCP] Cos[\[Alpha]GP]-Cos[lCP] Sin[\[Alpha]GP] Sin[\[Delta]GP],Cos[lCP] Cos[\[Delta]GP]},{Cos[lCP] Sin[\[Alpha]GP]-Sin[lCP] Cos[\[Alpha]GP] Sin[\[Delta]GP],-Cos[lCP] Cos[\[Alpha]GP]-Sin[lCP] Sin[\[Alpha]GP] Sin[\[Delta]GP],Sin[lCP] Cos[\[Delta]GP]},{Cos[\[Alpha]GP] Cos[\[Delta]GP],Sin[\[Alpha]GP] Cos[\[Delta]GP], Sin[\[Delta]GP]}};
Return[M.Inverse[P].vec]
]
Gal2Equat[vec_,T_:0]:=Module[{arcsecond,\[Zeta]A,zA,\[Theta]A,P,lCP,\[Alpha]GP,\[Delta]GP,M},
arcsecond=1/60^2 \[Degree];
\[Zeta]A:=2306.083227 arcsecond T + 0.298850 arcsecond T^2;
zA:=2306.077181 arcsecond T + 1.092735 arcsecond T^2;
\[Theta]A:=2004.191903 arcsecond T -0.429493 arcsecond T^2;
P:={{Cos[\[Zeta]A]Cos[\[Theta]A]Cos[zA]-Sin[\[Zeta]A]Sin[zA],-Sin[\[Zeta]A]Cos[\[Theta]A]Cos[zA]-Cos[\[Zeta]A]Sin[zA],-Sin[\[Theta]A]Cos[zA]},{Cos[\[Zeta]A]Cos[\[Theta]A]Sin[zA]+Sin[\[Zeta]A]Cos[zA],-Sin[\[Zeta]A]Cos[\[Theta]A]Sin[zA]+Cos[\[Zeta]A]Cos[zA],-Sin[\[Theta]A]Sin[zA]},{Cos[\[Zeta]A]Sin[\[Theta]A],-Sin[\[Zeta]A]Sin[\[Theta]A],Cos[\[Theta]A]}};
lCP=122.932 \[Degree];
\[Alpha]GP=192.85948 \[Degree];
\[Delta]GP=27.12825 \[Degree];
M={{-Sin[lCP] Sin[\[Alpha]GP]-Cos[lCP] Cos[\[Alpha]GP] Sin[\[Delta]GP],Sin[lCP] Cos[\[Alpha]GP]-Cos[lCP] Sin[\[Alpha]GP] Sin[\[Delta]GP],Cos[lCP] Cos[\[Delta]GP]},{Cos[lCP] Sin[\[Alpha]GP]-Sin[lCP] Cos[\[Alpha]GP] Sin[\[Delta]GP],-Cos[lCP] Cos[\[Alpha]GP]-Sin[lCP] Sin[\[Alpha]GP] Sin[\[Delta]GP],Sin[lCP] Cos[\[Delta]GP]},{Cos[\[Alpha]GP] Cos[\[Delta]GP],Sin[\[Alpha]GP] Cos[\[Delta]GP], Sin[\[Delta]GP]}};
Return[P.Inverse[M].vec]
]


(* ::Text::Initialization:: *)
(*(*(*(*(*Transformation between Lab and Equat Frame*)*)*)*)*)


(* ::Input::Initialization:: *)
Lab2Equat[vec_,{lat_,lon_},nJ2000_]:=Module[{\[Theta],\[Phi],N},
\[Theta]=\[Pi]/2-lat;
\[Phi]=(2\[Pi])/86400 LAST[nJ2000,lon];
N={{-Sin[\[Phi]],-Cos[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Cos[\[Phi]]},{Cos[\[Phi]],-Cos[\[Theta]]Sin[\[Phi]],Sin[\[Theta]]Sin[\[Phi]]},{0,Sin[\[Theta]],Cos[\[Theta]]}};
Return[N.vec]
]
Equat2Lab[vec_,{lat_,lon_},nJ2000_]:=Module[{\[Theta],\[Phi],N},
\[Theta]=\[Pi]/2-lat;
\[Phi]=(2\[Pi])/86400 LAST[nJ2000,lon];
N={{-Sin[\[Phi]],-Cos[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Cos[\[Phi]]},{Cos[\[Phi]],-Cos[\[Theta]]Sin[\[Phi]],Sin[\[Theta]]Sin[\[Phi]]},{0,Sin[\[Theta]],Cos[\[Theta]]}};
Return[Inverse[N].vec]
]


(* ::Text::Initialization:: *)
(*(*(*(*(*Transformation between lab and galactic frame*)*)*)*)*)


(* ::Input::Initialization:: *)
Lab2Gal[vec_,{lat_,lon_},nJ2000_]:=Equat2Gal[Lab2Equat[vec,{lat,lon},nJ2000],nJ2000/36525]
Gal2Lab[vec_,{lat_,lon_},nJ2000_]:=Equat2Lab[Gal2Equat[vec,nJ2000/36525],{lat,lon},nJ2000]



(* ::Subsection::Initialization:: *)
(*(*(*(*(*Sidereal Time*)*)*)*)*)


(* ::Input::Initialization:: *)
FractionalDays[{d_,m_,y_},{h_:0,min_:0,sec_:0}]:=Module[{n=0},
If[m==1||m==2,
n+=Floor[365.25 (y-1)]+Floor[30.61((m+12)+1)];
,
n+=Floor[365.25 (y)]+Floor[30.61(m+1)];
];
n+=d-730563.5+h/24+min/(24 60)+sec/(24 60^2);
Return[n]
]
GMST[n_]:=Module[{T},
T=n/36525;
Return[Mod[86400(0.7790572732640+0.00273781191135448n +Mod[n,1])+0.00096707+307.47710227T+0.092772113T^2,86400]];
]
EqEq[n_]:=Module[{T,\[CapitalOmega],L,\[CapitalDelta]\[Psi],\[Epsilon]A},
T=n/36525;
\[CapitalOmega]=125.04455501\[Degree]-n 0.05295376\[Degree];
L=280.47\[Degree]+n 0.98565\[Degree];
\[CapitalDelta]\[Psi]=-1.1484 Sin[\[CapitalOmega]]-0.0864Sin[2L];
\[Epsilon]A=23.4392794444\[Degree]-0.01301021361T^2;
Return[\[CapitalDelta]\[Psi] Cos[\[Epsilon]A]+0.000176Sin[\[CapitalOmega]]+0.000004Sin[2\[CapitalOmega]]];
]
GAST[n_]:=GMST[n]+EqEq[n]
LAST[n_,\[Lambda]_]:=Mod[GAST[n]+\[Lambda]/(2\[Pi]) 86400,86400];
SecondsToClock[sec_]:=Module[{h,m,s},
h=Floor[sec/3600];
m=Floor[(sec-3600h)/60];
s=sec-3600h-60m;
Return[{h,m,s}]
]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*Earth and laboratory velocity*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*Earth Velocity and fractional days from arxiv:1312:1355v2*)*)*)*)*)


(* ::Input::Initialization:: *)
EarthVelocity[n_:0]:=Module[{e,L,\[Omega],T,ex,ey,ve,uE},
e=0.01671; (*Ellipticity of earth's orbit*)
L=Mod[280.460 \[Degree] +n 0.9856474 \[Degree],2\[Pi]];
\[Omega]=Mod[282.932 \[Degree] +n 0.0000471 \[Degree] ,2\[Pi]];
T=n/36525;
ex={0.054876,-0.494109,0.867666}+{-0.024232,-0.002689,1.546 10^-6}T;
ey={0.993824,0.110992,0.000352}+{0.001316,-0.011851,0.021267}T;
ve=29.79 km/sec;
uE=-ve (Sin[L]+e Sin[2L-\[Omega]])ex+ve (Cos[L]+e Cos[2L-\[Omega]])ey;
Return[{0,220 ,0} km/sec+{11.1 ,12.2 ,7.3 } km/sec+uE]
]
RotationVelocity[n_:0,{latitude_,longitude_}]:=Module[{ex,ey,T,tLAST,vEq,\[Omega]Rot,GMST,\[Delta],vRot},
T=n/36525;
ex={1,0,0};
ey={0,1,0};
vEq=0.465 km/sec;
\[Omega]Rot = (2\[Pi])/86400;
tLAST=LAST[n,longitude];
GMST=Mod[18.697374558+24.06570982441908 n,24];
\[Delta]=\[Omega]Rot tLAST;
vRot=-vEq Cos[latitude](Sin[\[Delta]]Equat2Gal[ex,T]-Cos[\[Delta]]Equat2Gal[ey,T]);
Return[vRot]
]
LabVelocity[n_:0,{latitude_,longitude_}]:=EarthVelocity[n]+RotationVelocity[n,{latitude,longitude}]


(* ::Subsection::Initialization:: *)
(*(*(*(*(*Dark Matter Functions*)*)*)*)*)


(* ::Input::Initialization:: *)
VelocityDistribution[v_,vE_]:=1/Nesc Exp[-((v-vE).(v-vE)/v0^2)]
VelocityDistribution2[v_?NumericQ,cos\[Theta]_?NumericQ,\[Phi]_?NumericQ,vE_]:=Module[{ve},
ve=Cartesian2Spherical[vE];
Return[(1/(\[Pi] v0^2))^(3/2) Exp[-(v^2+ve[[1]]^2)/v0^2]Exp[-((2v ve[[1]])/v0^2)(cos\[Theta] Cos[ve[[2]]]+Cos[\[Phi]-ve[[3]]] Sqrt[1-cos\[Theta]^2] Sin[ve[[2]]])]]
]
SpeedDistribution[v_?NumericQ,vE_]:=Module[{cos\[Theta],\[Phi]},Return[NIntegrate[v^2 VelocityDistribution2[v,cos\[Theta],\[Phi],vE],{cos\[Theta],-1,1},{\[Phi],0,2\[Pi]}]];
]


(* ::Input::Initialization:: *)
AverageDMSpeed[date_,time_]:=Module[{angles,\[Theta],\[Phi],v,Cos\[Delta],vO,vEsc,int,int2,ve,vE},
vE=EarthVelocity[FractionalDays[date,time]];
angles=Cartesian2Spherical[vE];
int=Simplify[Integrate[Sin[\[Theta]]v^3 1/Nesc Exp[-((v^2+ve^2+2v ve Cos\[Delta])/vO^2)],{v,0,-Cos\[Delta] ve+Sqrt[vEsc^2-ve^2 (1-Cos\[Delta]^2)]}]];
int2=int/.{Cos\[Delta]->Cos[angles[[2]]]Cos[\[Theta]]+Cos[\[Phi]-angles[[3]]]Sin[angles[[2]]]Sin[\[Theta]]};
Return[NIntegrate[int2/.{ve->angles[[1]],vEsc->vesc,vO->v0},{\[Theta],0,\[Pi]},{\[Phi],0,2\[Pi]}]];
]


(* ::Input::Initialization:: *)
EtaFunction[vmin_?NumericQ,{v1_?NumericQ,v2_?NumericQ,v3_?NumericQ}]:=Module[{x,y,z,NEsc},
NEsc=Erf[vesc/v0]-2/Sqrt[\[Pi]] vesc/v0 Exp[-vesc^2/v0^2];
x=vmin/v0;
y=Norm[{v1,v2,v3}]/v0;
z=vesc/v0;
Which[
y+z<x,
Return[0]
,
Norm[y-z]<x<y+z,
Return[1/(2NEsc v0 y) (Erf[z]-Erf[x-y]-2/Sqrt[\[Pi]] (y+z-x)Exp[-z^2])]
,
z>y&&x<Norm[y-z],
Return[1/(2NEsc v0 y) (Erf[x+y]-Erf[x-y]-4/Sqrt[\[Pi]] y Exp[-z^2])]
,
z<y&&x<Norm[y-z],
Return[1/(v0 y)]
]
]


(* ::Input::Initialization:: *)
vMinimum[ER_?NumericQ,m\[Chi]_,A_:131]:=Sqrt[mNucleon A  ER/(2 ReducedMass[m\[Chi],A mNucleon ]^2)];


(* ::Subsection::Initialization:: *)
(*(*(*(*(*Detector Modules*)*)*)*)*)


(* ::Text::Initialization:: *)
(*(*(*(*(*Detector data for a given time*)*)*)*)*)


(* ::Input::Initialization:: *)
DetectorPosition[n_,{lat_,lon_},depth_]:=Module[{\[Theta],\[Phi],r},
r=r\[Earth]-depth;
\[Theta]=\[Pi]/2-lat;
\[Phi]=(2\[Pi])/86400 LAST[n,lon];
Return[Equat2Gal[SphericalCoordinates[r,\[Theta],\[Phi]],n/36525]]
]


(* ::Text:: *)
(*\[Theta] angle as a function of the LAST*)


(* ::Input::Initialization:: *)
DetectorAngle[last_,date_,{lat_?NumericQ,lon_?NumericQ}]:=Module[{\[Theta],\[Phi],x\[Earth],v\[Earth],n},
n=FractionalDays[date,{}];
v\[Earth]=EarthVelocity[n];
\[Theta]=\[Pi]/2-lat;
\[Phi]=(2\[Pi])/86400 last;
x\[Earth]=Equat2Gal[SphericalCoordinates[1,\[Theta],\[Phi]],n/36525];
Return[InUnits[ArcCos[x\[Earth].v\[Earth]/(Norm[x\[Earth]]Norm[v\[Earth]])],\[Degree]]]
]
DetectorAngle2[n_,{lat_?NumericQ,lon_?NumericQ}]:=Module[{\[Theta],\[Phi],x\[Earth],v\[Earth]},
v\[Earth]=EarthVelocity[n];
\[Theta]=\[Pi]/2-lat;
\[Phi]=(2\[Pi])/86400 LAST[n,lon];
x\[Earth]=Equat2Gal[SphericalCoordinates[1,\[Theta],\[Phi]],n/36525];
Return[InUnits[ArcCos[x\[Earth].v\[Earth]/(Norm[x\[Earth]]Norm[v\[Earth]])],\[Degree]]]
]


End[]


(* ::Input::Initialization:: *)
Print[
"The DaMaSCUStoolbox Mathematica package is part of DaMaSCUS (available at https://github.com/temken/DaMaSCUS). It is a collection of useful functions and serves as a little tool box for plotting DaMaSCUS results. For a list of functions type ?DaMaSCUS_Toolbox`*"]
EndPackage[]
