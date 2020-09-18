#!/usr/bin/env wolframscript


rootdir = Directory[];
intdif = rootdir <> "/Intdif2008.m";
job = "initialize";
data = rootdir <> "/e.m";
Get[intdif];
xref = Import[rootdir <> "/cartesians.csv"];
n = 3 Length[xref] - 6; 
fcTrunc = True;

mass = {mC,mO,mH,mH,mH,mH};
nder = 2;  (*quadratic etc*)
tder = 1;  (*tder=1--> Hii, tder=2-->Hij*)  
dispder = 0;
epres = 10;
favg = False;
sym = {};
an = (\[Pi]/180);
cb = 0.5291772085936;

u1 = 1/Sqrt[2];

re = Import[rootdir <> "/symVariables.csv"][[1]];

rdisp = 0.005;
adisp = 0.01;

rr = {rdisp,rdisp,rdisp,rdisp,rdisp,adisp,adisp,adisp,adisp,adisp,adisp,adisp};   
conv = HART;

A = Import[rootdir <> "/eigen.csv"]

Clear[INTC, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6];
sym = {};
U = Inverse[A];
normU = U;
For[i = 1, i <= Length[U], i++,
 normU[[i]] = Normalize[U[[i]]];
 ]


INTC[x1_, y1_, z1_, x2_, y2_, z2_, x3_, y3_, z3_, x4_, y4_, z4_, x5_, y5_, z5_, x6_, y6_, z6_] =
  Dot[normU,
   {Rab[{x2,y2,z2},{x1,y1,z1}],
    Rab[{x3,y3,z3},{x1,y1,z1}],
    Rab[{x4,y4,z4},{x1,y1,z1}],
    Rab[{x5,y5,z5},{x1,y1,z1}],
    Rab[{x6,y6,z6},{x2,y2,z2}],
    ArcCos[Cos\[Theta]abc[{x3,y3,z3},{x1,y1,z1},{x2,y2,z2}]],
    ArcCos[Cos\[Theta]abc[{x4,y4,z4},{x1,y1,z1},{x2,y2,z2}]],
    ArcCos[Cos\[Theta]abc[{x5,y5,z5},{x1,y1,z1},{x2,y2,z2}]],
    ArcCos[Cos\[Theta]abc[{x6,y6,z6},{x2,y2,z2},{x1,y1,z1}]],
    ArcTan[Tan\[Tau]abcd[{x4,y4,z4},{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}]]+\[Pi],
    ArcTan[Tan\[Tau]abcd[{x5,y5,z5},{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}]]+\[Pi],
    ArcTan[Tan\[Tau]abcd[{x6,y6,z6},{x2,y2,z2},{x1,y1,z1},{x3,y3,z3}]]+\[Pi]
    }];
(*Cartesian coordinates in bohr*)



(*Difference should be close to zero *)

job = "disp";
Get[intdif];
(* Print["Without symmetry ", Length[dispmin], " displaced geometries are required\n"] *)
job = "Cart";
Cartfunc = False;
geomDir = rootdir;
geomFile = "dispcart";
Get[intdif];

job = "disp";
Get[intdif];
(* Print["With symmetry ", Length[dispmin], " displaced geometries are required\n"] *)

job = "Cart";
Get[intdif];


(* Force Constant Code *)
Clear[Edata];
theory = "CCSD(T)";
Get[data];
job = "trans";
Get[intdif];

Edata = EdataS;
job = "loaddata";
Get[intdif];
Print["The energies at the displaced geometries:\n", Edispmin // MatrixForm];

job = "dif";
Get[intdif];

job = "fcout";
fcDir = rootdir;
fcFile = "fc";
fout = fres;
Get[intdif];
