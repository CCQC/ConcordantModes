#!/usr/bin/env wolframscript


rootdir = Directory[];
intdif = rootdir <> "/Intdif2008.m";
job = "initialize";
Get[intdif];

xref = Import[rootdir <> "/cartesians.csv"];
n = 3 Length[xref] - 6; 
fcTrunc = True;

mass = {mH,mO,mH};
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

rr = {rdisp,rdisp,adisp};
conv = HART;

A = Import[rootdir <> "/eigen.csv"];

Clear[INTC, x1, y1, z1, x2, y2, z2, x3, y3, z3];
sym = {};
U = Inverse[A];
normU = U;
For[i = 1, i <= Length[U], i++,
 normU[[i]] = Normalize[U[[i]]];
 ]


INTC[x1_, y1_, z1_, x2_, y2_, z2_, x3_, y3_, z3_] =
  Dot[normU,
   {Rab[{x2,y2,z2},{x1,y1,z1}],
    Rab[{x3,y3,z3},{x2,y2,z2}],
    ArcCos[Cos\[Theta]abc[{x3,y3,z3},{x2,y2,z2},{x1,y1,z1}]]
    }];
(*Cartesian coordinates in bohr*)



(*Difference should be close to zero *)
Print[INTC @@ Flatten[xref BOHR]]
Print[INTC @@ Flatten[xref BOHR] - re]

job = "disp";
Get[intdif];
Print["Without symmetry ", Length[dispmin], " displaced geometries are required\n"]
job = "Cart";
Cartfunc = False;
geomDir = rootdir;
geomFile = "dispcart";
Get[intdif];

job = "disp";
Get[intdif];
Print["With symmetry ", Length[dispmin], " displaced geometries are required\n"]

job = "Cart";
Get[intdif];
