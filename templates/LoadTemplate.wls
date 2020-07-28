#!/usr/bin/env wolframscript


rootdir = Directory[];
intdif = rootdir <> "/Intdif2008.m";
job = "initialize";
Get[intdif];

xref = Import[rootdir <> "/cartesians.csv"];
n = 3 Length[xref] - 6; 
fcTrunc = True;

mass = Import[rootdir <> "/masses.csv"][[1]]
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

rdisp = 0.01;
adisp = 0.02;
(* I may just have to write this one out explicitly as well, but this shouldn't be too big of a burden *)
rr = Import[rootdir <> "/disps.csv"]

conv = HART;

A = Import[rootdir <> "/eigen.csv"];

sym = {};
U = Inverse[A];
normU = U;
For[i = 1, i <= Length[U], i++,
 normU[[i]] = Normalize[U[[i]]];
 ]

INTC[] =
  Dot[normU,
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
