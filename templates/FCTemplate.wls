#!/usr/bin/env wolframscript


rootdir = Directory[];
intdif = rootdir <> "/Intdif2008.m";
job = "initialize";
Get[intdif];


n = 30; 
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

re = Import[rootdir <> "/variables.csv"]

rdisp = 0.01;
adisp = 0.02;

rr = Import[rootdir <> "/disps.csv"]
   
conv = HART;

A = Import[rootdir <> "/eigen.csv"]

sym = {};
U = Inverse[A];
normU = U;
For[i = 1, i <= Length[U], i++,
 normU[[i]] = Normalize[U[[i]]];
 ]

INTC[] =
  Dot[normU,
(*Cartesian coordinates in bohr*)


xref = Import[rootdir <> "/cartesians.csv"]
(*Difference should be close to zero *)
Print[INTC @@ Flatten[xref BOHR]]
Print[INTC @@ Flatten[xref BOHR] - re] 

(* Print[Norm[INTC @@ Flatten[xref BOHR] - re]]
norm_thresh = 1.0*^-4;
Print[1.0*^-4 > Norm[INTC @@ Flatten[xref BOHR] - re]]  
CheckNorm[x_] := If[1.0*^-4 > x, Print[Norm criteria met], Throw[Norm criteria not met]];
Check[Catch[CheckNorm[Norm[INTC @@ Flatten[xref BOHR] - re]]],err]
Check[Catch[Throw[Test]],err]
Check[Catch[CheckNorm[1.0*^-5]],err]
*)

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
