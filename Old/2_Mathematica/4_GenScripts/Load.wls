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

Clear[INTC, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6];
sym = {};
U = Inverse[A];
normU = U;
For[i = 1, i <= Length[U], i++,
 normU[[i]] = Normalize[U[[i]]];
 ]


INTC[x1_, y1_, z1_, x2_, y2_, z2_, x3_, y3_, z3_, x4_, y4_, z4_, x5_, y5_, z5_, x7_, y7_, z7_] =
  Dot[normU,
  {Rab[{x2, y2, z2}, {x1, y1, z1}],
    Rab[{x3, y3, z3}, {x2, y2, z2}],
    Rab[{x4, y4, z4}, {x2, y2, z2}],
    Rab[{x5, y5, z5}, {x2, y2, z2}],
    Rab[{x6, y6, z6}, {x2, y2, z2}],
    ArcCos[Cos\[Theta]abc[{x3, y3, z3}, {x2, y2, z2}, {x1, y1, z1}]],
    ArcCos[Cos\[Theta]abc[{x4, y4, z4}, {x2, y2, z2}, {x1, y1, z1}]],
    ArcCos[Cos\[Theta]abc[{x5, y5, z5}, {x2, y2, z2}, {x1, y1, z1}]],
    ArcCos[Cos\[Theta]abc[{x5, y5, z5}, {x2, y2, z2}, {x1, y1, z1}]],
    ArcTan[Tan\[Tau]abcd[{x4, y4, z4}, {x2, y2, z2}, {x1, y1, z1}, {x3, y3, z3}]] + \[Pi],
    ArcTan[Tan\[Tau]abcd[{x5, y5, z5}, {x2, y2, z2}, {x1, y1, z1}, {x3, y3, z3}]],
    ArcTan[Tan\[Tau]abcd[{x6, y6, z6}, {x2, y2, z2}, {x1, y1, z1}, {x3, y3, z3}]]
    }];
(*Cartesian coordinates in bohr*)


xref = Import[rootdir <> "/cartesians.csv"]
(*Difference should be close to zero *)
(* Print[INTC @@ Flatten[xref BOHR]] *)
(* Print[INTC @@ Flatten[xref BOHR] - re] *)
Print[Norm[INTC @@ Flatten[xref BOHR] - re]]
norm_thresh = 1.0*^-4;
Print[1.0*^-4 > Norm[INTC @@ Flatten[xref BOHR] - re]]

(* Giraffe_general: I need to figure out how to properly throw errors in mathematica... but don't spend too much time with it *)
CheckNorm[x_] := If[1.0*^-4 > x, Print[Norm criteria met], Throw[Norm criteria not met]];
Check[Catch[CheckNorm[Norm[INTC @@ Flatten[xref BOHR] - re]]],err]
Check[Catch[Throw[Test]],err]
Check[Catch[CheckNorm[1.0*^-5]],err]


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
