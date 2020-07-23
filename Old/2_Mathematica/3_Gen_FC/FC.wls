#!/usr/bin/env wolframscript

rootdir = Directory[];
data = rootdir <> "/e.m";
intdif = rootdir <> "/Intdif2008.m";
job = "initialize";
Get[intdif];

n = 30; 
fcTrunc = True;
mass = {mC, mC, mC, mO, mO, mH, mH, mH, mH, mH, mH, mH};  
nder = 2;  (*quadratic etc*)
tder = 1;  (*tder=1--> Hii, tder=2-->Hij*)  
dispder = 0;
epres = 10; 
favg = False;
sym = {};
an = (\[Pi]/180);
cb = 0.5291772085936;

u1 = 1/Sqrt[2];


(* Giraffe1 *)
re = {2.9727515231, 1.2917344808, 0.7924709411, 
   1.8016850869, 0.2528610525, 3.9115028969, 3.5068468496, 
   2.3432592342, 1.8612640334, -0.6661950188, 0.5646526723, 
   0.9432805483, 1.3564492954, 0.0794605796, 0.5662277751, 
   3.4938121890, -0.2501789570, 2.1671852390, 
   2.0960217173, -1.5117140996, 3.2612967640, 0.1578639860, 
   1.0890351379, 0.7848110744, 2.3242168046, 1.1188126015, 
   0.1651575280, -0.2265920988, 0.0548913631, 0.2047233579};


rdisp = 0.01;
adisp = 0.02;
rr = {(* Bond distances*) rdisp, rdisp, rdisp, rdisp, rdisp, rdisp, rdisp, rdisp, rdisp, rdisp, rdisp,
   (* Bond angles*) adisp, adisp, adisp, adisp, adisp, adisp, adisp, adisp, adisp, adisp,
   (* Bond dihedrals*)adisp, adisp, adisp, adisp, adisp, adisp, adisp, adisp, adisp};
conv = HART;

A = Import[rootdir <> "/eigen.csv"]

Clear[INTC, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, 
  z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, x9, y9, z9, x10, y10, z10, 
  x11, y11, z11, x12, y12, z12];
sym = {};
U = Inverse[A];
normU = U;
For[i = 1, i <= Length[U], i++,
 normU[[i]] = Normalize[U[[i]]];
 ]

INTC[x1_, y1_, z1_, x2_, y2_, z2_, x3_, y3_, z3_, x4_, y4_, z4_, x5_, 
   y5_, z5_, x6_, y6_, z6_, x7_, y7_, z7_, x8_, y8_, z8_, x9_, y9_, 
   z9_, x10_, y10_, z10_, x11_, y11_, z11_, x12_, y12_, z12_] =
  Dot[normU,
   {Rab[{x1, y1, z1}, {x2, y2, z2}],
    Rab[{x2, y2, z2}, {x3, y3, z3}],
    Rab[{x2, y2, z2}, {x4, y4, z4}],
    Rab[{x4, y4, z4}, {x5, y5, z5}],
    Rab[{x1, y1, z1}, {x6, y6, z6}],
    Rab[{x1, y1, z1}, {x7, y7, z7}],
    Rab[{x1, y1, z1}, {x8, y8, z8}],
    Rab[{x2, y2, z2}, {x9, y9, z9}],
    Rab[{x3, y3, z3}, {x10, y10, z10}],
    Rab[{x3, y3, z3}, {x11, y11, z11}],
    Rab[{x3, y3, z3}, {x12, y12, z12}],
    ArcCos[Cos\[Theta]abc[{x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3}]],
    ArcCos[Cos\[Theta]abc[{x3, y3, z3}, {x2, y2, z2}, {x4, y4, z4}]],
    ArcCos[Cos\[Theta]abc[{x2, y2, z2}, {x4, y4, z4}, {x5, y5, z5}]],
    ArcCos[Cos\[Theta]abc[{x2, y2, z2}, {x1, y1, z1}, {x6, y6, z6}]],
    ArcCos[Cos\[Theta]abc[{x2, y2, z2}, {x1, y1, z1}, {x7, y7, z7}]],
    ArcCos[Cos\[Theta]abc[{x2, y2, z2}, {x1, y1, z1}, {x8, y8, z8}]],
    ArcCos[Cos\[Theta]abc[{x4, y4, z4}, {x2, y2, z2}, {x9, y9, z9}]],
    ArcCos[
     Cos\[Theta]abc[{x2, y2, z2}, {x3, y3, z3}, {x10, y10, z10}]],
    ArcCos[
     Cos\[Theta]abc[{x2, y2, z2}, {x3, y3, z3}, {x11, y11, z11}]],
    ArcCos[
     Cos\[Theta]abc[{x2, y2, z2}, {x3, y3, z3}, {x12, y12, z12}]],
    ArcTan[
      Tan\[Tau]abcd[{x4, y4, z4}, {x2, y2, z2}, {x3, y3, z3}, {x1, y1,
         z1}]] + \[Pi],
    ArcTan[
      Tan\[Tau]abcd[{x5, y5, z5}, {x4, y4, z4}, {x2, y2, z2}, {x1, y1,
         z1}]] + \[Pi],
    ArcTan[
     Tan\[Tau]abcd[{x6, y6, z6}, {x1, y1, z1}, {x2, y2, z2}, {x3, y3, 
       z3}]],
    ArcTan[
      Tan\[Tau]abcd[{x7, y7, z7}, {x1, y1, z1}, {x2, y2, z2}, {x3, y3,
         z3}]] + \[Pi],
    ArcTan[
     Tan\[Tau]abcd[{x8, y8, z8}, {x1, y1, z1}, {x2, y2, z2}, {x3, y3, 
       z3}]],
    ArcTan[
     Tan\[Tau]abcd[{x9, y9, z9}, {x2, y2, z2}, {x4, y4, z4}, {x5, y5, 
       z5}]],
    ArcTan[
     Tan\[Tau]abcd[{x10, y10, z10}, {x3, y3, z3}, {x2, y2, z2}, {x1, 
       y1, z1}]],
    ArcTan[
     Tan\[Tau]abcd[{x11, y11, z11}, {x3, y3, z3}, {x2, y2, z2}, {x1, 
       y1, z1}]],
    ArcTan[
      Tan\[Tau]abcd[{x12, y12, z12}, {x3, y3, z3}, {x2, y2, z2}, {x1, 
        y1, z1}]] + \[Pi]
    }];
(*Cartesian coordinates in bohr*)

xref = {{-3.101453947, 1.586241108, 0.020127067}, {-0.864067001, -0.071005451, -0.656151092}, {-1.118044488, -2.779891761, 0.246215890}, {1.311487864, 1.104949528, 0.589001118}, 
        {3.449295345, 0.120060183, -0.251828404}, {-4.804867472, 0.847236857, -0.876647946}, {-2.806351805, 3.522970736, -0.618661947}, {-3.392813348, 1.597259180, 2.062741060}, 
        {-0.460013342, -0.004348894, -2.677343878}, {-2.697536421, -3.692575647, -0.717730545}, {-1.471602235, -2.823825065, 2.279176614}, {0.605243032, -3.830466668, -0.161330336}};

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
(* Giraffe2 *)





















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
