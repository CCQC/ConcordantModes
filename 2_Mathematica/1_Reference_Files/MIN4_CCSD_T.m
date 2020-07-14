(* ::Package:: *)

(* ::Input:: *)
(*(* intdif="/Users/marcusbartlett/Documents/Combustion/i-propyl+o2/MIN4/CCSD_T_cc-pVTZ/Intdif2008.m"; *)*)
(*(* rootdir = "/Users/marcusbartlett/Documents/Combustion/i-propyl+o2/MIN4/CCS\.00\.00\.00\.00\.00\.00\.00\.00\.00\.00\.00\.00D_T_cc-pVTZ"; *)*)
(*rootdir = NotebookDirectory[];*)
(*data=rootdir<>"/e.m";*)
(*intdif=rootdir<>"Intdif2008.m";*)
(*(*data=StringJoin[rootdir,"/e.m"];*)*)
(*job="initialize";*)
(*Get[intdif];*)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Input the data*)


(* ::Input:: *)
(*(*test run in MP2/cc-pvdz, harmonic frequencies*)(*dispmin to check to see if displacements are correct*)*)
(*n=30; *)
(*fcTrunc=True;*)
(*mass={mC,mC,mC,mO,mO,mH,mH,mH,mH,mH,mH,mH};  *)
(*nder=2;  (*quadratic etc*)*)
(*tder=1;  (*tder=1--> Hii, tder=2-->Hij*)  *)
(*dispder=0;*)
(*epres=10; *)
(*favg = False;*)
(*sym={};*)
(*an=(\[Pi]/180);*)
(*cb=0.5291772085936;*)
(**)
(*u1=1/Sqrt[2];*)
(**)
(**)
(*(*re={(* Bond distances*)-2.0040224329,3.6955720160,-1.3907280147,1.8408086874,2.6441105893,2.4398578819,3.3218796943,2.8279196519,1.5886734822,2.9130208604,-0.3416054802,*)
(*(*Bond angles*)137.03773859Degree,66.523852973Degree,-62.399915715Degree,23.979450658Degree,102.24726828Degree,174.99341355Degree,9.3072422303Degree,34.702507519Degree,177.95018545Degree,1.8206669122Degree,*)
(*(*Bond dihedrals*)134.883386Degree,72.660264569Degree,132.71947343Degree,99.117385732Degree,-30.362592124Degree,-19.820888758Degree,32.492860075Degree,85.908422487Degree,6.9041207819Degree};*)*)
(**)
(*re=(*{1.5162307539,1.5168878326,1.4651916691,1.3226490930,1.0911771331,1.0904208379,1.0918612577,1.0913018980,1.0918222415,1.0921919113,1.0895811521,113.7722993647Degree,109.2331516136Degree,110.6960404369Degree,109.5802853575Degree,111.0158522775Degree,110.3046905258Degree,105.8473828999Degree,109.8029079510Degree,110.2040798914Degree,110.2253597264Degree,117.1687233234Degree,192.1312386253Degree,-61.6029843004Degree,178.3394563794Degree,57.6717850786Degree,-49.7831691534Degree,64.6814194094Degree,-54.9559195907Degree,184.5493688843Degree}*){2.9727515231,1.2917344808,0.7924709411,1.8016850869,0.2528610525,3.9115028969,3.5068468496,2.3432592342,1.8612640334,-0.6661950188,0.5646526723,0.9432805483,1.3564492954,0.0794605796,0.5662277751,3.4938121890,-0.2501789570,2.1671852390,2.0960217173,-1.5117140996,3.2612967640,0.1578639860,1.0890351379,0.7848110744,2.3242168046,1.1188126015,0.1651575280,-0.2265920988,0.0548913631,0.2047233579};*)
(**)
(**)
(*rdisp=0.01;*)
(*adisp=0.02;*)
(*rr={(* Bond distances*)rdisp,rdisp,rdisp,rdisp,rdisp,rdisp,rdisp,rdisp,rdisp,rdisp,rdisp,*)
(*(* Bond angles*)       adisp,adisp,adisp,adisp,adisp,adisp,adisp,adisp,adisp,adisp,*)
(*(* Bond dihedrals*)adisp,adisp,adisp,adisp,adisp,adisp,adisp,adisp,adisp};*)
(*conv=HART;*)
(**)
(*A=Import[NotebookDirectory[]<>"eigen.csv"]*)
(*(*A={{-0.000364,-0.000022,-0.000225,-0.001211,0.001510,0.006131,0.052791,0.062252,0.147722,-0.003109,-0.140507,0.021094,-0.111035,-0.109345,0.154864,0.072509,-0.058402,0.063054,-0.214470,-0.068644,-0.017763,-0.043547,-0.020250,0.024317,-0.044492,-0.032321,-0.000855,0.007772,-0.003431,-0.003527},{0.000801,0.000305,-0.000126,0.005108,0.004533,0.021169,0.004165,0.101120,0.175914,0.072123,0.065773,0.011022,0.087170,0.216104,0.035810,0.001235,0.056948,-0.211670,0.048892,0.068607,-0.004161,-0.013898,-0.028323,-0.046287,-0.020292,-0.030411,-0.003332,0.006222,-0.001816,-0.007806},{-0.001633,0.001186,0.003319,0.000334,-0.009173,0.013149,0.063704,0.218994,-0.162363,-0.048161,0.063397,-0.186664,0.092626,-0.098205,-0.000595,-0.064690,-0.024166,0.022319,-0.000643,-0.001911,0.001227,0.013052,0.043882,0.000417,0.001671,-0.019771,-0.002054,0.005401,0.000645,-0.000266},{0.000551,-0.001199,-0.000063,-0.006626,-0.008306,0.005328,0.012288,0.016712,-0.011561,-0.009023,-0.011014,0.289449,0.178761,-0.082581,0.035226,-0.016006,-0.003343,0.002675,0.007105,-0.000859,-0.001541,0.000752,0.002133,-0.000037,0.000162,0.000258,-0.000297,0.000016,0.000033,0.000161},{0.000182,-0.000180,0.000580,-0.001231,-0.001064,0.000692,0.003176,-0.001135,0.005990,0.000301,0.004074,-0.005612,0.002550,-0.010274,0.005896,-0.000115,-0.008865,0.003866,0.005460,-0.000907,0.008859,-0.007145,-0.007698,-0.262276,0.512152,-0.030996,-0.300797,-0.560693,-0.580999,-0.023368},{-0.000073,-0.000368,0.000135,-0.000236,0.001287,0.000723,0.000896,0.001295,0.001356,-0.005361,-0.003916,0.001629,-0.003134,0.005235,0.009306,0.007447,-0.002069,-0.003093,0.005800,-0.010198,-0.001008,0.010666,0.005008,-0.237475,0.494664,-0.041777,-0.093366,-0.215056,0.844803,0.059620},{0.000083,0.000193,-0.000277,0.000133,-0.000196,-0.000191,-0.000532,0.003892,0.000727,0.004833,-0.005569,0.003860,-0.006425,-0.002431,-0.006695,0.006646,-0.008644,0.002957,0.003185,0.008873,-0.012067,-0.001602,0.005224,-0.279832,0.556452,0.139519,0.368938,0.697666,-0.208015,-0.038284},{-0.001682,-0.000279,0.000096,0.001024,0.001119,0.002299,0.002279,-0.003963,0.004057,0.000588,-0.001332,-0.002162,0.003324,0.000824,0.004372,0.010009,0.003112,-0.005805,-0.004746,0.001360,-0.002314,-0.004877,0.001460,-0.001535,-0.053364,1.015792,0.060019,-0.178245,0.037914,0.061560},{0.000201,-0.000101,0.000251,0.000834,-0.000658,0.002115,-0.000805,0.002250,0.002110,0.001287,-0.006604,-0.005659,0.011431,0.004278,0.000883,-0.005118,0.007562,0.005475,0.003694,-0.000156,-0.010524,0.003600,-0.009641,0.545738,0.267599,-0.065342,0.669733,-0.354376,0.002500,-0.356571},{-0.000194,0.000124,-0.000211,-0.000506,-0.000166,-0.000653,-0.000282,0.002564,0.005821,-0.003909,0.004105,0.003483,-0.003580,0.002664,-0.008320,0.001239,0.011191,0.002177,0.004041,-0.006361,0.004592,-0.011656,0.007937,0.582452,0.288459,0.128429,-0.643029,0.308582,0.033794,-0.355221},{0.001592,0.000261,-0.000187,0.001154,0.000873,0.000294,0.001924,0.001026,0.000312,0.006157,0.003886,0.002017,-0.001911,0.005897,0.008543,0.006556,0.002920,0.005527,-0.003289,0.007663,0.007437,0.005344,0.002541,0.426432,0.219343,-0.038086,0.020923,0.022965,-0.066981,0.915583},{0.005732,-0.003803,-0.002401,-0.018485,0.211810,0.109285,0.085390,-0.008471,-0.062706,-0.015045,0.013621,-0.051977,0.050447,-0.140945,-0.234542,-0.105044,-0.007334,0.027166,0.038418,0.003833,0.000855,0.002034,0.011191,0.004943,0.003256,0.066781,0.005383,0.003889,-0.048072,-0.056259},{-0.013114,0.013047,0.012039,0.103359,-0.078857,0.209349,-0.106952,-0.046092,0.015746,-0.010855,0.024618,0.142008,-0.173333,-0.091334,-0.008991,0.094614,-0.071955,0.053858,-0.066760,-0.042669,0.020654,-0.010108,-0.002648,0.000080,-0.004213,0.050813,-0.039650,0.013101,0.000857,0.032887},{-0.008707,0.011461,-0.019233,0.105563,-0.003275,-0.037997,0.222790,-0.191496,0.022812,0.030517,-0.005930,-0.037722,-0.075757,0.123591,-0.131912,-0.078431,0.066662,-0.042704,0.064162,0.037045,0.017300,0.021148,-0.005840,-0.001669,-0.000470,0.032595,0.002359,-0.006638,0.001044,0.002523},{-0.000809,0.009925,-0.008749,0.050930,0.065307,-0.033566,-0.054903,0.241001,-0.239248,-0.087584,-0.384437,0.224448,-0.247291,0.337546,-0.071993,0.259082,-0.198083,0.012260,0.533837,0.064564,0.160409,-0.077666,-0.191987,-0.015286,0.029874,0.030484,0.029876,0.043775,0.059970,0.003116},{-0.002902,0.020760,-0.002040,-0.024190,-0.081841,-0.051089,0.000967,0.008853,0.192188,0.465672,0.117962,-0.092598,0.061812,-0.350913,-0.281667,0.100915,-0.345873,0.105183,0.449190,-0.197100,0.013378,0.282813,0.155785,-0.013924,0.031618,0.039092,0.010097,0.011877,-0.078835,-0.003307},{0.004691,-0.032458,0.010704,-0.033528,0.011362,0.091157,0.064214,-0.214405,0.101294,-0.383455,0.206771,-0.141570,0.087307,-0.066370,0.496601,0.005657,-0.124400,0.000099,0.518307,0.253178,-0.219904,-0.027197,0.193455,-0.012197,0.029759,-0.059919,-0.042427,-0.061729,0.022839,-0.000181},{-0.007399,0.006059,-0.023349,0.021540,0.012788,-0.081053,0.003611,-0.135117,0.048779,0.008423,-0.027480,-0.044770,0.298010,-0.011753,0.087541,0.870759,0.195805,-0.127406,-0.081403,-0.008318,-0.039128,-0.084186,-0.117008,-0.001316,-0.000458,-0.063367,-0.003771,0.011889,-0.005176,-0.007523},{0.005379,-0.006610,-0.004732,-0.029347,0.051928,-0.110217,0.075136,0.108532,0.010937,0.005336,0.496288,0.271084,-0.398604,0.027614,0.031439,0.190673,0.298549,0.501216,0.058642,-0.109619,-0.162207,0.040423,-0.205357,0.027995,0.012459,0.033689,-0.060330,0.026813,-0.001662,0.046615},{0.004066,0.017683,0.023707,0.009376,0.017718,0.125274,0.009645,-0.121485,-0.206447,0.382784,-0.200523,-0.137499,0.251348,0.239720,0.356771,-0.042175,0.092243,0.533689,0.028258,-0.209202,0.020097,-0.231549,0.227772,0.023569,0.015276,-0.060126,0.063030,-0.019918,-0.007342,0.043012},{-0.005740,-0.010365,-0.019363,0.023613,-0.075076,-0.007366,-0.084561,0.055579,0.266107,-0.366540,-0.265235,-0.122762,0.202530,-0.076367,-0.355178,-0.077577,0.384333,0.455207,0.129760,0.174650,0.195287,0.192872,0.124260,0.036514,0.019784,0.033833,0.001435,-0.008512,0.006731,-0.075030},{0.011566,-0.022832,0.035603,-0.127155,-0.061686,0.270856,0.095004,-0.139143,-0.021927,0.010524,-0.006472,0.069232,-0.220390,0.125340,-0.313645,0.029219,0.082009,-0.046392,0.052112,0.026527,0.005234,0.026375,-0.010657,0.000463,-0.002699,0.135536,0.004418,0.044196,0.023929,-0.006864},{0.321505,-0.010607,-0.040878,-0.073166,0.095747,0.024909,0.047432,-0.031297,-0.031143,0.039151,-0.020340,0.003951,-0.020925,-0.060983,-0.135213,0.018307,-0.043819,0.038630,-0.018893,0.007233,-0.006316,-0.002387,0.000729,-0.001051,0.000452,0.051201,0.022242,0.029669,-0.039327,-0.001108},{0.018480,0.503976,-0.232768,-0.138767,-0.102135,0.158378,-0.050287,-0.000835,0.063307,0.247857,0.000852,0.157095,-0.187495,-0.074579,-0.412691,0.124684,-0.116251,0.138566,-0.079301,0.510972,-0.272430,-0.380081,0.036604,-0.001750,0.000688,0.096036,-0.010880,0.090019,-0.053494,0.004174},{0.019868,0.477326,-0.223155,-0.169033,-0.094459,0.233993,0.006511,-0.209970,0.149483,-0.115919,0.221629,-0.014877,-0.032714,-0.120307,0.125868,-0.047979,0.052840,-0.033292,0.114514,-0.222596,0.546336,-0.107645,-0.415918,0.000928,-0.002905,0.045736,-0.080977,-0.036886,-0.012660,0.008633},{0.018477,0.488920,-0.233717,-0.118109,-0.033399,0.200564,-0.049340,-0.000313,-0.092278,-0.192559,-0.136253,0.253754,-0.276570,0.312409,-0.024251,0.193281,0.165288,-0.168607,-0.035767,-0.295905,-0.226126,0.450500,0.357680,0.002114,-0.004357,0.068361,-0.027619,0.058850,0.096361,0.012796},{0.324704,0.004178,-0.041837,-0.091917,0.019053,-0.130575,0.016975,0.047948,-0.102112,0.338333,-0.170432,-0.050869,-0.101154,-0.200542,0.169915,-0.113277,0.701249,-0.352074,0.310185,0.225978,-0.001082,0.034155,-0.019616,-0.006997,0.002876,-0.015548,-0.003328,0.002752,-0.001702,-0.000842},{0.016493,0.271449,0.489380,0.124954,0.119625,-0.160558,0.057183,0.000916,-0.152775,0.197344,0.058497,-0.132973,0.169640,0.033672,0.407523,-0.022760,-0.183732,0.152606,-0.176523,0.419040,0.196368,0.481112,-0.119397,-0.008025,-0.001907,-0.094455,0.029367,-0.066743,-0.027402,0.072597},{0.010496,0.264043,0.478539,0.143156,0.056735,-0.165300,-0.016250,0.044015,0.069634,-0.158423,-0.201157,-0.274162,0.363641,-0.113209,-0.035351,-0.277990,0.056199,-0.009981,0.163020,-0.304701,-0.485273,-0.108825,-0.359325,0.007259,0.004913,-0.070239,0.028719,-0.075332,-0.014820,-0.083711},{0.012241,0.253117,0.460246,0.116648,0.101994,-0.263745,0.048658,0.123743,0.065321,-0.154831,0.252394,0.024418,-0.086710,-0.143113,-0.006777,0.029312,0.084674,-0.112706,0.044937,-0.122548,0.337935,-0.383958,0.516473,0.001694,0.001749,-0.044684,-0.087055,-0.014642,-0.017389,-0.004291}};*)*)
(**)


(* ::Input:: *)
(*Clear[INTC,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12];*)
(*sym={};*)
(*U=Inverse[A];*)
(*normU=U;*)
(*For[i=1,i<=Length[U],i++,*)
(*normU[[i]]=Normalize[U[[i]]];*)
(*]*)
(**)
(*INTC[x1_,y1_,z1_,x2_,y2_,z2_,x3_,y3_,z3_,x4_,y4_,z4_,x5_,y5_,z5_,x6_,y6_,z6_,x7_,y7_,z7_,x8_,y8_,z8_,x9_,y9_,z9_,x10_,y10_,z10_,x11_,y11_,z11_,x12_,y12_,z12_]=*)
(*Dot[normU,*)
(*{Rab[{x1,y1,z1},{x2,y2,z2}],*)
(*Rab[{x2,y2,z2},{x3,y3,z3}],*)
(*Rab[{x2,y2,z2},{x4,y4,z4}],*)
(*Rab[{x4,y4,z4},{x5,y5,z5}],*)
(*Rab[{x1,y1,z1},{x6,y6,z6}],*)
(*Rab[{x1,y1,z1},{x7,y7,z7}],*)
(*Rab[{x1,y1,z1},{x8,y8,z8}],*)
(*Rab[{x2,y2,z2},{x9,y9,z9}],*)
(*Rab[{x3,y3,z3},{x10,y10,z10}],*)
(*Rab[{x3,y3,z3},{x11,y11,z11}],*)
(*Rab[{x3,y3,z3},{x12,y12,z12}],*)
(*ArcCos[Cos\[Theta]abc[{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}]],*)
(*ArcCos[Cos\[Theta]abc[{x3,y3,z3},{x2,y2,z2},{x4,y4,z4}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x4,y4,z4},{x5,y5,z5}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x1,y1,z1},{x6,y6,z6}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x1,y1,z1},{x7,y7,z7}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x1,y1,z1},{x8,y8,z8}]],*)
(*ArcCos[Cos\[Theta]abc[{x4,y4,z4},{x2,y2,z2},{x9,y9,z9}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x3,y3,z3},{x10,y10,z10}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x3,y3,z3},{x11,y11,z11}]],*)
(*ArcCos[Cos\[Theta]abc[{x2,y2,z2},{x3,y3,z3},{x12,y12,z12}]],*)
(*ArcTan[Tan\[Tau]abcd[{x4,y4,z4},{x2,y2,z2},{x3,y3,z3},{x1,y1,z1}]]+\[Pi],*)
(*ArcTan[Tan\[Tau]abcd[{x5,y5,z5},{x4,y4,z4},{x2,y2,z2},{x1,y1,z1}]]+\[Pi],*)
(*ArcTan[Tan\[Tau]abcd[{x6,y6,z6},{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}]],*)
(*ArcTan[Tan\[Tau]abcd[{x7,y7,z7},{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}]]+\[Pi],*)
(*ArcTan[Tan\[Tau]abcd[{x8,y8,z8},{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}]],*)
(*ArcTan[Tan\[Tau]abcd[{x9,y9,z9},{x2,y2,z2},{x4,y4,z4},{x5,y5,z5}]],*)
(*ArcTan[Tan\[Tau]abcd[{x10,y10,z10},{x3,y3,z3},{x2,y2,z2},{x1,y1,z1}]],*)
(*ArcTan[Tan\[Tau]abcd[{x11,y11,z11},{x3,y3,z3},{x2,y2,z2},{x1,y1,z1}]],*)
(*ArcTan[Tan\[Tau]abcd[{x12,y12,z12},{x3,y3,z3},{x2,y2,z2},{x1,y1,z1}]]+\[Pi]*)
(*}];*)
(*(*Cartesian coordinates in bohr*)*)
(*xref={{-3.101453947,1.586241108,0.020127067},{-0.864067001,-0.071005451,-0.656151092},{-1.118044488,-2.779891761,0.246215890},{1.311487864,1.104949528,0.589001118},{3.449295345,0.120060183,-0.251828404},{-4.804867472,0.847236857,-0.876647946},{-2.806351805,3.522970736,-0.618661947},{-3.392813348,1.597259180,2.062741060},{-0.460013342,-0.004348894,-2.677343878},{-2.697536421,-3.692575647,-0.717730545},{-1.471602235,-2.823825065,2.279176614},{0.605243032,-3.830466668,-0.161330336}};*)
(*(*Difference should be close to zero *)*)
(*INTC @@ Flatten[xref BOHR]*)
(*INTC @@ Flatten[xref BOHR]-re*)
(**)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Generate the displaced geometries*)


(* ::Input:: *)
(*job="disp";*)
(*Get[intdif];*)
(*Print["Without symmetry ",Length[dispmin]," displaced geometries are required\n"]*)
(*job="Cart";*)
(*Cartfunc=False;*)
(*geomDir=rootdir;*)
(*geomFile="dispcart";*)
(*Get[intdif];*)
(**)
(*(*job = "getsym";*)
(*Get[intdif];*)*)
(**)
(*job = "disp";*)
(*Get[intdif];*)
(*Print["With symmetry ",Length[dispmin]," displaced geometries are required\n"]*)
(**)
(*job = "Cart";*)
(*Get[intdif];*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Read inmk\[RegisteredTrademark]\[Cent]txmomjkajkqjiqjujxjjq*)


(* ::Input:: *)
(*Clear[Edata];*)
(*theory="CCSD(T)";*)
(*Get[data];*)
(*job="trans";*)
(*Get[intdif];*)


(* ::Section:: *)
(*Set up the data in internal coordinates and set-up for numerical differentiation*)


(* ::Input:: *)
(*Edata=EdataS;*)
(*job="loaddata";*)
(*Get[intdif];*)
(*Print["The energies at the displaced geometries:\n",Edispmin//MatrixForm];*)


(* ::Input:: *)
(**)
(*Table[{i,Edispmin[[1,i]]},{i,Length[Edispmin[[1]]]}]//TableForm;*)


(* ::Section:: *)
(*Numerically differentiate the data*)


(* ::Input:: *)
(*job="dif";*)
(*Get[intdif];*)


(* ::Section:: *)
(*Output the force constants in internal coordinates, in a format consistent with IntDer*)


(* ::Input:: *)
(*job="fcout";*)
(*fcDir=rootdir;*)
(*fcFile="fc";*)
(*fout=fres;*)
(*Get[intdif];*)
(**)


(* ::Section:: *)
(*Compute the harmonic frequencies from the data*)


(* ::Input:: *)
(*xref={{-3.54713827, 1.24239874,0.46699913},{-2.05884150,-0.66160490,-1.11935566},{0.13711742,-1.82640140,0.32140209},{1.93954290,0.10132312,1.11681119},{3.15033286,1.00261703,-0.82714399},{-5.09908690,2.08498945,-0.63512447},{-2.30532991,2.78082365,1.11913259},{-4.37498794,0.32831600,2.14895940},{-1.29530184,0.24588228,-2.83362465},{-3.29738016,-2.21908885,-1.75236554},{-0.47332423,-2.69094312,2.11492726},{1.18219602,-3.21907489,-0.81851538}}BOHR;*)
(*fc2=fres[[2,1]];*)
(*job="freq";*)
(*Get[intdif];*)
(*Print["The harmonic frequencies:\n",freqs//TableForm];*)
(*Print["The eigenvectors of the hessian:\n",TableForm[{Chop[Transpose[Lvib],10^-5]}]];*)
(*(*Frequencies to be matched:*)
(*{503.1538i,87.1897,174.3683,238.7436,408.5695,457.0040,476.7031,544.9953,671.5019,714.3514,741.1315,773.9302,856.3962,876.7519,953.0060,1063.0984,1085.0887,1106.8924,1115.6133,1119.7319,1167.6427,1183.7379,1222.6566,1274.5673,1278.7711,1342.5129,1407.7729,1439.4815,1593.4453,1645.5047,1776.7239,1799.6682,2028.0466,3341.1847,3356.0716,3367.8995,3389.1910,3390.2619,4117.8583}*)*)


(* ::Print:: *)
(*SequenceForm["The eigenvectors of the hessian:\n", TableForm[CompressedData["*)
(*1:eJwtl1c8EI4XxSVkZ1W07Bkpq5CuPSNbMhJKZJNRyczIyJ5l7+wRCdfeZI+M*)
(*jJCfhpRQ4t/D/+G83vvwPfd+zmE1d9C+c5iAgODQP539vy6mKL/4nn4Rk7ai*)
(*9LP2DeDE5rdp1qvC8GxdSIFmQwIaX3uoizssIgWx9HXtpE/QL7IpSPezGYO3*)
(*xDa/r9WCagL7LyLGJoyqohD0OFWJT6npQx0ky+DurMQroaIMWBY3P2HW+ATo*)
(*/Wb9s/segNQWB8WSUykwhS6Njtwox/It9doHvH7gX7z3H88vN2Reit1zZguE*)
(*FlXXnxotbGDWMEQ8IeeDXgnfz0XfTUeDU2uF+uvOoK1PMPKJpgAJxPWI3EwS*)
(*sauQvcvz5Gto1PF836kcheRURX9H6YLQdy5H+93naCR7wJ07LZOAVwRL+ppb*)
(*XYHR7NxkrbQdMAf4km6/NgU5487nHZ9lkLTdd+/JRDvGXpQ2PX4lF+w/xT78*)
(*fK0egjVDr1lodeILF0tbRv4SoCjOpGt6sAwUcfdIy0c+oiKd2FShUSEI18T0*)
(*Uan54WbkAFNNXiFy833jFJwOga85p4Kp5UOwmjqxgd9FDssPjObE2bLwUjWx*)
(*LzNdBC5v54x2fyvBoUcpNpzbpaCoJHVDNrMSjvafbH6xXYmOvvdvhayHoNS8*)
(*EYuBygt06vnxfq5DDVvv6Za4XK0BNXupRBPBShBO+1KUzVWAtLvWtspOhaj5*)
(*OXPOxcwbJ9h4E8osTPD2c3ahvFljYC0s/rImToNzQmnfv5IswufKDzP+94eg*)
(*r2ktXM93GSeMW0mabN4Bi6RohqxyK3aKsEzO2raAQHPa1/IfHWj5X9xPHvdI*)
(*1NOM14ujCwDm0pr0TscQDDBuuirqmIcCA9YdfZaeeCAcOZN1oQxmC/44LySl*)
(*Ym9YiV1/eCRY+HxRM3fNwIrTk7bPomLAurj/RnKnL9hMCUW/eZSOxLes5a72*)
(*lSGd13aVZ1oxuNd0fdZKT0fhybTUS2KR8LR8jOoXxT9uLS5HF1aTUMXhoUkH*)
(*pOM+J3E3w+n7QCcVPp6row1E8lGJZYpOWGlIv5EMp6F4UnqSUrwJaMqC3c9u*)
(*9YCbtG86meEA9tUQJf4VWkXBuvzrY9fnwUWDyWmq/CWazZxeqShqgwvr+ZaF*)
(*//bIeVxKfXIhHdUFlBmbgzzgaqQXC1tUBG4dkd/PSCtCfybVhxO2lfCmVGF9*)
(*IqsY6HTG05IE85Hyy3DnQsErIEr74G0XEI7s+S6vCCgbcN84adTnWi6sH6I1*)
(*DhPKhGuk2QGnu+tBVM+mhaGpBp649izp3X8GtieXDKTr0qBvk7GUm7YGf11g*)
(*UrA4UoFB0nLz83M+6Hzz8lazSiGq+I4rVMMLvLZCdOFNgD4eIbW6tXnZBbu5*)
(*/Jx2zBxAgCHquX1UGGieVDzTvdEFYUKW9FqedTAgKmQdxr8CTvaqz6r8luFS*)
(*10xinp4dVLkb09y/Xw5cf7/4Mh2EoxM32YSYTDDYCxCeTB9SQTXlpzMmi+64*)
(*IEQ2sfcyBh/uxcFaZzZU8yd55ITcRq83Wjal0tFw/KZ/DUVgApSX5XSrpxqC*)
(*7kzccw2CDAw6F8PDrpWNM/EBu8cP0vHPFcVMgS8OMGobbhxp+Ritboax0a0X*)
(*gz3h8ooZUTzqGoqT8S39m3/xd6XrQBIYs3fPPCB8gOk1J3qdmflAwjfaZnY7*)
(*Eu7Ovw/X0S1HYZWKo9dKCrGxacNHumseNY4ti5t8WMM2AVeepKY6yHlWT/ei*)
(*uR3MfX05/B9U4jOrK349YtnwE/gSAxul8dFfSRIBlVws1Eyw+uSeAm9bf+y/*)
(*FHABb1V2LggIwuzmBZlmvWgIk/KfPngfBrkejTqKO9no3hWnxvzLE0wFo51u*)
(*BKXDsU9zFuAfibJk4SuJZ/zxIdOrg86lLCT+8fRLvkcZPvJ5HVogHwW8z9tW*)
(*rxfkALXdeYuB0wFQ0BL0zmEsH++SCoAp4x0Q+BrcI8wkCcfuZ5cOG9CCSGPE*)
(*6DLxApx2789sVJwHGn1C5TSbZWg+3H+M+XgUjl1nk73UVoYLL/+7GmxjhlLH*)
(*CyT57ByhrdrM40GNCnKSFF/YHzUDjThXwz6vg6sS7DTthOLpwNh+Pev1VhZw*)
(*iLefkRGMBU0RFaoWakcs2tRcthfQwBo3QhOnIS3YSSKHZiEl4Gn7uS806AfN*)
(*ZCfMCdZZoYC649d5dydYywsv35wLRpuPRvx/7/qgoYrztl/TJRyrfnFd/LUF*)
(*Cr1P3RF5EgAyIo/N9gT8oNz/2vRWfhioMRAfi68tRp7SAkFZJTfQbGZ8sF2i*)
(*DfQsVOHxspHg8Hv/O+ltBxTZZmg9OhuBN9xSEmM2dLEyt7NvzE8Wgk4QDG22*)
(*0qD28tJFimktXPrbQGnyxBWCxCY4F8wCIeDEmaNVeY6Yv7Tif2NmBP6Q09FX*)
(*bNXhacGw+49Vh5BiwvOCE1MT8JxfyE48iqDk+I1vxzcJ9ERYtSzFikCArH+E*)
(*VNAcn0nVvN5yssftq/eIlxj70IJ3OtPXfQXobaM8jg/M4SmzScuy8gK4Uuuo*)
(*MTxZgYtt15lrpcIg1Fq6aeh9BOrM2IQ3sz3EqMBYz5vnsrG1/4i4otQjvCRB*)
(*brJAE4M/OBsifQkcQabeVbv0WwSSLZav9twOglPN03VVHHehUdLsE9PqbZyK*)
(*C9va1zYEeCT+MbrfEUKoJzmCf7rAAEvT1ZS+CHjLcMHPPtwAv+YdDf9o+RYg*)
(*NFmlbn0Eo0Nyfg8TvoOnu2ZVDWRZ4Jc9RhC8/BoKbv7RWixvgpPygertewjb*)
(*PnopQhxPYNZIKERJ0gcKNzftby5/hPwxSmla/U6kEAw/pJQ8h0rzdfK7RQhu*)
(*XpxvlZZSsIDma93VgD5IId+tqy/pROkcprvNu9m4WkqsEJach4L6NdsRvXGY*)
(*1HnnkVFfHv4KdZrxyLXDh8W6b2yWPcHQgHWMkd4ZydTWvUuvOIF98Bb57vYd*)
(*uPi2Qt069ymQDjDpHsw+hW+boY+m3sSBopP6pvZsMVT8l/owlkMWh0tuLdxe*)
(*9EU32uj7R/ZaYZHoaGa1TDhucHJQqRHEw05ZTwq9ciPIED8RznfqhaPJr+gO*)
(*DEdAsp5X7lATP34+KqldkeWDio+6D41cHcLyJHXJk8TZyEZTomB4sgdIitYz*)
(*+p8OgG7PUYee4kGs4wmgOcG2CAx9Fi9JyedRpnS+ZfFhOcZr5TNQJvpCW8bL*)
(*Wtl0JjhUT+W+1PMcigd5JaJOO2CEkkYj9agw2tgHGWcXmMKdhEkF1+pgsDEL*)
(*alLZiIBFTwu2trEH6ETskSfG6oeZZombbt6VoN5z8tzB3WK0WA/zJCkNhTdF*)
(*xDXiy8FApMAaa+f5AohKTmWRG5eBtY1ax6ZpF9AOy1z58acDNrZ29G6tduJ+*)
(*ae2M1Gw5st3JOafM7wJuh49ecXfOQd/U98oVx/Ohn6chOOBeERqMJw0NW5WC*)
(*y1ju6M1rH4Glq5Rxk3kRZSd1uKusBnHmxIXxC6JDkH7yhG8jRQfiGZnHNhYR*)
(*MMb5YfhlXiqkThL/IDTJhgdV1kQZXwNgKNkhyYHMGzzz+8Q3b/pCltYZj/U1*)
(*Z6xYmKKZELyHQ92RX4mZdOC+TscsG706UpwoHfV0joWFOG923pQY/LJ7N2Zw*)
(*OAkrlwZpv8EccNVk7d0znAfnH39CGAY+wPfXa1Repa3YN8OQNeE3iBUuND4G*)
(*PQg0h+uoO6cGoXSF8lOSWAJEfE3JCRVPR4tlsq6Vj7FQ8UftynZoNdhO6frt*)
(*MOYjN09SlO7rMmAxj9sT/XMfU//sOF1W7cB8Jym+haONEMtfpUgmWQsZohOs*)
(*znFd+FaeLUEluhz+E26tNPSvgpeXCI7kvXTA4vUusYmuG/D52QyfTL4AbJzK*)
(*/i40YIWK7XZ6xbuA3mpmbfFWpvCWhl7Qn8IK8fspktDofPjy9efNu68y4Nj0*)
(*qkHa9xy0ey4TycmVCmt0wiR379XhxXRX0amtRgT/HIeozg8gT8PMEzxcC3eY*)
(*jc7k6S+BX8m54Xtmw/DJ7VfNs/kXqM+rKveUpgw8k+MU5g6XwPNMg1U5sloU*)
(*IHyf7jrTDVqHi50b6eoxZMzHckO7CNR5AxSNGMew2FyusLezDT6f5fEl8CmH*)
(*wdhG2vLxCqTQ2HUu30jANm4vJtuqYjxz+CWl0xEPkHGieq+x4wsMYuw6z2ee*)
(*wmgeb/S4owkmFpso2JA548+fuWofXR7jFokd8Z2rN7H0r/X0ekkQDjmT8kt+*)
(*SYcpvbAuUp50+KMRyU8oNAFiwTHratXdMGPk9e6QfQeY93LwxlK9hy4qASK1*)
(*7Dk4keqqaczQiBdZZqiGRmewvNDKsO9hFqZm6NHIJ5ZCCdOZu4XkSbgj8etb*)
(*198mLHT+mPIuqhqYB9kJHs4NoWrDm9tHSJrBKu4WK+erNqB6+uE3T0MTemVl*)
(*STldjIIPOqK5ykZvccxSa5ZeOxl12gRNhK3yQdaLiC/IwwqzDNTOGkvaweJy*)
(*6TsqfxlkifaPuqohi860j/vnzH2Q6bqoPUesI+4t9VtvKD9Gyg7R3zJa5fDh*)
(*69ScaUghNHewGAUSZmJA+MD3bO5m6NB64urXN4KXny2QsnN0wQTH7ORPygk0*)
(*9AgMjvz+AV48aJU+b1yHZw06LpPsDIHpBeZzM88zUMzI7dahSVa0IkzymHPs*)
(*xdKQ8zqaRdlg/3E5yeXgHVCcjjhBUzAEcSx7A62xUzC0L/GKyLUSieRJSdoV*)
(*quCusvNFe9JYNNsvf38ouQzdxF9LijHkwOwPWl9a2S78Kvlk8cZZR3S4NMFT*)
(*PAFIsyHw0OLov7xGPGE0l+yFJMZPMmd17yGV6er7OoWLuKrEKWf6LgijLr+2*)
(*zvtUCpQbZ++z/fzHWfCbpG1yEzz2YSLOm2yBwbIwxWXHDLjbba6upNCL395c*)
(*Va560QaTYv/MyfECx3VnTCmVhtFyi9vHvHgW5iseU8utlaAhiazyXE0JXFG7*)
(*KUnCUweniy8ymBB0YVq5nbTDvxzzRTtpVKJuGId+UySm2bbjjAaesP/3pwj6*)
(*tTftSxfg6+u1yQQpDzQt5qX8TucH5EPXreSI2jFnsdDo2OY4iLGWH9KIVUPW*)
(*ebPZc3N3kZScRq3F2wMeGMhZXKCUw6/9ZBO8tLfxVZ7opIWKF0QPZb0tmZNF*)
(*xr7hyjCSQJio2yhKoqmDI6YnD3d+KsAN91atU2NvoDH4U4+cWxuMz7eYvxUZ*)
(*RCeulScyU6Vw5rjZHSLtZlBgaGEhXpzGgIFiRdEbNUA6+GmjwagHXB6tP5zW*)
(*zEXargBDe48qoNoJOPPSeAjtnvpzaU81491HbWwTdq2YSj7v0u5UgOWHfx32*)
(*9FiE5VQDhabxGpBf5D8Qv5cB550dRxWKm8Ex7Mv6j+IqtGu5ohq1lQLHNK8v*)
(*VpS4goW/In1/ozWQM/m13ynwQI1fIhI/6O3RhJ4sSumXHThMS9g6UWkCpeak*)
(*mM0lLfg6FWHHGtmN9LKCF2YSKtHkaog6x412+L0tQvUm5R1yPdPXvEVQj4dT*)
(*KR8PtU7BI8d5xmqdcbjBHnlZcKMReQeYHq8qD8KknqFcz0oHXuCJe/2TcBTP*)
(*FUrckbLphmNnbEMJ2Ypwzq9P07V6ClR3j9HL9tWCwJ3yuHnWKky6HZDnsVEP*)
(*XUqPHkiI9cKfwyzEFUIDYMhLy+5OPYTkn6Tz+PVbYMeQK+uybxxk/vb+rpfz*)
(*HHp0yK972Prg2ONxTjI7HeCi5+Tg/eaC4dHSIXaMzhjCPOPArGEJzakj+RPs*)
(*rlDEnH3rb7AlDMzxKBO2d4Nok9t0fEwPVgfFFBzLCsQvr6iddl6Nw9q49cDt*)
(*D//4XUj3438Ugk0an8UheQikbtRfEQ3owTtsPeF0M4i6hJJcjlm18F8vlwVr*)
(*RCYuacd+KCFoQ1nDS4pO98rQi2bhFVVPGiYWKWc59U9CyFCxh/D5dmDYrjoj*)
(*PdMOSTNC7m1cgxhsNV3NcHEQS/IH1DPftqB+zDvL/cfd4BtfTcQo2oF8JseC*)
(*RE5m43BOMtW9QCP4r5K3cTDmMe6ynWTu/8qEgsEfZSlK7EAtaei5sYY/xhXm*)
(*ih4+ZQRUi7QdzMM6YGVDG6jztxYzbvp83Fp6DezlRV4veAvgr9L1gJt8gyjF*)
(*I36+7fkkbFqXSE6bd+D9ex9FzNnbMSIiPvj4z0VQybt6u+x2CxyKMb5ZzFwB*)
(*RKq9L7SjWmDk9NeoxJh0fICGG1f9c6Gb+XpdKH8TXrrg+8ZZexKPTa7Ueez3*)
(*Ih7efMEymQgkFMNvTdp9gMwjN3KDsxZ232KxKucbTBb/mThwugq/i5yniGQa*)
(*hI8f3d6XtSQC02hUsciALUTrilKcL3XBGlY386m/Lph6biWo4aEN8s3951lT*)
(*x4wIxzgbR+TgHMmbPwGyzphps7/4dz8BCDzFnt7fyccbB90SVD5tqEaAQlfT*)
(*F5Gd+o/hBl8PuJlVaeYcjIKOdtd169J+0FNofrSa0Qhjb1RWSyLeIcubSPN3*)
(*ckNgd+pZCb9cLnj3KBTnEZehpBWRhEVMDPz8IBPu0D4Gi5/9Ok5pdePzFDKL*)
(*lOpo/L2yft5ivw3o9bx+RnPVY2AQ+bMU4Sx067QcEB3qgRe2nZd6GFqAPpNb*)
(*VlG6DleWi2x0+oNQVcaDru7MA2Alpku6Z+AGsdzIA4ZyKHQm2pvm4D4kD5NO*)
(*RVTdAppJkqRc/vvgfLNIp+KeLVCfPZYMlBkYXzrcblv9FEsHRbzFaSswVCnw*)
(*4vPzr9FSvY1hPKMT9dyKinO4xmA6h688maEHxP55mK8sFmZsj1WRHJ7BxPFj*)
(*dcclZuHczMEff9qmf733qVbFciPGHRMdXbPsR77Ub3Jz//hfdNomnjBNRQt+*)
(*Qu3elGr03q2b+ioeAxwCI0cy2UbR2kpiOYVnAoMb8k687miD+nWz2pGmPiAd*)
(*m3UfUR6B+OX5n9GZ8ZAntFM9T3QHvqWci+mRNUe65LT4+rrHoHLTV5mTQQdq*)
(*hgb5dfARfJz/AtnSduCowDgtfc8Coj/ZFxBKF2FnhazfL4sBFNSWf8p2vBaV*)
(*UtXqvM70gvqPX2ohw0PQHk7PUyY5ixIM7oWvKQvQ5qD7jDxPEzDqW1iUvGyB*)
(*Us4bNPUa//qRHbv4r6B0pCuJ/UCxFwB5DHQJJnsTUN+v+E5QpQq/HbxPJmV9*)
(*A7dUV+f+mg2BmlWbYtjcCPwcczn+KjAD5dx2ksct+0A6wmt3/FAgWqm6dz/x*)
(*LoJ7WoSVj8jroYCJT/bgvA9E7ZyRmLPXwdGiqGepB/fg3ngd5REuc8z/r4H1*)
(*G7cbpBXsVClxO2Ju/jj3n8a9Rnenz58n0vRQaZNgJLGyFPlCDfd7vcrhk/5y*)
(*o4lQKd5NUy7TJm+CoV4Dq5KlMdyVYZgP3OoFtYrlItaDfvw8/cMymq4bCt7y*)
(*GzEctEGGPL5Vf9KAP3x5R9Myc6EfFU3083KBeHj+QvfVHrzx1HLyS4Dzvzz7*)
(*xFiWsQF0buuacqbPYbxyLSOB3Dwy3zC6fUGkFypYDwdfTh4EKzkNOrvKKOQ0*)
(*m5l+v54Acspk6pkNiUi+vz3I5u8GF6e3sxaYDZEw+PtqOI0urCTHutp3ncdS*)
(*IudAZToXTGZ6dyWBww78/c/3eMgyAV1KDfMctz12RU3Ghko3o21gyJH5xRFI*)
(*HdEc1vV4j3QhKvRUVW/ABU3ulkdl4Jc0YTdDzjLUr5VxuKLfA1HM2+5jxp0o*)
(*60D+6blbEwZn2UQlN1Rj09NHp0112yBRl6LM170MOoeFqCcIu7BarylKJXQA*)
(*NYcKtPYetmKg/eyxn8HjcGDEJ3VkcRLyNx5vvLKchGwC0rgD6U54u79DP8cX*)
(*DfvW7AKLHyzxr0lYxdx8JtJ7u4yQj1ngbvavN3nN1iA1USx5VlIGQ5xeMZhf*)
(*UUDB9/rcB9HmcE3bk+v4dxG4/cvJ812sAQbmSDiGHzJDf9+jbllfJ6Cu8BqJ*)
(*fecgZqR0j5JGv0OivpmjC2Xt+ORGz2MK50rIJKvKUqqqBY2S20bPA9txhAul*)
(*TEpiwVqh04T0Syfcvb3ZsllXg3QaJ9ltC3qAnNFbZ11kCikFeo61r2XDcY9D*)
(*ptJXEBTcy1QMZpPwRV7wFMeVVjj8hPz+bHc8JIquU5dHTsDda4/a82U+AOtt*)
(*tot+ha6oyfykcmo1DL+mRR4Z5XPD3Itmwmf8zLFFY6g+qPAaBrPUWG7MqOGB*)
(*xvHfBTl64DmX8Xr8jzEKnzt/s1jcDLsoPYPeCkmgO3vSKzFmA1xrUc9bD2/F*)
(*xFWR20nUiFb2KS5rrhV4dkSu8frZPNgKz41zVbaEZ/qqCyTXcjHEpU6a62Q8*)
(*dAfxcGh4FQF3nAxpde8jXMkP9h1MKwebTjG54oJuPDAWCy2QLIDfSee+RZrO*)
(*g2nqZ7aUtVlo+IjtmrrzoPXQhesbQT/+bk9bTVhqQ+UA3q3c0X4MKPApWqts*)
(*x4o/ZobuFcKQQyr5LdvLBfhedIz6+CtjB2Gy4KdPfFDqMqAXua0CeZy+ZOYW*)
(*aiA9X9PlcU4Oom5RxV3nZMPmetK4RSMT4KvgQa4nx2G1h4qajMUIaK9veAv7*)
(*14KVkslYrP8bcNO/KW/JEwyqchdbIgniUT2LnP1yWBIk8nDPPzPUhzxTyUPx*)
(*mbFgRkEl29hdgU5LCuyWgiG4Mp9EQlach5dydjHqVQNeCZtKPlzXBdm5Mz2e*)
(*OQswr7ThuOS+AN4h1vL/HVmAxvufI522s/EVbtB9sanGU66joWI8Vbjy50OS*)
(*sWkLxuW673JEaoHAD0HpLnlFmFUeYfToc0fmQyMJ1kK6kDl+z/hAnxPnkvaC*)
(*052NsJe+dCM2Ug4qOn5zqZeawK0E4YY1WnkY4u40lqq4DL5FrT03L+mi2LT9*)
(*r5MkURjTkiDvuemGpHSjkt2//WGUkSm2e8AMs3MNrd6NqMNuXegngkuycMqv*)
(*9tRUWCrseR08SacNRJu22BVdhkRc/S80XRXDwFO47fDGv7ush5ap/rj3OH+e*)
(*0PpP0RA4enJSfVAdAmuq9QGRf71I/JqaHinDOPhb1Nmynp2AxXzhXwcuY7Ax*)
(*zlLllDIG/BeCTdgW2dCkOI5fcocZ/tYduvXwrSmU0k3KsRLoo0ZwVPg+ORnq*)
(*hD9sbyO7hdf2uCZd/wgDzhsQ0aVaow8lSzyrzGmc3TKISU9VRai4odwL1tAz*)
(*Qr181dASTplHCnFvxELsVIn6kdFgcIooj9aMDkSe96ePZ2VIYtJ8cG6ZThCE*)
(*Dcubsp2MBYXN8aXLVlHo/c3FjeiOD7xtOG5CSu2F7Ox89Cm748C8l0HmKr8E*)
(*ngmlA+4qDcC5JuM0b4fwjclsa6W1AW79JY6/QTUOfE5UO5H14wAWlC6JHWNw*)
(*pHkqQadiFP4He9U0Aw==*)
(*"]]]*)


(* ::Print:: *)
(*SequenceForm["The harmonic frequencies:\n", TableForm[{3197.419058980466, 3193.82296683702, 3183.975931556111, 3159.4047655572745`, 3116.1868565063514`, 3097.960114075379, 3087.874418753801, 1505.4682353229414`, 1499.3357867275754`, 1480.8471670990864`, 1475.7921491990953`, 1418.062214564708, 1405.047599078378, 1370.4190597781346`, 1303.7981627202223`, 1292.8004924189427`, 1242.588122981129, 1186.8543699911263`, 1125.8829568563744`, 1074.8266212362694`, 957.091181493338, 910.6502142577558, 882.5704116705181, 766.3570840766524, 562.6814132910889, 435.80483936838505`, 298.9235053907368, 237.73070515976806`, 158.03958625118773`, 81.55426306101235}]]*)


(* ::Print:: *)
(**)


(* ::Print:: *)
(**)


(* ::Print:: *)
(**)


(* ::Print:: *)
(**)


(* ::Print:: *)
(**)


(* ::Print:: *)
(**)


(* ::Input:: *)
(*fres[[1,1]] (*Gradients*)(*Should all be very close to 0*)*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
