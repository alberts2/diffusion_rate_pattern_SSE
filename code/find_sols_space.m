(* ::Package:: *)

(*Note that maximum/minimum may be attained at the boundary of the parameter, so if you use open interval, mathematica
will return a warning message. Can change all < and > to <= and >= to check if the output is still consistent*)

(*For each parameter, I give a maximum value so mathematica won't return infinity*)

(*Ruleset 1*)
constraints = {lambdaA == 1/3*(3*muA + 3*qAB - 2*qBA), lambdaB == 1/2*(2*muB - 3*qAB + 2*qBA),
				muA > 1/3*(-3*qAB + 2*qBA), muB > 0, qAB > 0, qBA > 3*qAB/2,
				lambdaA > 0, lambdaB > 0, muA > 0, muB > 0, qAB > 0, qBA > 0,
				lambdaA <= 10, lambdaB <= 10, muA <= 10, muB <= 10, qAB <= 10, qBA <= 10};

(*maxLambdaA = Maximize[{lambdaA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minLambdaA = Minimize[{lambdaA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxLambdaB = Maximize[{lambdaB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minLambdaB = Minimize[{lambdaB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxMuA = Maximize[{muA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minMuA = Minimize[{muA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxMuB = Maximize[{muB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minMuB = Minimize[{muB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxQAB = Maximize[{qAB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minQAB = Minimize[{qAB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxQBA = Maximize[{qBA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minQBA = Minimize[{qBA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
*)



(* ::Input:: *)
(**)


(*Note that maximum/minimum may be attained at the boundary of the parameter, so if you use open interval, mathematica
will return a warning message. Can change all < and > to <= and >= to check if the output is still consistent*)

(*For each parameter, I give a maximum value so mathematica won't return infinity*)

(*Ruleset 2*)
constraints = {lambdaA == 1/3*(3*muA + 3*qAB - 2*qBA), lambdaB == 1/2*(2*muB - 3*qAB + 2*qBA),
				0< qBA < 3*qAB/2, muB > 1/2*(3*qAB - 2*qBA),
				lambdaA > 0, lambdaB > 0, muA > 0, muB > 0, qAB > 0, qBA > 0,
				lambdaA <= 10, lambdaB <= 10, muA <= 10, muB <= 10, qAB <= 10, qBA <= 10};

maxLambdaA = Maximize[{lambdaA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minLambdaA = Minimize[{lambdaA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxLambdaB = Maximize[{lambdaB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minLambdaB = Minimize[{lambdaB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxMuA = Maximize[{muA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minMuA = Minimize[{muA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxMuB = Maximize[{muB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minMuB = Minimize[{muB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxQAB = Maximize[{qAB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minQAB = Minimize[{qAB, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]

maxQBA = Maximize[{qBA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]
minQBA = Minimize[{qBA, constraints}, {lambdaA, lambdaB, muA, muB, qAB, qBA}]


(* ::Input:: *)
(**)


(*stationary = [0.6,0.4]*)
(*Solution surface for lambdaA, qBA and muA*)

(*Fixed values for qAB, muB, and lambdaB*)

(*this scenario satisfies ruleset 1 but not ruleset 2*)
(*one example will gives us under ruleset 1 qBA = 0.2, muA = 0.1, lambdaA =0.2/3*)
qABVal = 0.1;
muBVal = 0.1;
lambdaBVal = 0.2;

(*this scenario satisfies ruleset 2 but not ruleset 1*)
(*one example will give us under ruleset 2 qBA=0.1, muA = 0.1 and lambdaA = 0.4/3*)
(*qABVal = 0.1;
muBVal = 0.2;
lambdaBVal = 0.1 ;*)
(*tolerance for the lambdaA != muA (and for B)*)
tol = 0.00001;
(*Note LambdaB cannot be equal to MuA. Similarly for LambdaA and MuA.*)

(*Ruleset 1*)

(*plot the region bounded by the inequalities *)
region1 = RegionPlot3D[
  (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0 &&
   muA > 1/3*(-3*qABVal + 2*qBA) &&
   qBA > 3*qABVal/2 &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for lambdaA equality *)
surf1 = ContourPlot3D[
  lambdaA - 1/3*(3*muA + 3*qABVal - 2*qBA) == 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
    muA > 1/3*(-3*qABVal + 2*qBA) &&
    qBA > 3*qABVal/2 &&
    (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
    (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0 &&
     Abs[lambdaA - muA] > tol &&
     Abs[lambdaBVal - muBVal] > tol
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for lambdaB equality *)
surf2 = ContourPlot3D[
  lambdaBVal - 1/2*(2*muBVal - 3*qABVal + 2*qBA) == 0,
   {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
    muA > 1/3*(-3*qABVal + 2*qBA) &&
    qBA > 3*qABVal/2 &&
    (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
    (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0 &&
     Abs[lambdaA - muA] > tol &&
     Abs[lambdaBVal - muBVal] > tol
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

(*Ruleset 2*)

(*plot the region bounded by the inequalities *)
region2 = RegionPlot3D[
  (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0 &&
   qBA < 3*qABVal/2 &&
   muBVal > 1/2*(3*qABVal - 2*qBA)&&
   Abs[lambdaA - muA] > tol,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  PlotStyle -> Opacity[0.2, Blue],
  Mesh -> None
];

(* plot the plan for lambdaA equality *)
surf3 = ContourPlot3D[
  lambdaA - 1/3*(3*muA + 3*qABVal - 2*qBA) == 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
   (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0 &&
   qBA < 3*qABVal/2 &&
   muBVal > 1/2*(3*qABVal - 2*qBA)&&
   Abs[lambdaA - muA] > tol
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for lambdaB equality *)
surf4 = ContourPlot3D[
  lambdaBVal - 1/2*(2*muBVal - 3*qABVal + 2*qBA) == 0,
   {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
   RegionFunction -> Function[{lambdaA, qBA, muA},
   (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
   (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0 &&
   qBA < 3*qABVal/2 &&
   muBVal > 1/2*(3*qABVal - 2*qBA)&&
   Abs[lambdaA - muA] > tol
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

Show[region1,surf1,surf2,AxesLabel -> {"\[Lambda]A", "qBA", "muA"}]
Show[region2,surf3,surf4,AxesLabel -> {"\[Lambda]A", "qBA", "muA"}]



(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(*stationary = [0.5,0.5]*)
(*Solution surface for lambdaA, qBA and muA*)

(*Fixed values for qAB, muB, and lambdaB*)

(*this scenario satisfies ruleset 1 but not ruleset 2*)
(*one example will gives us under ruleset 1 qBA = 0.3, muA = 0.2, lambdaA = 0.1 *)
(*In fact any muA > -qAB + qBA and lambdaA given by the equality will do. For example:
qBA = 0.3, muA = 0.3, lambdaA= 0.2 *)
(*In any case, both correspond to the same net diversification in A, so for inferences,
perhaps it's better to parameterize events using net diversification and turnover (see my slack)*

can design an experiment where we do inference with likelihood function written using actual rates
and using net diversification and turnover instead*)

(*qABVal = 0.2;
muBVal = 0.1;
lambdaBVal = 0.2;
*)
(*this scenario satisfies ruleset 2 but not ruleset 1*)
(*one example will give us under ruleset 2 qBA=0.1, muA = 0.2 and lambdaA = 0.3*)
(*In fact, any muA > 0 will do, for example: muA = 0.1, which will give us lambda A = 0.2 with qBA = 0.1*)

qABVal = 0.2;
muBVal = 0.2;
lambdaBVal = 0.1 ;
(*tolerance for the lambdaA != muA (and for B)*)
tol = 0.00001;
(*Note LambdaB cannot be equal to MuA. Similarly for LambdaA and MuA.*)

(*Ruleset 1*)

(*plot the region bounded by the inequalities *)
region1 = RegionPlot3D[
  muA + qABVal - qBA > 0 &&
  muBVal - qABVal + qBA > 0 &&
  muA > -qABVal + qBA &&
   qBA > qABVal &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for lambdaA equality *)
surf1 = ContourPlot3D[
  lambdaA - muA - qABVal + qBA == 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
    muA + qABVal - qBA > 0 &&
    muBVal - qABVal + qBA > 0 &&
    muA > -qABVal + qBA &&
    qBA > qABVal &&
    Abs[lambdaA - muA] > tol &&
    Abs[lambdaBVal - muBVal] > tol
    lambdaA > 0 &&
    qBA > 0 &&
    muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for lambdaB equality *)
surf2 = ContourPlot3D[
  lambdaBVal - muBVal + qABVal - qBA == 0,
   {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
    muA + qABVal - qBA > 0 &&
    muBVal - qABVal + qBA > 0 &&
    muA > -qABVal + qBA &&
    qBA > qABVal &&
    Abs[lambdaA - muA] > tol &&
    Abs[lambdaBVal - muBVal] > tol&&
    lambdaA > 0 &&
   qBA > 0 &&
   muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

(*Ruleset 2*)

(*plot the region bounded by the inequalities *)
region2 = RegionPlot3D[
  muA + qABVal - qBA > 0 &&
  muBVal - qABVal + qBA > 0 &&
  qBA < qABVal &&
  muBVal > qABVal - qBA &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for lambdaA equality *)
surf3 = ContourPlot3D[
  lambdaA - muA - qABVal + qBA == 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
  muA + qABVal - qBA > 0 &&
  muBVal - qABVal + qBA > 0 &&
  qBA < qABVal &&
  muBVal > qABVal - qBA &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for lambdaB equality *)
surf4 = ContourPlot3D[
  lambdaBVal - muBVal + qABVal - qBA == 0,
   {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
   RegionFunction -> Function[{lambdaA, qBA, muA},
  muA + qABVal - qBA > 0 &&
  muBVal - qABVal + qBA > 0 &&
  qBA < qABVal &&
  muBVal > qABVal - qBA &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

Show[region1,surf1,surf2,AxesLabel -> {"\[Lambda]A", "qBA", "muA"}]
Show[region2,surf3,surf4,AxesLabel -> {"\[Lambda]A", "qBA", "muA"}]


(*stationary = [0.4,0.6]*)
(*Solution surface for lambdaA, qBA and muA*)

(*Fixed values for qAB, muB, and lambdaB*)

(*this scenario satisfies ruleset 1 but not ruleset 2*)
(*one example will gives us under ruleset 1 qBA = 0.7/3, muA = 0.2, lambdaA = 0.1/2 *)
(*In fact, any values where mu_A > 1/2*(-2qAB + 3qBA) and \lambda_A = 1/2(2mu_A+2qAB-3qBA) 
will lead to the same stationary frequencies.*)
(*For example: qBA = 0.7/3, muA = 0.4, lambdaA = 0.25*)
(*Again, the net diversification is A is preserved in these two examples*)

(*qABVal = 0.2;
muBVal = 0.1;
lambdaBVal = 0.2;*)


(*this scenario satisfies ruleset 2 but not ruleset 1*)
(*one example will give us under ruleset 2 qBA=0.1/3, muA = 0.1 and lambdaA = 0.25*)
(*In fact, any values mu_A > 0 and \lambda_A = 1/2(2mu_A+2qAB-3qBA) will do*)
(*For example: qBA=0.1/3, muA = 0.3 and lambdaA = 0.45 )
(*Again, the net diversification is A is preserved in these two examples*)
*)

qABVal = 0.2;
muBVal = 0.2;
lambdaBVal = 0.1;

(*the polynomials in scenario 1 and 2 differ but they still lead to the same frequency*)

(*tolerance for the lambdaA != muA (and for B)*)
tol = 0.00001;
(*Note LambdaB cannot be equal to MuA. Similarly for LambdaA and MuA.*)

(*Ruleset 1*)

(*plot the region bounded by the inequalities *)
region1 = RegionPlot3D[
  1/2*(2*muA + 2*qABVal - 3*qBA) > 0 &&
  1/3*(3*muBVal - 2*qABVal + 3*qBA) > 0 &&
  (**)
   muA > 1/2*(-2*qABVal + 3*qBA) &&
   qBA > 2*qABVal/3 &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for lambdaA equality *)
surf1 = ContourPlot3D[
  lambdaA - 1/2*(2*muA + 2*qABVal - 3*qBA) == 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
    1/2*(2*muA + 2*qABVal - 3*qBA) > 0 &&
    1/3*(3*muBVal - 2*qABVal + 3*qBA) > 0 &&
    muA > 1/2*(-2*qABVal + 3*qBA) &&
    qBA > 2*qABVal/3 &&
    Abs[lambdaA - muA] > tol &&
    Abs[lambdaBVal - muBVal] > tol &&
    lambdaA > 0 &&
    qBA > 0 &&
    muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for lambdaB equality *)
surf2 = ContourPlot3D[
  lambdaBVal - 1/3*(3*muBVal - 2*qABVal + 3*qBA) == 0,
   {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
    1/2*(2*muA + 2*qABVal - 3*qBA) > 0 &&
    1/3*(3*muBVal - 2*qABVal + 3*qBA) > 0 &&
    muA > 1/2*(-2*qABVal + 3*qBA) &&
    qBA > 2*qABVal/3 &&
    Abs[lambdaA - muA] > tol &&
    Abs[lambdaBVal - muBVal] > tol &&
    lambdaA > 0 &&
   qBA > 0 &&
   muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

(*Ruleset 2*)

(*plot the region bounded by the inequalities *)
region2 = RegionPlot3D[
  1/2*(2*muA + 2*qABVal - 3*qBA) > 0 &&
  1/3*(3*muBVal - 2*qABVal + 3*qBA) > 0 &&
  qBA < 2*qABVal/3 &&
  muBVal > 1/3*(2*qABVal - 3*qBA) &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for lambdaA equality *)
surf3 = ContourPlot3D[
  lambdaA - 1/2*(2*muA + 2*qABVal - 3*qBA) == 0,
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  RegionFunction -> Function[{lambdaA, qBA, muA},
  1/2*(2*muA + 2*qABVal - 3*qBA) > 0 &&
  1/3*(3*muBVal - 2*qABVal + 3*qBA) > 0 &&
  qBA < 2*qABVal/3 &&
  muBVal > 1/3*(2*qABVal - 3*qBA) &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for lambdaB equality *)
surf4 = ContourPlot3D[
  lambdaBVal - 1/3*(3*muBVal - 2*qABVal + 3*qBA) == 0,
   {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
   RegionFunction -> Function[{lambdaA, qBA, muA},
  1/2*(2*muA + 2*qABVal - 3*qBA) > 0 &&
  1/3*(3*muBVal - 2*qABVal + 3*qBA) > 0 &&
  qBA < 2*qABVal/3 &&
  muBVal > 1/3*(2*qABVal - 3*qBA) &&
   Abs[lambdaA - muA] > tol &&
   Abs[lambdaBVal - muBVal] > tol &&
   lambdaA > 0 &&
   qBA > 0 &&
   muA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

Show[region1,surf1,surf2,AxesLabel -> {"\[Lambda]A", "qBA", "muA"}]
Show[region2,surf3,surf4,AxesLabel -> {"\[Lambda]A", "qBA", "muA"}]


(*stationary = [piA = 0.2, piB = 0.1, piAB = 0.7]*)
(*Solution surface for wA, dBA and eA*)

(*Fixed values for dAB, wB, eB and bAB*)

(*this scenario satisfies ruleset 1 but not ruleset 2*)
(*one example will gives us under ruleset 1 qBA = 0.7/3, muA = 0.2, lambdaA = 0.1/2 *)
(*In fact, any values where mu_A > 1/2*(-2qAB + 3qBA) and \lambda_A = 1/2(2mu_A+2qAB-3qBA) 
will lead to the same stationary frequencies.*)
(*For example: qBA = 0.7/3, muA = 0.4, lambdaA = 0.25*)
(*Again, the net diversification is A is preserved in these two examples*)

dABVal = 0.1;
eBVal = 0.1;
wBVal = 0.075;
bABVal = 0.01;


(*Ruleset 1*)

(*plot the region bounded by the inequalities *)
region1 = RegionPlot3D[
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal <= 49*eBVal/15 &&
   dBA > 1/2*(-18*dABVal + 63*eBVal) &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0,
  {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for wA equality *)
surf1 = ContourPlot3D[
  wA - 1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) == 0,
  {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  RegionFunction -> Function[{wA, dBA, eA},
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal <= 49*eBVal/15 &&
   dBA > 1/2*(-18*dABVal + 63*eBVal) &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0
  ],
  PlotPoints -> 60, MaxRecursion -> 3,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for wB equality *)
surf2 = ContourPlot3D[
  wBVal - 1/8*(-2*dABVal + 8*eBVal) == 0, (*this will be a fixed value since wB, dAB, and eB are fixed, so this will not be plotted*)
  {wA, 0, 1}, {dBA, 0, 3}, {eA, 0, 1},
  RegionFunction -> Function[{wA, dBA, eA},
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal <= 49*eBVal/15 &&
   dBA > 1/2*(-18*dABVal + 63*eBVal) &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0
  ],
  PlotPoints -> 60, MaxRecursion -> 3,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

(* plot the plan for eA equality *)
surf3 = ContourPlot3D[
  eA - 1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) == 0,
  {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  RegionFunction -> Function[{wA, dBA, eA},
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal <= 49*eBVal/15 &&
   dBA > 1/2*(-18*dABVal + 63*eBVal) &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Blue]
];

(*Ruleset 2*)

(*plot the region bounded by the inequalities *)
region2 = RegionPlot3D[
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal < 4*eBVal &&
   dABVal > 49*eBVal/15 &&
   dBA > 9*dABVal/14 &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0,
  {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  PlotStyle -> Opacity[0.2, Gray],
  Mesh -> None
];

(* plot the plan for wA equality *)
surf4 = ContourPlot3D[
  wA - 1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) == 0,
  {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  RegionFunction -> Function[{wA, dBA, eA},
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal < 4*eBVal &&
   dABVal > 49*eBVal/15 &&
   dBA > 9*dABVal/14 &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Red]
];

(* plot the plan for wB equality *)
surf5 = ContourPlot3D[
  wBVal - 1/8*(-2*dABVal + 8*eBVal) == 0,
   {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  RegionFunction -> Function[{wA, dBA, eA},
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal < 4*eBVal &&
   dABVal > 49*eBVal/15 &&
   dBA > 9*dABVal/14 &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Yellow]
];

(* plot the plan for eA equality *)
surf6 = ContourPlot3D[
  eA - 1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) == 0,
   {wA, 0, 3}, {dBA, 0, 3}, {eA, 0, 3},
  RegionFunction -> Function[{wA, dBA, eA},
  1/9*(-7*bABVal + 2*dABVal + 2/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal)- 7*eBVal) > 0&&
  1/8*(-2*dABVal + 8*eBVal) > 0 &&
  1/7*(-7*bABVal + 2*dABVal + dBA - 7*eBVal) > 0 &&
  (**)
   bABVal < 1/63*(18*dABVal + 2*dBA - 63*eBVal) &&
   dABVal < 4*eBVal &&
   dABVal > 49*eBVal/15 &&
   dBA > 9*dABVal/14 &&
   wA > 0 &&
   dBA > 0 &&
   eA > 0
  ],
  PlotPoints -> 30, MaxRecursion -> 2,
   Mesh -> None,
  ContourStyle -> Opacity[1, Blue]
];

Show[surf1,AxesLabel -> {"wA", "dBA", "eA"}]
Show[surf2,AxesLabel -> {"wA", "dBA", "eA"}]
Show[surf3,AxesLabel -> {"wA", "dBA", "eA"}]
(*Show[region1,surf1,surf2,surf3,AxesLabel -> {"wA", "dBA", "eA"}]
Show[region2,surf4,surf5,surf6,AxesLabel -> {"wA", "dBA", "eA"}]*)



