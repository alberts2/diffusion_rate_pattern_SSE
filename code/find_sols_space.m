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


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(*Solution surface for lambdaA, qAB a,d qBA*)
(* Fixed values for muA and muB *)
muAVal = 1;
muBVal = 1;

(* RegionPlot3D with lambdaB replaced by its formula *)
rule1 = RegionPlot3D[
  (* Constraints *)
  muAVal > 1/3*(-3*qAB + 2*qBA) &&
  qBA > 3*qAB/2 &&
  (1/3*(3*muAVal + 3*qAB - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qAB + 2*qBA)) > 0,
  
  (* Plot variables *)
  {lambdaA, 0, 1}, {qAB, 0, 1}, {qBA, 0, 1},
(*  {lambdaA, 0, 5}, {qAB, 0, 5}, {qBA, 0, 5},*)
  
  PlotPoints -> 50,
  Mesh -> None,
  PlotStyle -> Opacity[0.3, Gray]
];

rule2 = RegionPlot3D[
  (* Constraints *)
  qBA < 3*qAB/2 &&
  muBVal > 1/2*(3*qAB - 2*qBA) &&
  (1/3*(3*muAVal + 3*qAB - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qAB + 2*qBA)) > 0,
  
  (* Plot variables *)
  {lambdaA, 0, 1}, {qAB, 0, 1}, {qBA, 0, 1},
 (* {lambdaA, 0, 5}, {qAB, 0, 5}, {qBA, 0, 5},*)
  
  PlotPoints -> 50,
  Mesh -> None,
  PlotStyle -> Opacity[0.3, Blue]
];

combinedPlot = Show[rule1, rule2, AxesLabel -> {"\[Lambda]A", "qAB", "qBA"},ViewPoint -> {0, -5, 1}]



(*Solution surface for lambdaA, qBA a,d muA*)
(* Fixed values for muA and muB *)
qABVal = 0.1;
muBVal = 0.1;

(* RegionPlot3D with lambdaB replaced by its formula *)
rule1 = RegionPlot3D[
  (* Constraints *)
  muA > 1/3*(-3*qABVal + 2*qBA) &&
  qBA > 3*qABVal/2 &&
  (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0,
  
  (* Plot variables *)
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  
  PlotPoints -> 50,
  Mesh -> None,
  PlotStyle -> Opacity[0.3, Gray]
];

rule2 = RegionPlot3D[
  (* Constraints *)
  qBA < 3*qABVal/2 &&
  muBVal > 1/2*(3*qABVal - 2*qBA) &&
  (1/3*(3*muA + 3*qABVal - 2*qBA) > 0) &&
  (1/2*(2*muBVal - 3*qABVal + 2*qBA)) > 0,
  
  (* Plot variables *)
  {lambdaA, 0, 1}, {qBA, 0, 1}, {muA, 0, 1},
  
  PlotPoints -> 50,
  Mesh -> None,
  PlotStyle -> Opacity[0.3, Blue]
];

combinedPlot = Show[rule1, rule2, AxesLabel -> {"\[Lambda]A", "qBA", "muA"},ViewPoint -> {0, -5, 1}]



(*Solution surface for lambdaB, qAB a,d muB*)
(* Fixed values for muA and muB *)
qBAVal = 0.1;
muAVal = 0.1;

(* RegionPlot3D with lambdaB replaced by its formula *)
rule1 = RegionPlot3D[
  (* Constraints *)
  muAVal > 1/3*(-3*qAB + 2*qBAVal) &&
  qBAVal > 3*qAB/2 &&
  (1/3*(3*muAVal + 3*qAB - 2*qBAVal) > 0) &&
  (1/2*(2*muB - 3*qAB + 2*qBAVal)) > 0,
  
  (* Plot variables *)
  {lambdaB, 0, 1}, {qAB, 0, 1}, {muB, 0, 1},
  
  PlotPoints -> 50,
  Mesh -> None,
  PlotStyle -> Opacity[0.3, Gray]
];

rule2 = RegionPlot3D[
  (* Constraints *)
  qBAVal < 3*qAB/2 &&
  muB > 1/2*(3*qAB - 2*qBAVal) &&
  (1/3*(3*muAVal + 3*qAB - 2*qBAVal) > 0) &&
  (1/2*(2*muB - 3*qAB + 2*qBAVal)) > 0,
  
  (* Plot variables *)
  {lambdaB, 0, 1}, {qAB, 0, 1}, {muB, 0, 1},
  
  PlotPoints -> 50,
  Mesh -> None,
  PlotStyle -> Opacity[0.3, Blue]
];

combinedPlot = Show[rule1, rule2, AxesLabel -> {"\[Lambda]B", "qAB", "muB"},ViewPoint -> {0, -5, 1}]




