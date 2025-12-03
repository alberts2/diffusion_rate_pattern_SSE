(* ::Package:: *)

(*Finding rate relationships under BiSSE model where rates must be positive and 
  net diversifation rates are non-zero (for both states)
*)

(*Import random matrix from R*)
(*M = {{0.6,0.4}};*)
(*M = {{0.5,0.5}};*)
(*M = {{0.4,0.6}};*)
(*M = {{0.39,0.61}}; *)
(*this has the similar rulesets as to above and below, only differ in coefficients*)
(*M = {{0.3999,0.6001}};*)
(*M = {{0.4,0.8}};*)
M = {{0.1,0.9}};
M = Rationalize[M];

(*Define the functions*)
PiA[pia_,pib_] := (la-ma-qab)*pia + qba*pib;
PiB[pia_,pib_] := (lb-mb-qba)*pib + qab*pia;

(*Define the total rate functions*)
PhiA := (la+qba)-(ma+qab);
PhiB := (lb+qab)-(mb+qba);

(*store the solutions*)
solsset = {sols1};
	If[M[[1,1]]==M[[1,2]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,2]]]==0 && PiB[M[[1,1]],M[[1,2]]]==0 && PhiA==PhiB && la!=ma && lb!=mb && la >0 && lb>0 && ma>0 && mb>0 && qab>0 && qba>0},{la,lb,ma,mb,qab,qba},MaxExtraConditions->All]]
	]
	If[M[[1,1]]>M[[1,2]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,2]]]==0 && PiB[M[[1,1]],M[[1,2]]]==0 && PhiA>PhiB && la!=ma && lb!=mb && la >0 && lb>0 && ma>0 && mb>0 && qab>0 && qba>0},{la,lb,ma,mb,qab,qba},MaxExtraConditions->All]]
	]
	If[M[[1,1]]<M[[1,2]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,2]]]==0 && PiB[M[[1,1]],M[[1,2]]]==0 && PhiA<PhiB && la!=ma && lb!=mb && la >0 && lb>0 && ma>0 && mb>0 && qab>0 && qba>0},{la,lb,ma,mb,qab,qba},MaxExtraConditions->All]]
	]
