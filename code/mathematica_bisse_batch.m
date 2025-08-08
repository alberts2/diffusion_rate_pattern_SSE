(* ::Package:: *)

(*Finding rate relationships under BiSSE model where rates must be positive and 
  net diversifation rates are non-zero (for both states)
*)

(*Import random matrix from R*)
M = Import["/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/patterns_statio_batch.csv"];
(*Rationalize the matrix elements in order to use exact solution from Solve*)
M = Rationalize[M];
(*M = {{0.4,0.6}}*)

(*Define the functions*)
PiA[pia_,pib_] := (la-ma-qab)*pia + qba*pib;
PiB[pia_,pib_] := (lb-mb-qba)*pib + qab*pia;

(*Define the total rate functions*)
PhiA := (la+qba)-(ma+qab);
PhiB := (lb+qab)-(mb+qba);

(*Number of rows*)
nrow := Length[M];

(*store the solutions*)
solsset = Array[ans,nrow];
	Do[If[M[[i,1]]==M[[i,2]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,2]]]==0 && PiB[M[[i,1]],M[[i,2]]]==0 && PhiA==PhiB && la!=ma && lb!=mb && la >0 && lb>0 && ma>0 && mb>0 && qab>0 && qba>0},{la,lb,ma,mb,qab,qba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]>M[[i,2]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,2]]]==0 && PiB[M[[i,1]],M[[i,2]]]==0 && PhiA>PhiB && la!=ma && lb!=mb && la >0 && lb>0 && ma>0 && mb>0 && qab>0 && qba>0},{la,lb,ma,mb,qab,qba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]<M[[i,2]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,2]]]==0 && PiB[M[[i,1]],M[[i,2]]]==0 && PhiA<PhiB && la!=ma && lb!=mb && la >0 && lb>0 && ma>0 && mb>0 && qab>0 && qba>0},{la,lb,ma,mb,qab,qba},MaxExtraConditions->All]]
	],{i,nrow}]
(*Save the output as .txt*)
SetDirectory["/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/"];
(*Print[solsset[[1]]]*)
Do[Export["output_sim_"<>ToString[i]<>".txt",solsset[[i]]],{i,1,nrow}];


