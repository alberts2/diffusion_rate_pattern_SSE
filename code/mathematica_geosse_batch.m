(* ::Package:: *)

(*Import random matrix from R*)
M = Import["/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/patterns_statio_batch.csv"];
(*Rationalize the matrix elements in order to use exact solution from Solve*)
M = Rationalize[M];
(*Define the functions*)
PiA[pia_,piab_] := wa*(pia+piab) + piab*bab + eb*piab - dab*pia - ea*pia;
PiB[pib_,piab_] := wb*(pib+piab) + piab*bab + ea*piab - dba*pib - eb*pib;
PiAB[pia_,pib_,piab_] := dab*pia + dba*pib - piab*(ea+eb+bab);
(*Define the total rate functions*)
PhiA := 2*wa + bab + eb - dab - ea;
PhiB := 2*wb + bab + ea - dba - eb;
PhiAB := dab + dba - ea - eb - bab;
(*Number of rows*)
nrow := Length[M];
(*store the solutions*)
solsset = Array[ans,nrow];
Do[If[M[[i,1]]==M[[i,2]] && M[[i,2]]==M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA==PhiB && PhiB==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]>M[[i,2]] && M[[i,2]]>M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA>PhiB && PhiB>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]<M[[i,2]] && M[[i,2]]<M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA<PhiB && PhiB<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]>M[[i,2]] && M[[i,2]]==M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA>PhiB && PhiB==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]<M[[i,2]] && M[[i,2]]==M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA<PhiB && PhiB==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]==M[[i,2]] && M[[i,2]]<M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA==PhiB && PhiB<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]==M[[i,2]] && M[[i,2]]>M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA==PhiB && PhiB>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,1]]==M[[i,3]] && M[[i,3]]>M[[i,2]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA>PhiB && PhiB<PhiAB && PhiA==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,2]]<M[[i,1]] && M[[i,1]]<M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA>PhiB && PhiB<PhiAB && PhiA<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,2]]<M[[i,3]] && M[[i,3]]<M[[i,1]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA>PhiB && PhiB<PhiAB && PhiA>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,3]]<M[[i,2]] && M[[i,1]]==M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA<PhiB && PhiB>PhiAB && PhiA==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,2]]>M[[i,3]] && M[[i,3]]>M[[i,1]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA<PhiB && PhiB>PhiAB && PhiA<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[i,2]]>M[[i,1]] && M[[i,1]]>M[[i,3]],
	solsset[[i]] = Quiet[Solve[{PiA[M[[i,1]],M[[i,3]]]==0 && PiB[M[[i,2]],M[[i,3]]]==0 && PiAB[M[[i,1]],M[[i,2]],M[[i,3]]]==0 && PhiA<PhiB && PhiB>PhiAB && PhiA>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	],{i,nrow}]
(*Save the output as .txt*)
SetDirectory["/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/data/"];
Do[Export["output_sim_"<>ToString[i]<>".txt",solsset[[i]]],{i,1,nrow}];
