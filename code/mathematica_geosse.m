(* ::Package:: *)

(*Import random matrix from R*)
M = Import["/Users/albertsoewongsono/Documents/Code Testing/SSA/diffusion_Rpackage/patterns_statio.csv"];
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

solsset = {sols1};
	If[M[[1,1]]==M[[1,2]] && M[[1,2]]==M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA==PhiB && PhiB==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]>M[[1,2]] && M[[1,2]]>M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA>PhiB && PhiB>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]<M[[1,2]] && M[[1,2]]<M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA<PhiB && PhiB<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]>M[[1,2]] && M[[1,2]]==M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA>PhiB && PhiB==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]<M[[1,2]] && M[[1,2]]==M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA<PhiB && PhiB==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]==M[[1,2]] && M[[1,2]]<M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA==PhiB && PhiB<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]==M[[1,2]] && M[[1,2]]>M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA==PhiB && PhiB>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]>M[[1,2]] && M[[1,2]]<M[[1,3]] && M[[1,1]]==M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA>PhiB && PhiB<PhiAB && PhiA==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]>M[[1,2]] && M[[1,2]]<M[[1,3]] && M[[1,1]]<M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA>PhiB && PhiB<PhiAB && PhiA<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]>M[[1,2]] && M[[1,2]]<M[[1,3]] && M[[1,1]]>M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA>PhiB && PhiB<PhiAB && PhiA>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]<M[[1,2]] && M[[1,2]]>M[[1,3]] && M[[1,1]]==M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA<PhiB && PhiB>PhiAB && PhiA==PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]<M[[1,2]] && M[[1,2]]>M[[1,3]] && M[[1,1]]<M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA<PhiB && PhiB>PhiAB && PhiA<PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
	If[M[[1,1]]<M[[1,2]] && M[[1,2]]>M[[1,3]] && M[[1,1]]>M[[1,3]],
	solsset[[1]] = Quiet[Solve[{PiA[M[[1,1]],M[[1,3]]]==0 && PiB[M[[1,2]],M[[1,3]]]==0 && PiAB[M[[1,1]],M[[1,2]],M[[1,3]]]==0 && PhiA<PhiB && PhiB>PhiAB && PhiA>PhiAB && wa>0 && wb>0 && ea>0 && eb>0 && dab>0 && dba>0 && bab>0},{wa,wb,ea,eb,bab,dab,dba},MaxExtraConditions->All]]
	];
(*Save the output as .txt*)
SetDirectory["/Users/albertsoewongsono/Documents/Code Testing/SSA/diffusion_Rpackage/"];
Export["output_sim"<>".txt",solsset[[1]]]


