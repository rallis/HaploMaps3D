(* ::Package:: *)

BeginPackage["GenResearch`"]


linearCGRgrayscale::usage="linearCGRgrayscale[seq_String, {a1_, a2_, a3_, a4_}, k_Integer?Positive] returns FCGR of a DNA sequence using 2^k resolution. assumes linear, non-circular dna string"
convertFCGRtoCGR::usage="convertFCGRtoCGR[x_] converts the FCGR matrix to the CGR matrix by making all elements greater than 0, equal to 1"
visualizeFCGR::usage="visualizeFCGR[fcgrplot_,colors_,colorbounds_,imgsize_] colours FCGR plot"
SSIMmax::usage="SSIMmax[img1_,img2_,max_] = non compiled version of SSIMmax"
HamDist::usage="HamDist[x_List, y_List] = Hamming Distance of two vectors"
mydescriptor::usage="mydescriptor[img_,winsizes_,steps_,histbins_] = Returns the descriptor of a  matrix/image, given the windowsizes, windowsteps and histogram bins"
mymail::usage="mymail sends an email"
Kruskal::usage="kruskal code from wolfram project"
SSIMcompiled::usage = "SSIMcompiled[img1, img2] gives the Structural Similarity distance between two images"
ApproxInfoDist::usage = "ApproxInfoDist[img1, img2] gives the Approximate Information distance between two images"
mydescriptorVer2::usage = "Descriptor[img, winsizes, steps, histbins] gives the descriptor of a matrix, given the windowsizes, windowsteps and histogram bins"
mds::usage="mds[delta_List,dimensions_,dbg_] gets as input the distances delta_{ij}, and performs MDS in 2D"
mdsVer2::usage = 
	"mdsVer2[delta, dim] performs MDS in 2 or 3 dimensions on the given distance matrix.\nMDS[delta, dim, accuracy] does the same, with a given argument to N."
getUpperHalf::usage="getUpperHalf[matr_List] accepts as input an nxn matrix and returns the upper half of it (excluding the diagonal), flattened"
linearRegr::usage="linearRegr[x_List,y_List] gets as input one list X and one list Y and returns the coefficients of linear regression in these two sets"
MDSError::usage="MDSError[deltaMatrix_List,dMatrix_List] gets d and delta matrices and prints the kruskal error and (optionally) the raw error too"
swapMatrix::usage="swapMatrix[matr_List,i_,j_] gets a matrix and two indices and swap these columns and rows in that matrix"
stressPerPoint::usage="stressPerPoint[delta_,d_] returns the stress per point given delta (input) distances and d (output) distances"
noFmdsError::usage="noFmdsError[deltaMatrix_List,dMatrix_List] returns the error, without functioni f (coefs={0,1})"


Begin["`Private`"]



linearCGRgrayscale[seq_String, {a1_, a2_, a3_, a4_}, k_Integer?Positive]:=
	Module[{shifts, pts, arrayrules, resimg},
		shifts = StringCases[seq,a1|a2|a3|a4] /. {a1->{0, 0}, a2->{0, 2^k}, a3->{2^k, 2^k}, a4->{2^k, 0}};
		pts = FoldList[IntegerPart[(#+#2)/2]&, 2^(k-1), #]& /@ Transpose[shifts];
		pts = Drop[#, k]&/@pts;
		arrayrules = Rule @@@ Tally[Transpose[{2^k-pts[[2]], 1+pts[[1]]}]];
		resimg = SparseArray[arrayrules, {2^k, 2^k}, 0];
		Return[Normal[resimg]];
	];


convertFCGRtoCGR[x_]:=If[x>0,Return[1],Return[0]];
SetAttributes[convertFCGRtoCGR,Listable];


visualizeFCGR[fcgrplot_,colors_,colorbounds_,imgsize_]:=
	Module[{clrules={},i,j},
		If[Length@colors!=Length@colorbounds,Print["Error in dimensions!"];Return[]];
		For[i=1,i<=Length@colorbounds-1,i++,
			For[j=colorbounds[[i]],j<colorbounds[[i+1]],j++,
				AppendTo[clrules,{j->colors[[i]]}];
		]];
		AppendTo[clrules,{_->colors[[i]]}];
		clrules=Flatten@clrules;
		(*Print[{Min[fcgrplot],Max[fcgrplot]}];*)
		Return[ArrayPlot[fcgrplot,ColorRules->clrules,ImageSize->imgsize]];
]


SSIMmax[img1_,img2_,max_]:=
	Module[{w, c1, c2, m1, m2, m1sq, m2sq, m1m2, sigma1sq, sigma2sq,sigma12, ssimmap, mssim},
		w = GaussianMatrix[{5, 1.5}]; c1 = (0.01*max)^2; c2 = (0.03*max)^2;
		m1 = ListCorrelate[w, img1]; 
		m2 = ListCorrelate[w, img2];
		m1sq = m1*m1; m2sq = m2*m2; m1m2 = m1*m2;
		sigma1sq = ListCorrelate[w, img1*img1] - m1sq;
		sigma2sq = ListCorrelate[w, img2*img2] - m2sq;
		sigma12 = ListCorrelate[w, img1*img2] - m1m2; 
		ssimmap = ((c1 + 2*m1m2)*(c2 + 2*sigma12))/((c1 + m1sq + m2sq)*(c2 + sigma1sq + sigma2sq));
		mssim = Mean[Mean[ssimmap]];
		Return[mssim]
	];


HamDist[x_List, y_List] := 
	Module[{xfl = Flatten[x], yfl = Flatten[y]},
		Return[N[HammingDistance[xfl, yfl]/Length[xfl], 10]]
	]; 


mydescriptor[img_,winsizes_,steps_,histbins_]:=
	Module[{dim=Length[img],allwinds={},winds,i,j,i1,nowsize,nowstep},
		If[Length@winsizes!=Length@steps,Print["Error in winsizes and steps lists."];Return[]];
		For[i1=1,i1<=Length[winsizes],i1++,
			nowsize=winsizes[[i1]];
			nowstep=steps[[i1]];
			winds=Table[Take[img,{i,nowsize+i-1},{j,nowsize+j-1}],{i,1,dim-nowsize+1,nowstep},{j,1,dim-nowsize+1,nowstep}];
			AppendTo[allwinds,Flatten[Table[(Length/@BinLists[Sort[Flatten[winds[[i,j]]]],{histbins}])/(nowsize^2),{i,1,Length[winds]},{j,1,Length[winds]}]]];
		];
		Return[Flatten[allwinds]]
	];


mymail[to_, bd_, sub_: "Results from Mathematica!"] :=
 SendMail["To" -> to, "Subject" -> StringJoin[ToString@# <> "_" & /@ IntegerPart@DateList[]]<>sub, 
  "Body" -> 
   "MyMailFunc Mathem " <> ToString[DateList[][[3]]] <> "-" <> 
    ToString[DateList[][[2]]] <> "-" <> ToString[DateList[][[1]]] <> 
    " at " <> ToString[DateList[][[4]]] <> ":" <> 
    ToString[DateList[][[5]]] <> " \n \n \n" <> ToString[bd], 
  "Server" -> "smtp.live.com", 
  "From" -> "research_acc_for_mathematica@hotmail.com", 
  "UserName" -> "research_acc_for_mathematica@hotmail.com", 
  "Password" -> "researchLEME15", "PortNumber" -> 587, 
  "EncryptionProtocol" -> "StartTLS"];


(* READY CODE FROM WOLFRAM PROJECTS *)
Kruskal[pts_] := Module[{n = Length[pts], vpairs, jj = 0, hh, pair, dist, c1, c2, c1c2}, 
Do[hh[k] = {k}, {k, n}]; 
vpairs = Sort[Flatten[Table[{Norm[pts[[k]] - pts[[l]]], {k, l}}, {k, 1, n - 1}, {l, k + 1, n}], 1]]; 
First[Last[Reap[While[jj < Length[vpairs], jj++; {dist, pair} = vpairs[[jj]]; {c1, c2} = {hh[pair[[1]]], hh[pair[[2]]]}; 
If[c1 =!= c2, Sow[Rule @@ vpairs[[jj,2]]]; c1c2 = Union[c1, c2]; Do[hh[c1c2[[k]]] = c1c2, {k, Length[c1c2]}]; If[Length[hh[pair[[1]]]] == n, Break[]]; ]; 
]]]]]; 


mydescriptorVer2[img_List?MatrixQ, winsizes_List, steps_List, histbins_List]:=
	Module[{allwinds, currsize, currstep, winds},
		If[Length[winsizes] != Length[steps],
			Message[Descriptor::lengtherr];
			Return[$Failed];
		];
		
		allwinds = Map[(
			{currsize, currstep} = #;
			
			winds = Partition[img, {currsize, currsize}, currstep];
			
			(* note: Length/@BinLists is ~2x faster than BinCounts here *)
			Flatten@Map[(Length/@BinLists[Flatten[#], {histbins}]) / (currsize^2)&, winds, {2}]
		)&, Transpose[{winsizes, steps}]];
		
		Return[Flatten[allwinds]];
	];


SSIMcompiled = Compile[{{img1, _Integer, 2}, {img2, _Integer, 2}},
	Module[{w, c1, c2, m1, m2, m1sq, m2sq, m1m2, sigma1sq, sigma2sq, sigma12, ssimmap},
		(* like w = N[GaussianMatrix[{5, 1.5}], 6] - this is 15-20% faster *)
		w = {{0.00000380292, 0.0000175945, 0.0000663611, 0.000194557, 0.000412241, 0.000560994, 0.000412241, 0.000194557, 0.0000663611, 0.0000175945, 0.00000380292}, 
			 {0.0000175945, 0.0000814022, 0.000307025, 0.000900134, 0.00190726, 0.00259548, 0.00190726, 0.000900134, 0.000307025, 0.0000814022, 0.0000175945}, 
			 {0.0000663611, 0.000307025, 0.001158, 0.00339503, 0.00719362, 0.00978936, 0.00719362, 0.00339503, 0.001158, 0.000307025, 0.0000663611}, 
			 {0.000194557, 0.000900134, 0.00339503, 0.00995356, 0.0210903, 0.0287005, 0.0210903, 0.00995356, 0.00339503, 0.000900134, 0.000194557}, 
			 {0.000412241, 0.00190726, 0.00719362, 0.0210903, 0.0446874, 0.0608124, 0.0446874, 0.0210903, 0.00719362, 0.00190726, 0.000412241}, 
			 {0.000560994, 0.00259548, 0.00978936, 0.0287005, 0.0608124, 0.0827559, 0.0608124, 0.0287005, 0.00978936, 0.00259548, 0.000560994}, 
			 {0.000412241, 0.00190726, 0.00719362, 0.0210903, 0.0446874, 0.0608124, 0.0446874, 0.0210903, 0.00719362, 0.00190726, 0.000412241}, 
			 {0.000194557, 0.000900134, 0.00339503, 0.00995356, 0.0210903, 0.0287005, 0.0210903, 0.00995356, 0.00339503, 0.000900134, 0.000194557}, 
			 {0.0000663611, 0.000307025, 0.001158, 0.00339503, 0.00719362, 0.00978936, 0.00719362, 0.00339503, 0.001158, 0.000307025, 0.0000663611}, 
			 {0.0000175945, 0.0000814022, 0.000307025, 0.000900134, 0.00190726, 0.00259548, 0.00190726, 0.000900134, 0.000307025, 0.0000814022, 0.0000175945}, 
			 {0.00000380292, 0.0000175945, 0.0000663611, 0.000194557, 0.000412241, 0.000560994, 0.000412241, 0.000194557, 0.0000663611, 0.0000175945, 0.00000380292}};
		c1 = 0.01^2;
		c2 = 0.03^2;
		
		m1 = ListCorrelate[w, img1];
		m2 = ListCorrelate[w, img2];
		
		m1sq = m1^2;
		m2sq = m2^2;
		m1m2 = m1*m2;
		
		sigma1sq = ListCorrelate[w, img1^2] - m1sq;
		sigma2sq = ListCorrelate[w, img2^2] - m2sq;
		sigma12 = ListCorrelate[w, img1*img2] - m1m2;
		
		ssimmap = ((c1 + 2*m1m2)*(c2 + 2*sigma12)) / ((c1 + m1sq + m2sq)*(c2 + sigma1sq + sigma2sq));
		
		Mean[Mean[ssimmap]]
	]
, RuntimeOptions -> "Speed"];


ApproxInfoDist[img1_List?MatrixQ, img2_List?MatrixQ]:=
	Module[{x, y, xy},
		x = Total[Unitize[img1], 2]; 
		y = Total[Unitize[img2], 2];
		xy = Total[Unitize[img1 + img2], 2];

		Return[N[(2*xy - x - y)/xy]];
	];



mds[delta_List,dimensions_,dbg_,dim1_:1,dim2_:2,dim3_:3]:= (* DEBUG options= -1: PRINT NOTHING 0: all - 1: eigenvals+sol+MDSError *)
	Module[{n,deltasq,deltatotals,sumOfDelta,bMatr,eigensys,diagMatr,solpts,dMatr,mdserrors,pointstress,nofmdserrors,pos,step,finalpos},
		Monitor[
			n=Length[delta];
			deltasq=delta^2;
			step=1;
			deltatotals=Plus@@deltasq;
			sumOfDelta=Plus@@deltatotals;
			step=2;
			
			bMatr = -0.5*(deltasq - (1/n)*ConstantArray[deltatotals, n] - (1/n)*(ConstantArray[#, n] & /@deltatotals) + (1/n^2)* ConstantArray[sumOfDelta, {n, n}]);

			(*bMatr=Table[-(1/2)*(deltasq[[i, j]]-(1/n)*deltatotals[[i]]-(1/n)*deltatotals[[j]]+(1/(n^2))*sumOfDelta), {i, 1, n}, {j, 1, n}];*)			
			(*jMatr=IdentityMatrix[Length[delta]]-(1/Length[delta])*Table[1,{k,1,Length[delta]},{l,1,Length[delta]}];*)
			(*bMatr=-1/2*jMatr.deltasq.jMatr;*)
			
			step=3;
			eigensys=N[Eigensystem[N[bMatr,20],Min[10,n]],20];  (*eigensys=N[Eigensystem[N[bMatr,15]],15];*)
			(*Print["10_First_Eigenvalues= ",eigensys[[1]]];*)
			step=4;
			pos=Select[Transpose[eigensys],#[[1]]>0&];
			(*Print["xa=",Length[pos]];*)
			If[dimensions==1,finalpos={pos[[dim1]]};];
			If[dimensions==2,finalpos={pos[[dim1]],pos[[dim2]]};];
			If[dimensions==3,finalpos={pos[[dim1]],pos[[dim2]],pos[[dim3]]};];
			If[dimensions>3,finalpos=Take[pos,dimensions];];
			(*Print[dimensions];*)
			pos=Transpose[finalpos]; (* backward compatibility *)
			step=5;
			diagMatr=Sqrt[DiagonalMatrix[pos[[1]]]];  (* "Dimensions" largest positive eigenvalues *)
			step=6;
			solpts=Transpose[pos[[2]]].diagMatr;
			step=7;
			dMatr={};(*Table[EuclideanDistance[solpts[[i1]],solpts[[i2]]],{i1,1,n},{i2,1,n}];*)
			step=8;
			mdserrors=0;(*MDSError[delta,dMatr]*)
			step=9;
			nofmdserrors=0;(*noFmdsError[delta,dMatr];*)
			step=10;
			pointstress={};(*stressPerPoint[delta,dMatr];*)

			(* ALL dbg *)
			If[dbg==0,Print["delta=",delta//MatrixForm]];
			If[dbg==0,Print["bMatr=",N[bMatr,15]//MatrixForm]];
			If[dbg==0,Print["eigensys=",eigensys]];
			If[0<=dbg<=1,Print["eigenvals=",eigensys[[1]]]];
			If[dbg==0,Print["diagMatr=",diagMatr]];
			If[0<=dbg<=1,Print["sol=",solpts]];
			If[dbg==0,Print["dMatr=",dMatr//MatrixForm]];
			If[0<=dbg<=1,Print[mdserrors//MatrixForm]];
			If[dbg==0,Print["pointStress=",pointstress//MatrixForm]];
			If[0<=dbg<=1,Print[nofmdserrors//MatrixForm]];
			(* RETURN *)
			Return[{solpts,mdserrors,Select[eigensys[[1]],#>0&],pointstress,nofmdserrors}]
		,{step,{i,j},{i1,i2}}]
	];


mdsVer2[delta_List?MatrixQ, dim_Integer?Positive, accuracy_Integer: 20]:=
	Module[{n, deltasq, deltatotals, sumOfDelta, bMatr, eigenvals, eigenvecs, solpts},
		n = Length[delta];
		
		deltasq = delta^2;
		deltatotals = Total[deltasq]/n;
		sumOfDelta = Total[deltasq, 2]/(n^2);
		
		bMatr = -0.5*(deltasq - ConstantArray[deltatotals, n] - (ConstantArray[#, n] & /@ deltatotals) + ConstantArray[sumOfDelta, {n, n}]);

		{eigenvals, eigenvecs} = Eigensystem[N[bMatr, accuracy], dim];

		If[!VectorQ[eigenvals, Positive],
			Message[MDS::dimerr];
			Return[$Failed];
		];

		solpts = Transpose[eigenvecs].Sqrt[DiagonalMatrix[eigenvals]];
		Return[solpts];
	];

getUpperHalf[matr_List]:=
	Module[{i,res={}},
		For[i=1,i<=Length[matr]-1,i++,AppendTo[res,Take[matr[[i]],{i+1,Length[matr]}]]];Return[Flatten[res]]
	];

linearRegr[x_List,y_List]:=
	Module[{q},
		q={Table[1,{i,1,Length[x]}],x};
		Return[Flatten[Inverse[q.Transpose[q]].q.Transpose[{y}]]]
	];

MDSError[deltaMatrix_List,dMatrix_List]:=
	Module[{delta,d,coefs,fdelta},
		delta=getUpperHalf[deltaMatrix];
		d=getUpperHalf[dMatrix];
		coefs=linearRegr[delta,d];
		fdelta=Map[coefs[[1]]+coefs[[2]]*# &,delta];
		Return[{{"Regr.Coefs=",coefs},{"Kruskal Error=",Sqrt[(Plus@@((fdelta-d)^2))/(Plus@@(delta^2))]},{"Raw Error=",Plus@@(Plus@@((d-delta)^2))}}]
		(*
		{"Raw Error=",(Plus@@(Plus@@((d-delta)^2)))/2}
		Sqrt[Underscript[\[Sum], i<j](a*Subscript[\[Delta], ij]+b - Subscript[d, ij])^2/Underscript[\[Sum], i<j]\(Subscript[d, ij]^2)]
		Underscript[\[Sum], i<j](|Subscript[x, i]-Subscript[x, j]|-Subscript[\[Delta], ij])^2
		*)
	];

noFmdsError[deltaMatrix_List,dMatrix_List]:=  (* Simply manually assign coefs to {0,1} *)
	Module[{delta,d,coefs,fdelta},
	delta=getUpperHalf[deltaMatrix];
	d=getUpperHalf[dMatrix];
	coefs={0,1};
	fdelta=Map[coefs[[1]]+coefs[[2]]*# &,delta];
	Return[{{"Regr.Coefs=",coefs},{"\!\(\*SqrtBox[FractionBox[\(\*UnderscriptBox[\(\[Sum]\), \(i < j\)]\*SuperscriptBox[\((a*\*SubscriptBox[\(\[Delta]\), \(ij\)] + b\\\  - \\\ \*SubscriptBox[\(d\), \(ij\)])\), \(2\)]\), \(\*UnderscriptBox[\(\[Sum]\), \(i < j\)]\\\ \*SuperscriptBox[SubscriptBox[\(d\), \(ij\)], \(2\)]\)]]\)=",Sqrt[(Plus@@((fdelta-d)^2))/(Plus@@(d^2))]},{"\!\(\*UnderscriptBox[\(\[Sum]\), \(i < j\)]\)(|\!\(\*SubscriptBox[\(x\), \(i\)]\)-\!\(\*SubscriptBox[\(x\), \(j\)]\)|-\!\(\*SubscriptBox[\(\[Delta]\), \(ij\)]\)\!\(\*SuperscriptBox[\()\), \(2\)]\)=",(Plus@@(Plus@@((d-delta)^2)))/2}}]
	];


swapMatrix[matr_List,i_,j_]:=
	Module[{res=matr,ch1,ch2},
		ch1=res[[i]];ch2=res[[j]];
		res[[i]]=ch2;res[[j]]=ch1;
		res=Transpose[res];
		ch1=res[[i]];ch2=res[[j]];
		res[[i]]=ch2;res[[j]]=ch1;
		Return[res]
	];

stressPerPoint[delta_,d_]:= (* Coefficients of Regression - Fitting points - Kruskal stress *)
	Module[{data,i},
		data=Table[{},{i,1,Length[delta]}];
		For[i=1,i<=Length[delta],i++,AppendTo[data[[i]],linearRegr[Drop[delta[[i]],{i}],Drop[d[[i]],{i}]]]];
		For[i=1,i<=Length[delta],i++,AppendTo[data[[i]],Map[data[[i,1,1]]+data[[i,1,2]]*# &,Drop[delta[[i]],{i}]]]];
		For[i=1,i<=Length[delta],i++,AppendTo[data[[i]],Sqrt[(Plus@@((data[[i,2]]-Drop[delta[[i]],{i}])^2))/(Plus@@((Drop[delta[[i]],{i}])^2))]]];
		(*Print[data//MatrixForm];*)
		Return[Table[data[[i,3]],{i,1,Length[delta]}]]
	];



End[]

EndPackage[]




