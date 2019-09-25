(* ::Package:: *)

<<AceGEN`
<<NDSolve`FEM`

SMSInitialize["ShellElement","Environment"->"AceFEM"];

SMSTemplate[
	"SMSTopology"->"S2S",
	"SMSNodeID"->"DFi",
	"SMSDOFGlobal"->6,
	"SMSDefaultIntegrationCode"->{22,22,21},
	"SMSDomainDataNames"->{"E","\[Nu]","h"},
	"SMSDefaultData"->{20000,0.3,0.001}
];

SMSStandardModule["Tangent and residual"];

{Em, \[Nu], h}\[RightTee]SMSReal@SMSIO["All domain data"];

XIO\[RightTee]SMSReal@SMSIO["All coordinates"];
peIO\[RightTee]SMSReal@SMSIO["All DOFs"];
pe=Flatten[peIO];

uIO=peIO[[All,{1,2,3}]];
\[Phi]IO=peIO[[All,{4,5,6}]];

SMSDo[Ig,1,SMSInteger@SMSIO["No. integration points"]];
	\[CapitalXi]={\[Xi],\[Eta],\[Zeta]}\[RightTee]SMSReal@SMSIO["Integration point", Ig];
	wgp\[RightTee]SMSReal@SMSIO["Integration weight", Ig];

	Ni\[DoubleRightTee]ElementShapeFunction[QuadElement,2][\[Xi],\[Eta]];
	XM\[DoubleRightTee]Ni.XIO;
	uM\[DoubleRightTee]Ni.uIO;
	\[Phi]\[DoubleRightTee]Ni.\[Phi]IO;

	{G\[Xi],G\[Eta]}\[DoubleRightTee]SMSD[XM,{\[Xi],\[Eta]}]\[Transpose];
	N\[Zeta]=G\[Xi]\[Cross]G\[Eta];

	El3\[DoubleRightTee]N\[Zeta]/SMSSqrt[N\[Zeta].N\[Zeta]];
	El1\[DoubleRightTee]G\[Xi]/SMSSqrt[G\[Xi].G\[Xi]];
	El2\[DoubleRightTee]El3\[Cross]El1;

	T\[DoubleRightTee]{El1,El2,El3}\[Transpose];
	X\[DoubleRightTee]XM+h/2\[Zeta] El3;
	J\[DoubleRightTee]SMSD[X,\[CapitalXi]];
	Jd\[DoubleRightTee]Det[J];

	R\[DoubleRightTee]RollPitchYawMatrix[\[Phi]];
	d\[DoubleRightTee]R.El3;

	Ft\[DoubleRightTee]SMSD[XM+uM+h/2\[Zeta] d,\[CapitalXi]].SMSInverse[J].T;
	Et\[DoubleRightTee]1/2(Ft\[Transpose].Ft-IdentityMatrix[3]);
	Clear[E33];
	Et[[3,3]]=E33;

	{\[Lambda],\[Mu]}\[DoubleRightTee]SMSHookeToLame[Em,\[Nu]];
	St=\[Lambda] Tr[Et]IdentityMatrix[3]+2\[Mu] Et;
	\[Sigma]cond=First@Solve[St[[3,3]]==0,E33];

	Wshell\[DoubleRightTee]Simplify[1/2\[Lambda] Tr[Et]^2+\[Mu] Tr[Et.Et]/.\[Sigma]cond];
	Q\[DoubleRightTee]T\[Transpose].R\[Transpose].Ft;
	Wdrill\[DoubleRightTee]1/2 Em h^3 (Q[[1,2]]-Q[[2,1]])^2;

	SMSDo[m,1,SMSNoDOFGlobal];

		Rgm\[DoubleRightTee]Jd SMSD[Wshell+Wdrill,pe,m];
		SMSIO[wgp Rgm,"Add to","Residual"[m]];		

		SMSDo[n,m,SMSNoDOFGlobal];

			Kgmn\[DoubleRightTee]SMSD[Rgm,pe,n];
			SMSIO[wgp Kgmn,"Add to","Tangent"[m,n]];

		SMSEndDo[];

	SMSEndDo[];

SMSEndDo[];

SMSWrite[];
SMTMakeDll[]
