(* ::Package:: *)

<<MeshTools`
<<AceGEN`

SMSInitialize[
	"TangentBoundaryCondition",
	"Environment"->"AceFEM",
	"Mode"->"Optimal"
];

SMSTemplate[
	"SMSTopology"->"C2",
	"SMSNoNodes"->6,
	"SMSDOFGlobal"->{3,3,3,1,1,1},
	"SMSDefaultIntegrationCode"->5,
	"SMSAdditionalNodes"->"{#1,#2,#3}&",
	"SMSNodeID"->{"D","D","D","Lagr","Lagr","Lagr"},
	"SMSReferenceNodes"->{
		{-1.,0.,0.},{0.,0.,0.},{1.,0.,0.},
		{-1.,0.,0.},{0.,0.,0.},{1.,0.,0.}
	},
	"SMSSymmetricTangent"->True,
	"SMSDomainDataNames"->{
		"X0 -sphere center x",
		"Y0 -sphere center y",
		"\[Alpha] -normal angle",
		"\[Rho] -regularization parameter"
	},
	"SMSDefaultData"->{0.,0.,\[Pi]/6.,1.}
];

SMSStandardModule["Tangent and residual"];

SMSDo[IpIndex,1,SMSInteger@SMSIO["No. integration points"]];

	(*Element data*)
	{X0,Y0,\[Alpha],\[Rho]}\[DoubleRightTee]SMSReal@SMSIO["All domain data"];
	{Xi,Yi,Zi}\[DoubleRightTee]SMSReal@SMSIO["All coordinates"][[{1,2,3}]];
	{ui,vi,wi}\[DoubleRightTee]SMSReal@SMSIO["All DOFs"][[{1,2,3}]];
	\[Lambda]Ni\[DoubleRightTee]SMSReal@SMSIO["All DOFs"][[{4,5,6}]];
	at=Join[Transpose[{ui,vi,wi}],\[Lambda]Ni]//Flatten;
	(*Numerical integration*)
	\[Xi]\[RightTee]SMSReal@SMSIO["Integration point", IpIndex];
	wGauss\[RightTee]SMSReal@SMSIO["Integration weight", IpIndex];
	(*Shape functions*)
	Ni\[DoubleRightTee]ElementShapeFunction[LineElement,2][\[Xi]];
	{X,Y,Z,u,v,w,\[Lambda]N}\[DoubleRightTee]{Xi,Yi,Zi,ui,vi,wi,\[Lambda]Ni}.Ni;

	g\[Xi]\[DoubleRightTee]SMSD[{X,Y,Z},\[Xi]];
	Jd\[DoubleRightTee]SMSSqrt[g\[Xi].g\[Xi]];

	(*Tangent gap*)
	r0\[DoubleRightTee]SMSSqrt[X^2+Y^2];
	dr\[DoubleRightTee]SMSSqrt[(X+u)^2+(Y+v)^2]-r0;
	dz\[DoubleRightTee]Z+w;
	gN\[DoubleRightTee](dz Csc[\[Alpha]]-dr Cot[\[Alpha]] Csc[\[Alpha]])/(1+Cot[\[Alpha]]^2)

	(*Augmented Lagrangian of contact*)
	\[Lambda]Naug\[DoubleRightTee]\[Lambda]N+\[Rho] gN;
	Lagr\[DoubleRightTee](\[Lambda]N+\[Rho]/2 gN) gN;

	(*Tangent and residual*)
	SMSDo[i,1,SMSNoDOFGlobal];

		dLagr\[DoubleRightTee]Jd wGauss SMSD[Lagr,at,i];
		SMSIO[dLagr,"Add to","Residual"[i]];

		SMSDo[j,i,SMSNoDOFGlobal];

			ddLagr\[DoubleRightTee]SMSD[dLagr,at,j];
			SMSIO[ddLagr,"Add to","Tangent"[i,j]]

		SMSEndDo[];

	SMSEndDo[];

SMSEndDo[];

SMSWrite[];
SMTMakeDll[]
