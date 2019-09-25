(* ::Package:: *)

<<MeshTools`
<<AceGEN`

SMSInitialize[
	"PressureElement",
	"Environment"->"AceFEM"
];

SMSTemplate[
	"SMSTopology"->"S2S",
	"SMSSymmetricTangent"->False,
	"SMSNodeID"->"D",
	"SMSDomainDataNames"->{"p0 -pressure","freeze -freeze multiplier dependance"},
	"SMSDefaultData"->{1.,0.}
];

SMSStandardModule["Tangent and residual"];

X\[RightTee]SMSReal@SMSIO["All coordinates"];
u\[RightTee]SMSReal@SMSIO["All DOFs"];
\[DoubleStruckP]e=Flatten[u];
{p0,freeze}\[RightTee]SMSReal@SMSIO["All domain data"];
\[Lambda]\[DoubleRightTee]SMSReal@rdata$$["Multiplier"];

SMSDo[Ig,1,SMSInteger@SMSIO["No. integration points"]];
	
	{\[Xi],\[Eta],\[Zeta]}\[RightTee]SMSReal@SMSIO["Integration point", Ig];
	wg\[RightTee]SMSReal@SMSIO["Integration weight", Ig];
	\[Psi]=ElementShapeFunction[QuadElement,2][\[Xi],\[Eta]];
	\[DoubleStruckCapitalX]\[RightTee]SMSFreeze[\[Psi].X];
	\[DoubleStruckU]\[DoubleRightTee]\[Psi].u;
	\[DoubleStruckX]=\[DoubleStruckCapitalX]+\[DoubleStruckU];
	{g\[Xi],g\[Eta]}\[DoubleRightTee]{SMSD[\[DoubleStruckX],\[Xi]],SMSD[\[DoubleStruckX],\[Eta]]};
	pn\[RightTee]SMSFreeze[p0 g\[Xi]\[Cross]g\[Eta]];
	WP\[DoubleRightTee]-(freeze+(1-freeze)\[Lambda]) pn.\[DoubleStruckX];
	
	SMSDo[i,1,Total@SMSDOFGlobal];

		Rgi\[DoubleRightTee]SMSD[WP,\[DoubleStruckP]e,i,"Constant"->pn];
		SMSIO[wg Rgi,"Add to","Residual"[i]];
	
		SMSDo[j,1,Total@SMSDOFGlobal];
	
			Kgij\[DoubleRightTee]SMSD[Rgi,\[DoubleStruckP]e,j];
			SMSIO[wg Kgij,"Add to","Tangent"[i,j]];
	
		SMSEndDo[];

	SMSEndDo[];

SMSEndDo[];

SMSWrite[];
SMTMakeDll[]
