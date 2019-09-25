(* ::Package:: *)

<<MeshTools`
<<AceGEN`

SMSInitialize["GapElement","Environment"->"AceFEM"];

SMSTemplate[
	"SMSTopology"->"C1" ,
	"SMSNoNodes"->3,
	"SMSDOFGlobal"->{3,3,1},
	"SMSDefaultIntegrationCode"->0,
	"SMSAdditionalNodes"->Hold[{Null}&],
	"SMSNodeID"->{"D","D","Lagrange -LP -L"},
	"SMSDomainDataNames"->{"k -spring constant"},
	"SMSDefaultData"->{1}
];

SMSStandardModule["Tangent and residual"];

X\[RightTee]SMSIO["All coordinates"][[{1,2}]];
u\[RightTee]SMSIO["All DOFs"][[{1,2}]];
\[Lambda]\[RightTee]SMSIO["All DOFs"][[3,1]];
\[DoubleStruckP]=SMSIO["Nodal DOFs"]//Flatten;
{k}\[RightTee]SMSIO["All domain data"];

x\[DoubleRightTee]X+u;
\[DoubleStruckN]\[DoubleRightTee]X[[2]]-X[[1]];
\[DoubleStruckT]\[DoubleRightTee]x[[2]]-x[[1]];
gap\[DoubleRightTee]\[DoubleStruckT].\[DoubleStruckN]*\[DoubleStruckT].\[DoubleStruckT];

constraint=\[Lambda]+k gap;
\[CapitalPi]constraint\[DoubleRightTee]SMSIf[constraint<-10^-8,(\[Lambda]+k/2 gap)gap,-(1/(2k)) \[Lambda]^2];

\[DoubleStruckCapitalR]\[DoubleRightTee]SMSD[\[CapitalPi]constraint,\[DoubleStruckP]];
SMSIO[\[DoubleStruckCapitalR],"Add to","Residual"];
\[DoubleStruckCapitalK]=SMSD[\[DoubleStruckCapitalR],\[DoubleStruckP]];
SMSIO[\[DoubleStruckCapitalK],"Add to","Tangent"];
SMSWrite[]; 
SMTMakeDll[]
