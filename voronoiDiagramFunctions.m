(* ::Package:: *)

(* ::Input::Initialization:: *)
(* FUNCTIONS FOR CALCULATING CELL PARAMETERS *)
getCellPerimeter[cellNum_]:=(
subVertexList=cellularVoronoiDiagram[[1,cellularVoronoiDiagram[[2,cellNum,2]] ]];
perimeter=EuclideanDistance[subVertexList[[1]],subVertexList[[-1]]];
Do[perimeter +=EuclideanDistance[subVertexList[[i-1]],subVertexList[[i]]];,{i,2,Length[subVertexList]}];
perimeter
)

getCellArea[cellNum_] := (
subVertexList=cellularVoronoiDiagram[[1,cellularVoronoiDiagram[[2,cellNum,2]] ]];
area=0.5*Det[{subVertexList[[-1]],subVertexList[[1]]}];
Do[area +=0.5*Det[{subVertexList[[i-1]],subVertexList[[i]]}];,{i,2,Length[subVertexList]}];
area
)

getMoments[cellNum_]:=(
subVertexList=cellularVoronoiDiagram[[1,cellularVoronoiDiagram[[2,cellNum,2]] ]];
{x1,y1}=subVertexList[[-1]];
{x2,y2}=subVertexList[[1]];
momentX=(x1+x2)(x1 y2 - x2 y1)/6;
momentY=(y1+y2)(x1 y2 - x2 y1)/6;
Do[
{x1,y1}=subVertexList[[i-1]];
{x2,y2}=subVertexList[[i]];
momentX+=(x1+x2)(x1 y2 - x2 y1)/6;
momentY+=(y1+y2)(x1 y2 - x2 y1)/6;
,{i,2,Length[subVertexList]}];
{momentX,momentY}
)

getCentroid[cellNum_]:=getMoments[cellNum]/getCellArea[cellNum]
(* versions that return symbolic responses for use in ODEs *)
Clear[i];
getCellPerimeter[cellNum_, cellularVoronoiDiagram_, vertexList_]:=(
subVertexList=vertexList[[cellularVoronoiDiagram[[2,cellNum,2]] ]];
(*If[Length[Intersection[subVertexList[[All,1]][[All,0]][[All,2]](*accessing the subtext on the vertex list. Same with all the other strings of numbers after subVertexList*),Flatten[noContractilityList]]]>1,
perimeter=If[MemberQ[noContractilityList,subVertexList[[1]][[1]][[0]][[2]]]&&MemberQ[noContractilityList,subVertexList[[-1]][[1]][[0]][[2]]],0,EuclideanDistance[subVertexList[[1]],subVertexList[[-1]]]];
Do[perimeter +=If[Or[Evaluate[Sequence@@With[{i = i},(MemberQ[#,subVertexList[[i-1]][[1]][[0]][[2]]]&&MemberQ[#,subVertexList[[i]][[1]][[0]][[2]]]&/@noContractilityList)]]],0,EuclideanDistance[subVertexList[[i-1]],subVertexList[[i]]]];,{i,2,Length[subVertexList]}];
perimeter,*)
perimeter=EuclideanDistance[subVertexList[[1]],subVertexList[[-1]]];
Do[perimeter +=EuclideanDistance[subVertexList[[i-1]],subVertexList[[i]]];,{i,2,Length[subVertexList]}];
perimeter
)


(*this was for an attempt at fusion as described in White et al. It is buggy and crashes during full simulations but the mechanism for area sharing should be ok?*)
fusionList = {};

getCellArea[cellNum_, cellularVoronoiDiagram_, vertexList_] :=(
If[Or[Evaluate[Sequence@@MemberQ[cellNum]/@fusionList]],
area = 0;
fusedCellPosition = Position[fusionList,cellNum][[1]][[1]];
Do[
subVertexList=vertexList[[cellularVoronoiDiagram[[2,fusionList[[fusedCellPosition]][[j]],2]] ]];
area+=0.5*Det[{subVertexList[[-1]],subVertexList[[1]]}];
Do[area +=0.5*Det[{subVertexList[[i-1]],subVertexList[[i]]}];,{i,2,Length[subVertexList]}],
{j,Length[fusionList[[fusedCellPosition]]]}];
area,

subVertexList=vertexList[[cellularVoronoiDiagram[[2,cellNum,2]] ]];
area=0.5*Det[{subVertexList[[-1]],subVertexList[[1]]}];
Do[area +=0.5*Det[{subVertexList[[i-1]],subVertexList[[i]]}];,{i,2,Length[subVertexList]}];
area]
)

getMoments[cellNum_, cellularVoronoiDiagram_, vertexList_]:=(
subVertexList=vertexList[[cellularVoronoiDiagram[[2,cellNum,2]] ]];
{x1,y1}=subVertexList[[-1]];
{x2,y2}=subVertexList[[1]];
momentX=(x1+x2)(x1 y2 - x2 y1)/6;
momentY=(y1+y2)(x1 y2 - x2 y1)/6;
Do[
{x1,y1}=subVertexList[[i-1]];
{x2,y2}=subVertexList[[i]];
momentX+=(x1+x2)(x1 y2 - x2 y1)/6;
momentY+=(y1+y2)(x1 y2 - x2 y1)/6;
,{i,2,Length[subVertexList]}];
{momentX,momentY}
)

getCentroid[cellNum_, cellularVoronoiDiagram_, vertexList_]:=getMoments[cellNum, cellularVoronoiDiagram, vertexList]/getCellArea[cellNum, cellularVoronoiDiagram, vertexList]

A0s = Table[getCellArea[i],{i,Length[cellularVoronoiDiagram[[2]]]}];
P0s = Table[Sqrt[A0s[[i]]]*shP,{i,Length[cellularVoronoiDiagram[[2]]]}];
originalScale = Sqrt[Mean[A0s]];
originalPerimeterScale = Mean[P0s]; 
(* create an array of functions for time-dependent vertex positions *)
(* these will feed into NDSolve to find solutions of system of ODEs *)
(* be careful of 'vertexList' getting redefined *)
vertexList=Array[{Subscript[x, #][t], Subscript[y, #][t]}&,Length[cellularVoronoiDiagram[[1]]]];
oldVertexList = Array[{Subscript[x, #][t-\[Tau]], Subscript[y, #][t-\[Tau]]}&,Length[cellularVoronoiDiagram[[1]]]];
positionSubstitutionList = Flatten[Array[{Subscript[x, #][t]->cellularVoronoiDiagram[[1,#,1]], Subscript[y, #][t]->cellularVoronoiDiagram[[1,#,2]]}&,Length[cellularVoronoiDiagram[[1]]]]];

Clear[isOnBoundary];
(* general form for the force terms *)

(*Force2on1[point1_,point2_,k_,L0_,cellNum_,method_,cellularVoronoiDiagram_]:=
(k*(EuclideanDistance[point1,point2] -L0)+Fc[cellNum, method,cellularVoronoiDiagram])*Normalize[point2-point1 ];*)

isOnBoundary[vIndex_,cellularVoronoiDiagram_]:=Count[Flatten[cellularVoronoiDiagram[[2,All,2]]],vIndex]<3;

isOnBoundary[vIndex_]:=Count[Flatten[cellularVoronoiDiagram[[2,All,2]]],vIndex]<3;

cellIsOnBoundary[cell_]:=Or[Evaluate[Delete[isOnBoundary[#,cellularVoronoiDiagram]&/@cellularVoronoiDiagram[[2]][[cell]][[2]],0]]];

getCellEdges[cell_] := (cellVertices = cellularVoronoiDiagram[[2]][[cell]][[2]];
Flatten[Position[edgeData,x_/;MemberQ[cellVertices,x[[2]][[1]]]&&MemberQ[cellVertices,x[[2]][[2]]],1,Heads->False]]);

Clear[edgeData]
tempEdgeData = Flatten[Table[Table[{EuclideanDistance[vertexList[[cellularVoronoiDiagram[[2]][[i]][[2]][[j]]]],vertexList[[cellularVoronoiDiagram[[2]][[i]][[2]][[If[j==-1,j+2,j+1]]]]]],{cellularVoronoiDiagram[[2]][[i]][[2]][[j]],cellularVoronoiDiagram[[2]][[i]][[2]][[If[j==-1,j+2,j+1]]]}},{j,Join[{-1},Range[Length[cellularVoronoiDiagram[[2]][[i]][[2]]]-1]]}],{i,Length[cellularVoronoiDiagram[[2]]]}],1];
tempSorted = Table[Sort[tempEdgeData[[i]][[2]]],{i,Length[tempEdgeData[[All,2]]]}];
tempLocations = GatherBy[Range@Length[tempSorted],tempSorted[[#]]&];
deletions = List/@Flatten[Table[Delete[tempLocations[[i]],1],{i,Length[tempLocations]}]];
edgeData=Delete[tempEdgeData,deletions];

Clear[energy]
Clear[forceOnVertex]
(*the order of the edges carries information on orientation. See notebook for deets before changing this function 5/31/23*)

negativePressure[vIndex_,vertexList_,cellularVoronoiDiagram_,edgeData_,A0s_,outerTension_] := 
(vertices = DeleteDuplicates[Flatten[Select[edgeData[[All,2]],MemberQ[#,vIndex]&]]];
middleEdgeVertex = Select[vertices, If[outerTension==woundTension,(!MemberQ[cellularVoronoiDiagram[[2]][[woundNum]][[2]],#]),(!isOnBoundary[#,cellularVoronoiDiagram])]&][[1]];
cells = Select[cellularVoronoiDiagram[[2]],MemberQ[#[[2]],vIndex]&];
vertex1CCW = Complement[Intersection[Select[cells,#[[2]][[If[Position[#[[2]],vIndex][[1]][[1]]+1>Length[#[[2]]],1,Position[#[[2]],vIndex][[1]][[1]]+1]]]==middleEdgeVertex&][[1]][[2]],vertices],{vIndex,middleEdgeVertex}][[1]];
vertex3CCW = Complement[vertices,{vIndex,vertex1CCW,middleEdgeVertex}][[1]];
outerTension*0.5*(r[vertexList[[vertex1CCW]]-vertexList[[vertex3CCW]]])/originalScale)

energy[cellularVoronoiDiagram_,vertexList_, A0s_, contractileTension_,lineTensions_,edgeData_,health_] :=
(Sum[0.5*health[[i]]((getCellArea[i,cellularVoronoiDiagram,vertexList]-A0s[[i]])/((originalScale)^2))^2,{i, Length[cellularVoronoiDiagram[[2]]]}]+Sum[0.5*contractileTension[[i]]*(getCellPerimeter[i,cellularVoronoiDiagram,vertexList]/originalScale)^2,{i,Length[cellularVoronoiDiagram[[2]]]}]+Sum[lineTensions[[i]]*(edgeData[[i]][[1]]/originalScale),{i,Length[edgeData]}])

forceOnVertex[vIndex_,cellularVoronoiDiagram_,vertexList_,A0s_,contractileTension_,lineTensions_,edgeData_,health_] := 
((({-D[energy[cellularVoronoiDiagram,vertexList,A0s,contractileTension,lineTensions,edgeData,health],vertexList[[vIndex]][[1]]],-D[energy[cellularVoronoiDiagram,vertexList,A0s,contractileTension,lineTensions,edgeData,health],vertexList[[vIndex]][[2]]]}+If[isOnBoundary[vIndex,cellularVoronoiDiagram],
negativePressure[vIndex,vertexList,cellularVoronoiDiagram,edgeData,A0s,outerTension],{0,0}])*(k*originalScale^2/\[Mu]))/.{Abs'->Sign,k->timeConstant})

(*This is just the pre differentiated energy for faster NDSolve startup after intercalation*)

forceOnVertexSolved[vIndex_,cellularVoronoiDiagram_,vertexList_, A0s_, contractileTension_,lineTensions_,edgeData_,health_] :=
(adjacentCells = Flatten[Position[cellularVoronoiDiagram[[2]][[All,2]],x_/;MemberQ[x, vIndex]]];
adjacentVertices = Table[cellularVoronoiDiagram[[2]][[i]][[2]][[If[(pos = Position[cellularVoronoiDiagram[[2]][[i]][[2]],vIndex][[1]][[1]]+1)>Length[cellularVoronoiDiagram[[2]][[i]][[2]]],1,pos]]],{i,adjacentCells}];
listHolder = {2,3,1};
(If[isOnBoundary[vIndex,cellularVoronoiDiagram],(*boundary case*)
If[!MemberQ[cellularVoronoiDiagram[[2]][[adjacentCells[[1]]]][[2]],adjacentVertices[[2]]],
temp = adjacentCells[[2]];
adjacentCells[[2]] = adjacentCells[[1]];
adjacentCells[[1]] = temp;
temp = adjacentVertices[[2]];
adjacentVertices[[2]] = adjacentVertices[[1]];
adjacentVertices[[1]] = temp
];
AppendTo[adjacentVertices,cellularVoronoiDiagram[[2]][[adjacentCells[[2]]]][[2]][[If[(pos = Position[cellularVoronoiDiagram[[2]][[adjacentCells[[2]]]][[2]],vIndex][[1]][[1]]-1)==0,-1,pos]]]];
adjacentEdges = Table[Position[edgeData[[All,2]],x_/;MemberQ[x,vIndex]&&MemberQ[x,adjacentVertices[[i]]],1,Heads->False][[1]][[1]],{i,Length[adjacentVertices]}];
{(*xForce = *)-Sum[(0.5*(health[[adjacentCells[[i]]]]/(originalScale)^4)*((getCellArea[adjacentCells[[i]],cellularVoronoiDiagram,vertexList]-A0s[[adjacentCells[[i]]]])*(vertexList[[adjacentVertices[[i]]]][[2]]-vertexList[[adjacentVertices[[listHolder[[i]]]]]][[2]])))+(((contractileTension[[adjacentCells[[i]]]]/(originalScale^2))*getCellPerimeter[adjacentCells[[i]],cellularVoronoiDiagram,vertexList])*(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(vertexList[[vIndex]][[1]]-vertexList[[adjacentVertices[[i]]]][[1]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])+If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[listHolder[[i]]]]}]],0,1]*(vertexList[[vIndex]][[1]]-vertexList[[adjacentVertices[[listHolder[[i]]]]]][[1]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[listHolder[[i]]]]]]]))),{i,Length[adjacentCells]}]-Sum[((If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*lineTensions[[adjacentEdges[[i]]]]/originalScale)*((vertexList[[vIndex]][[1]]-vertexList[[adjacentVertices[[i]]]][[1]])/EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])),{i,Length[adjacentEdges]}],(*yForce = *)-Sum[(0.5*(health[[adjacentCells[[i]]]]/(originalScale)^4)*((getCellArea[adjacentCells[[i]],cellularVoronoiDiagram,vertexList]-A0s[[adjacentCells[[i]]]])*(vertexList[[adjacentVertices[[listHolder[[i]]]]]][[1]]-vertexList[[adjacentVertices[[i]]]][[1]])))+(((contractileTension[[adjacentCells[[i]]]]/originalScale^2)*getCellPerimeter[adjacentCells[[i]],cellularVoronoiDiagram,vertexList])*(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(vertexList[[vIndex]][[2]]-vertexList[[adjacentVertices[[i]]]][[2]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])+If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[listHolder[[i]]]]}]],0,1]*(vertexList[[vIndex]][[2]]-vertexList[[adjacentVertices[[listHolder[[i]]]]]][[2]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[listHolder[[i]]]]]]]))),{i,Length[adjacentCells]}]-Sum[(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(lineTensions[[adjacentEdges[[i]]]]/originalScale)*((vertexList[[vIndex]][[2]]-vertexList[[adjacentVertices[[i]]]][[2]])/EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])),{i,Length[adjacentEdges]}]}+negativePressure[vIndex,vertexList,cellularVoronoiDiagram,edgeData,A0s,outerTension]/originalScale
,(*normal case*)
If[!MemberQ[cellularVoronoiDiagram[[2]][[adjacentCells[[1]]]][[2]],adjacentVertices[[2]]],
temp = adjacentCells[[3]];
adjacentCells[[3]] = adjacentCells[[2]];
adjacentCells[[2]] = temp;
temp = adjacentVertices[[3]];
adjacentVertices[[3]] = adjacentVertices[[2]];
adjacentVertices[[2]] = temp
];
adjacentEdges = Table[Position[edgeData[[All,2]],x_/;MemberQ[x,vIndex]&&MemberQ[x,adjacentVertices[[i]]],1,Heads->False][[1]][[1]],{i,Length[adjacentVertices]}];
{(*xForce = *)-Sum[(0.5*(health[[adjacentCells[[i]]]]/(originalScale)^4)*((getCellArea[adjacentCells[[i]],cellularVoronoiDiagram,vertexList]-A0s[[adjacentCells[[i]]]])*(vertexList[[adjacentVertices[[i]]]][[2]]-vertexList[[adjacentVertices[[listHolder[[i]]]]]][[2]])))+(((contractileTension[[adjacentCells[[i]]]]/originalScale^2)*getCellPerimeter[adjacentCells[[i]],cellularVoronoiDiagram,vertexList])*(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(vertexList[[vIndex]][[1]]-vertexList[[adjacentVertices[[i]]]][[1]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])+If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[listHolder[[i]]]]}]],0,1]*(vertexList[[vIndex]][[1]]-vertexList[[adjacentVertices[[listHolder[[i]]]]]][[1]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[listHolder[[i]]]]]]])))+((lineTensions[[adjacentEdges[[i]]]]/originalScale)*(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(vertexList[[vIndex]][[1]]-vertexList[[adjacentVertices[[i]]]][[1]])/EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])),{i,Length[adjacentCells]}],(*yForce = *)-Sum[(0.5*(health[[adjacentCells[[i]]]]/(originalScale)^4)*((getCellArea[adjacentCells[[i]],cellularVoronoiDiagram,vertexList]-A0s[[adjacentCells[[i]]]])*(vertexList[[adjacentVertices[[listHolder[[i]]]]]][[1]]-vertexList[[adjacentVertices[[i]]]][[1]])))+(((contractileTension[[adjacentCells[[i]]]]/originalScale^2)*getCellPerimeter[adjacentCells[[i]],cellularVoronoiDiagram,vertexList])*(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(vertexList[[vIndex]][[2]]-vertexList[[adjacentVertices[[i]]]][[2]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])+If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[listHolder[[i]]]]}]],0,1]*(vertexList[[vIndex]][[2]]-vertexList[[adjacentVertices[[listHolder[[i]]]]]][[2]])/(EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[listHolder[[i]]]]]]])))+((lineTensions[[adjacentEdges[[i]]]]/originalScale)*(If[MemberQ[Sort/@noContractilityList,Sort[{vIndex,adjacentVertices[[i]]}]],0,1]*(vertexList[[vIndex]][[2]]-vertexList[[adjacentVertices[[i]]]][[2]])/EuclideanDistance[vertexList[[vIndex]],vertexList[[adjacentVertices[[i]]]]])),{i,Length[adjacentCells]}]}+If[NumberQ[woundNum],If[MemberQ[cellularVoronoiDiagram[[2]][[woundNum]][[2]],vIndex],negativePressure[vIndex,vertexList,cellularVoronoiDiagram,edgeData,A0s,woundTension]/originalScale,0],0]
])*(k*originalScale^3))/.k->(timeConstant/originalScale)

(* function for finding vertices connecting to any vertex in a cell *)
getAdjVertices[cellAdjList_,index_]:=(indexInAdjList=Flatten[Position[cellAdjList,index]][[1]];
prevVertex=If[indexInAdjList>1,cellAdjList[[indexInAdjList-1]],cellAdjList[[-1]] ];
nextVertex=If[indexInAdjList<Length[cellAdjList],cellAdjList[[indexInAdjList+1]],cellAdjList[[1]] ];
{prevVertex, nextVertex})

(* returns indices of vertices shared by two cells *)
getCommonVertices[cellNum1_,cellNum2_, cellularVoronoiDiagram_]:=Intersection[cellularVoronoiDiagram[[2,cellNum1,2]],cellularVoronoiDiagram[[2,cellNum2,2]]]

(* returns a list of cell indices for given cell's neighbors *)
getNeighboringCells[cellNum_,cellularVoronoiDiagram_]:=Select[cellularVoronoiDiagram[[2]],(Length[Intersection[#[[2]],cellularVoronoiDiagram[[2,cellNum,2]]]]>0)&&(#[[1]]!=cellularVoronoiDiagram[[2,cellNum,1]])&   ][[All,1]]

(* returns a list of cell indices for cells that surround a wound between cellNum1 and cellNum2. List is ordered to reflect degree of contact with the wounded edge. *)
getOrderedNeighboringCells[cellNum1_,cellNum2_,cellularVoronoiDiagram_]:=(nList1=getNeighboringCells[cellNum1,cellularVoronoiDiagram];
nList2=getNeighboringCells[cellNum2,cellularVoronoiDiagram];
firstNeighbors=Intersection[nList1,nList2];
allNeighbors=Select[Union[nList1,nList2],(#!=cellNum1) && (#!=cellNum2)&];
secondNeighbors=Flatten[(Intersection[allNeighbors,getNeighboringCells[#,cellularVoronoiDiagram]])&/@firstNeighbors];
otherNeighbors=Complement[allNeighbors,Join[firstNeighbors,secondNeighbors]];
orderedNeighbors = Join[firstNeighbors,secondNeighbors, otherNeighbors]
)

(* returns a list of pairs of indices of vertices for all the edges that connect a given cell to the rest of the tissue *)
getConnectingVertexPairs[cellNum_,cellularVoronoiDiagram_]:=(surroundingCellPairs=Table[Select[cellularVoronoiDiagram[[2]],(MemberQ[#[[2]],cellularVoronoiDiagram[[2,cellNum,2,i]] ])&&(#[[1]]!=cellularVoronoiDiagram[[2,cellNum,1]])&   ][[All,1]]
,{i,Length[cellularVoronoiDiagram[[2,cellNum,2]]]}];
getCommonVertices[#[[1]],#[[2]],cellularVoronoiDiagram]&/@surroundingCellPairs)
 
(* returns a list of cells that borders a given edge *)
getCellsAdjacentToEdge[edge_,cellularVoronoiDiagram_]:=Select[cellularVoronoiDiagram[[2]],Length[Intersection[#[[2]],edge]]==2&][[All,1]]

r = RotationTransform[Pi/2];
intercolationLength = intercolationProportion*originalScale;

Clear[intercolate];
(*input edge in the same form as edgeData holds edges*)
intercolate[edge_]:=
(ClearAll[cellList,p1,p2,tempPos,tempPos2];
cellList = Join[Select[cellularVoronoiDiagram[[2]],MemberQ[#[[2]], edge[[2]][[1]]]&&MemberQ[#[[2]],edge[[2]][[2]]]&],Select[cellularVoronoiDiagram[[2]],Xor[MemberQ[#[[2]], edge[[2]][[1]]],MemberQ[#[[2]],edge[[2]][[2]]]]&]];
p1 = If[MemberQ[edge[[2]],cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]][[-1]]]&&MemberQ[edge[[2]],cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]][[1]]],cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]][[-1]],Select[cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]],MemberQ[edge[[2]],#]&][[1]]];
p2 = Complement[edge[[2]],{p1}][[1]];
didRun = False;
boundaryEdge = False;
(*make a fake cell that represents the boundary so that code doesn't crash--we'll delete it later*)

If[isOnBoundary[p1,cellularVoronoiDiagram]&&isOnBoundary[p2,cellularVoronoiDiagram],didRun= True;boundaryEdge = True;
If[!MemberQ[cellList[[2]][[2]],p1],swap = cellList[[2]];cellList[[2]] = cellList[[3]];cellList[[3]] = swap];
cell = {Length[cellularVoronoiDiagram[[2]]]+1,{cellList[[3]][[2]][[If[(tPos = Position[cellList[[3]][[2]],p2][[1]][[1]]+1)>Length[cellList[[3]][[2]]],1,tPos]]],p2,p1,cellList[[2]][[2]][[If[(tPos = Position[cellList[[2]][[2]],p1][[1]][[1]]-1)==0,-1,tPos]]]}};
cellList = Insert[cellList,cell,2];AppendTo[cellularVoronoiDiagram[[2]],cell]];

If[isOnBoundary[p1,cellularVoronoiDiagram],didRun = True;
cell = {Length[cellularVoronoiDiagram[[2]]]+1,{cellList[[2]][[2]][[If[(tPos = Position[cellList[[2]][[2]],p1][[1]][[1]]+1)>Length[cellList[[2]][[2]]],1,tPos]]],p1,cellList[[1]][[2]][[If[(tPos = Position[cellList[[1]][[2]],p1][[1]][[1]]-1)==0,-1,tPos]]]}};
AppendTo[cellList,cell];AppendTo[cellularVoronoiDiagram[[2]],cell]];

If[isOnBoundary[p2,cellularVoronoiDiagram],didRun = True;
cell = {Length[cellularVoronoiDiagram[[2]]]+1,{cellList[[1]][[2]][[If[(tPos = Position[cellList[[1]][[2]],p2][[1]][[1]]+1)>Length[cellList[[1]][[2]]],1,tPos]]],p2,cellList[[2]][[2]][[If[(tPos = Position[cellList[[2]][[2]],p2][[1]][[1]]-1)==0,-1,tPos]]]}};
AppendTo[cellList,cell];AppendTo[cellularVoronoiDiagram[[2]],cell]];

(*Using a whole different procedure if need to delete a cell*)
cellDeleted = False;
If[Length[cellList[[1]][[2]]]==3||Length[cellList[[2]][[2]]]==3,
(*additional information needed for deleted Cell*)
cellDeleted = True;
cellToDelete = Select[cellList,Length[#[[2]]]==3&][[1]];
If[cellList[[1]]==cellToDelete,switch = cellList[[1]];cellList[[1]] = cellList[[2]];cellList[[2]] = switch;switch = p1;p1 = p2;p2 = switch];
thirdVertex = Complement[cellToDelete[[2]],{p1,p2}][[1]];
positionOfEdgesToDelete = Join[Position[edgeData,x_/;MemberQ[x[[2]],p1]&&MemberQ[cellList[[2]][[2]],Complement[x[[2]],{p1}][[1]]],1,Heads->False],Position[edgeData,x_/;MemberQ[x[[2]],p2]&&MemberQ[x[[2]],thirdVertex],1,Heads->False]];
(*changing connectivity of edgeData*)
tempPos = Position[edgeData,x_/;MemberQ[x[[2]],cellList[[1]][[2]][[If[(tPos =Position[cellList[[1]][[2]],p1][[1]][[1]]-1)==0,-1,tPos]]]]&&MemberQ[x[[2]],p1],1,Heads->False][[1]][[1]];
If[Position[edgeData[[tempPos]][[2]],p1][[1]][[1]] == 2,
edgeData[[tempPos]] = {EuclideanDistance[vertexList[[edgeData[[tempPos]][[2]][[1]]]],vertexList[[p2]]],{edgeData[[tempPos]][[2]][[1]],p2}},edgeData[[tempPos]] = {EuclideanDistance[vertexList[[p2]],vertexList[[edgeData[[tempPos]][[2]][[2]]]]],{p2, edgeData[[tempPos]][[2]][[2]]}}];

tempPos2 =Position[edgeData,x_/;MemberQ[x[[2]],thirdVertex]&&!MemberQ[x[[2]],p1]&&!MemberQ[x[[2]],p2],1,Heads->False][[1]][[1]] ;
If[Position[edgeData[[tempPos2]][[2]],thirdVertex][[1]][[1]] == 2,
edgeData[[tempPos2]] = {EuclideanDistance[vertexList[[edgeData[[tempPos2]][[2]][[1]]]],vertexList[[p2]]],{edgeData[[tempPos2]][[2]][[1]],p2}},edgeData[[tempPos2]] = {EuclideanDistance[vertexList[[p2]],vertexList[[edgeData[[tempPos2]][[2]][[2]]]]],{p2, edgeData[[tempPos2]][[2]][[2]]}}];

cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]],
p1];

If[MemberQ[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p1],

cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p1];
cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p2,Position[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],thirdVertex][[1]][[1]]+1];
(*cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],Complement[edgeData[[tempPos2]][[2]],{p2}][[1]],Position[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],thirdVertex][[1]][[1]]+1];*)
cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],thirdVertex];
cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],thirdVertex],

cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],p1];
cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],p2,Position[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],thirdVertex][[1]][[1]]+1];
(*cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],Complement[edgeData[[tempPos2]][[2]],{p2}][[1]],Position[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],thirdVertex][[1]][[1]]+1];*)
cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],thirdVertex];
cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],thirdVertex]];

A0s = Delete[A0s,cellToDelete[[1]]];
cellularVoronoiDiagram[[1]] = Delete[cellularVoronoiDiagram[[1]],p1];
cellularVoronoiDiagram[[1]] = Delete[cellularVoronoiDiagram[[1]],If[p1<thirdVertex,thirdVertex-1,thirdVertex]];
cellularVoronoiDiagram[[2]] = Delete[cellularVoronoiDiagram[[2]],cellToDelete[[1]]];
healths = Delete[healths,cellToDelete[[1]]];
contractileTension = Delete[contractileTension,cellToDelete[[1]]];
larger = If[p1>thirdVertex,p1,thirdVertex];
smaller = Complement[{p1,thirdVertex},{larger}][[1]];
Do[Do[If[larger>cellularVoronoiDiagram[[2]][[i]][[2]][[j]]>smaller,cellularVoronoiDiagram[[2]][[i]][[2]][[j]] = cellularVoronoiDiagram[[2]][[i]][[2]][[j]]-1,If[cellularVoronoiDiagram[[2]][[i]][[2]][[j]]>larger,cellularVoronoiDiagram[[2]][[i]][[2]][[j]] = cellularVoronoiDiagram[[2]][[i]][[2]][[j]]-2]],{j,Length[cellularVoronoiDiagram[[2]][[i]][[2]]]}],{i,Length[cellularVoronoiDiagram[[2]]]}];

Do[cellularVoronoiDiagram[[2]][[i]][[1]] =cellularVoronoiDiagram[[2]][[i]][[1]]-1,{i,cellToDelete[[1]],Length[cellularVoronoiDiagram[[2]]]}];
vertexList = Array[{Subscript[x, #][t], Subscript[y, #][t]}&,Length[cellularVoronoiDiagram[[1]]]];
edgeData = Delete[edgeData,positionOfEdgesToDelete];
edgeLeftover = tempPos2 - Count[positionOfEdgesToDelete[[All,1]],x_/;x<tempPos2 ];
lineTensions = Delete[lineTensions,positionOfEdgesToDelete];
Do[edgeData[[i]] = {EuclideanDistance[If[larger>(tPos1 = edgeData[[i]][[2]][[1]])>smaller,vertexList[[tPos1-1]],If[larger<(tPos1 = edgeData[[i]][[2]][[1]]),vertexList[[tPos1-2]],vertexList[[tPos1]]]],If[larger>(tPos2 = edgeData[[i]][[2]][[2]])>smaller,vertexList[[tPos2-1]],If[larger<(tPos2 = edgeData[[i]][[2]][[2]]),vertexList[[tPos2-2]],vertexList[[tPos2]]]]],{If[larger>tPos1>smaller,tPos1-1,If[larger<tPos1,tPos1-2,tPos1]],If[larger>tPos2>smaller,tPos2-1,If[larger<tPos2,tPos2-2,tPos2]]}},{i,Length[edgeData]}];
If[didRun,cellularVoronoiDiagram[[2]] = Delete[cellularVoronoiDiagram[[2]],-1]]
,
(*normal intercolation--no deleted cell*)
(*changing connectivity of edgeData*)
tempPos = Position[edgeData,x_/;MemberQ[x[[2]],cellList[[1]][[2]][[If[(tPos =Position[cellList[[1]][[2]],p1][[1]][[1]]-1)==0,-1,tPos]]]]&&MemberQ[x[[2]],p1],1,Heads->False][[1]][[1]];
If[Position[edgeData[[tempPos]][[2]],p1][[1]][[1]] == 2,
edgeData[[tempPos]] = {EuclideanDistance[vertexList[[edgeData[[tempPos]][[2]][[1]]]],vertexList[[p2]]],{edgeData[[tempPos]][[2]][[1]],p2}},edgeData[[tempPos]] = {EuclideanDistance[vertexList[[p2]],vertexList[[edgeData[[tempPos]][[2]][[2]]]]],{p2, edgeData[[tempPos]][[2]][[2]]}}];

tempPos2 = Position[edgeData,x_/;MemberQ[x[[2]],cellList[[2]][[2]][[If[(tPos =Position[cellList[[2]][[2]],p2][[1]][[1]]-1)==0,-1,tPos]]]]&&MemberQ[x[[2]],p2],1,Heads->False][[1]][[1]];
If[Position[edgeData[[tempPos2]][[2]],p2][[1]][[1]] == 2,
edgeData[[tempPos2]] = {EuclideanDistance[vertexList[[edgeData[[tempPos2]][[2]][[1]]]],vertexList[[p1]]],{edgeData[[tempPos2]][[2]][[1]],p1}},edgeData[[tempPos2]] = {EuclideanDistance[vertexList[[p1]],vertexList[[edgeData[[tempPos2]][[2]][[2]]]]],{p1, edgeData[[tempPos2]][[2]][[2]]}}];

cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[1]][[1]]]][[2]],
p1];
cellularVoronoiDiagram[[2]][[cellList[[2]][[1]]]][[2]] = DeleteCases[cellularVoronoiDiagram[[2]][[cellList[[2]][[1]]]][[2]],
p2];
(*changing connectivity of cellularVoronoiDiagram*)

If[MemberQ[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p1],

cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p2,Position[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p1][[1]][[1]]+1];
cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],p1,Position[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],p2][[1]][[1]]+1],

cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p1,Position[cellularVoronoiDiagram[[2]][[cellList[[3]][[1]]]][[2]],p2][[1]][[1]]+1];
cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]] = Insert[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],p2,Position[cellularVoronoiDiagram[[2]][[cellList[[4]][[1]]]][[2]],p1][[1]][[1]]+1]];
(*deleting the fake cell we made in case the edge is on the boundary*) 
If[didRun,cellularVoronoiDiagram[[2]] = Delete[cellularVoronoiDiagram[[2]],-1]]];
)

(*input cellnumbers as flat list*)
borderAblation[cells_] := (
healths[[cells[[1]]]] = 0;
healths[[cells[[2]]]] = 0;
contractileTension[[cells[[1]]]] = 0;
contractileTension[[cells[[2]]]] = 0;
meanTensionCell[[cells[[1]]]] = 0;
meanTensionCell[[cells[[2]]]] = 0;
(*edgesAffected = Flatten[Position[edgeData[[All,2]],x_/;Length[Intersection[Flatten[{cellularVoronoiDiagram[[2]][[cells[[1]]]][[2]],cellularVoronoiDiagram[[2]][[cells[[2]]]][[2]]}],x]]==2,1,Heads->False]];*)
meanTensionEdge = Table[Total[meanTensionCell[[getCellsAdjacentToEdge[edgeData[[All,2]][[i]],cellularVoronoiDiagram]]] ],{i,Length[edgeData[[All,2]]]}];
)

(*input cellnumbers as flat list*)
(*ablateCells[cNums_,tAblated_]:=
(
cNumsSorted = Sort[cNums];
woundNum = cNumsSorted[[1]];
layers = {{cNumsSorted[[1]]}};
ablationCounter = 2;
While[Length[Flatten[layers]]<Length[cNums],AppendTo[layers,DeleteCases[DeleteDuplicates[Flatten[Table[Select[getNeighboringCells[layers[[ablationCounter-1]][[i]],cellularVoronoiDiagram],MemberQ[cNums,#]&],{i,Length[layers[[ablationCounter-1]]]}]]],x_/;Or[Evaluate[Sequence@@MemberQ[x]/@layers]]]];
ablationCounter += 1];
layersSorted = Sort/@layers;
verticesWDupes = Flatten[Table[cellularVoronoiDiagram[[2]][[cNums[[i]]]][[2]],{i,Length[cNums]}]];
verticesDeleted = Select[verticesWDupes,Count[verticesWDupes,#]>1&];
Do[Do[borderBreakdown[{layersSorted[[1]][[1]],layersSorted[[l]][[-i]]}],{i,Length[layersSorted[[l]]]}];Do[Do[layersSorted[[i]][[k]] = layersSorted[[i]][[k]]-Count[layersSorted[[l]],x_/;x<layersSorted[[i]][[k]]],{k,Length[layersSorted[[i]]]}],{i,l+1,Length[layersSorted]}],{l,2,Length[layersSorted]}];
contractileTension[[cNumsSorted[[1]]]] = 0;
healths[[cNumsSorted[[1]]]] = 0;
(*no contractility so return line tensions to mean of 0*)
leadingEdgesPositions = Flatten[Position[edgeData,x_/;MemberQ[cellularVoronoiDiagram[[2]][[cNumsSorted[[1]]]][[2]],x[[2]][[1]]]&&MemberQ[cellularVoronoiDiagram[[2]][[cNumsSorted[[1]]]][[2]],x[[2]][[2]]],1,Heads->False]];
meanTensionCell = Table[(contractileTension[[i]]/originalScale)*(originalPerimeterScale-P0s[[i]]),{i,Length[P0s]}];
meanTensionEdge = Table[Total[meanTensionCell[[getCellsAdjacentToEdge[edgeData[[All,2]][[i]],cellularVoronoiDiagram]]] ],{i,Length[edgeData[[All,2]]]}];
lineTensionsLE=Flatten[Table[proc = ItoProcess[\[DifferentialD]lineTension[t]== -1/decayRate*(lineTension[t]-meanTensionEdge[[i]])\[DifferentialD]t+\[DifferentialD]w[t],lineTension[t],{lineTension,meanTensionEdge[[i]]},t,w\[Distributed]WienerProcess[0,lineTensionVariability]];
Interpolation[MovingAverage[RandomFunction[proc,{0,tMax,0.5}],5]][t+2],{i,leadingEdgesPositions}]];
increasingLineTension =Piecewise[{{0,t<tAblated},{pull*t/tResponse-pull*tAblated/tResponse,tAblated<=t<tAblated+tResponse},{pull,t>=tAblated+tResponse}}];
Do[lineTensions[[leadingEdgesPositions[[i]]]] = lineTensionsLE[[i]]+increasingLineTension,{i,Length[leadingEdgesPositions]}];
)
createContractileWave[tAblated_]:=
(
(*Parameters for Contractile Wave*)
damageTable = Table[damageFn[Norm[getCentroid[i]-getCentroid[woundNum]]],{i,Length[cellularVoronoiDiagram[[2]]]}];
timeTable = damageTable*tResponse;
cellsToRespond = Flatten[Position[timeTable,x_/;x>timeWaveStart,1,Heads->False]];
(*layers = {cNums};
contractileCounter = 2;
While[Count[layers[[contractileCounter-1]],x_/;Count[cellularVoronoiDiagram[[2]][[x]][[2]],y_/;isOnBoundary[y,cellularVoronoiDiagram]]>0]==0,AppendTo[layers,DeleteCases[DeleteDuplicates[Flatten[Table[getNeighboringCells[layers[[contractileCounter-1]][[i]],cellularVoronoiDiagram],{i,Length[layers[[contractileCounter-1]]]}]]],x_/;Or[Evaluate[Sequence@@MemberQ[x]/@layers]]]];
contractileCounter += 1];*)
(*4 is approximate/intended number of layers minus 1;*)
tContracted = (tResponse-tCrossover*4)/4;

contractileWave = Piecewise[{{-0.5*nonDimensionalizedContractileConstant/tCrossover*t,0<t<tCrossover},{-0.5*nonDimensionalizedContractileConstant,tCrossover<t<tAblated+timeToWave},{t*1.5*nonDimensionalizedContractileConstant/tCrossover-0.5*nonDimensionalizedContractileConstant-1.5*nonDimensionalizedContractileConstant/tCrossover*(tAblated+timeToWave),tAblated+timeToWave+tCrossover>t>tAblated+timeToWave},{nonDimensionalizedContractileConstant,tContracted+tAblated+timeToWave+tCrossover>t>tAblated+timeToWave+tCrossover},{-t*nonDimensionalizedContractileConstant/tCrossover+nonDimensionalizedContractileConstant+nonDimensionalizedContractileConstant/tCrossover*(tContracted+tAblated+timeToWave+tCrossover),tContracted+tAblated+timeToWave+tCrossover<t<tContracted+tAblated+timeToWave+2*tCrossover},{0,t>tContracted+tAblated+timeToWave+2*tCrossover}}];
Do[contractileTension[[j]] = contractileTension[[j]]+contractileWave/.timeToWave->timeTable[[j]],{j,cellsToRespond}];
)*)
Clear[borderBreakdown];
borderBreakdown[cNums_]:=
(
allVerticesWRepeats = Join[cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]],cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]]];
verticesToDelete = DeleteDuplicates[Flatten[Select[allVerticesWRepeats,Count[allVerticesWRepeats,#]==2&]]];
p1 = If[MemberQ[verticesToDelete,cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]][[-1]]]&&MemberQ[verticesToDelete,cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]][[1]]],cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]][[-1]],Select[cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]],MemberQ[verticesToDelete,#]&][[1]]];
p2 = Complement[verticesToDelete,{p1}][[1]];
edgesToDelete = Join[Position[edgeData,x_/;MemberQ[verticesToDelete,x[[2]][[1]]]&&MemberQ[verticesToDelete,x[[2]][[2]]],1,Heads->False],Position[edgeData,x_/;MemberQ[x[[2]],p1]&&MemberQ[cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]],Complement[x[[2]],{p1}][[1]]]&&!MemberQ[x[[2]],p2],1,Heads->False],Position[edgeData,x_/;MemberQ[x[[2]],p2]&&MemberQ[cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]],Complement[x[[2]],{p2}][[1]]]&&!MemberQ[x[[2]],p1],1,Heads->False]];

(*changing edge Connectivity*)
tempPos = Position[edgeData,x_/;MemberQ[x[[2]],p2]&&MemberQ[cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]],Complement[x[[2]],{p2}][[1]]]&&!MemberQ[x[[2]],p1],1,Heads->False][[1]][[1]];
If[Position[edgeData[[tempPos]][[2]],p2][[1]][[1]] == 2,
edgeData[[tempPos]] = {EuclideanDistance[vertexList[[Complement[edgeData[[edgesToDelete[[3]][[1]]]][[2]],{p2}][[1]]]],vertexList[[edgeData[[tempPos]][[2]][[1]]]]],{Complement[edgeData[[edgesToDelete[[3]][[1]]]][[2]],{p2}][[1]],edgeData[[tempPos]][[2]][[1]]}},edgeData[[tempPos]] = {EuclideanDistance[vertexList[[edgeData[[tempPos]][[2]][[2]]]],vertexList[[Complement[edgeData[[edgesToDelete[[3]][[1]]]][[2]],{p2}][[1]]]]],{Complement[edgeData[[edgesToDelete[[3]][[1]]]][[2]],{p2}][[1]], edgeData[[tempPos]][[2]][[2]]}}];

tempPos2 = Position[edgeData,x_/;MemberQ[x[[2]],p1]&&MemberQ[cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]],Complement[x[[2]],{p1}][[1]]]&&!MemberQ[x[[2]],p2],1,Heads->False][[1]][[1]];
If[Position[edgeData[[tempPos2]][[2]],p1][[1]][[1]] == 2,
edgeData[[tempPos2]] = {EuclideanDistance[vertexList[[Complement[edgeData[[edgesToDelete[[2]][[1]]]][[2]],{p1}][[1]]]],vertexList[[edgeData[[tempPos2]][[2]][[1]]]]],{Complement[edgeData[[edgesToDelete[[2]][[1]]]][[2]],{p1}][[1]],edgeData[[tempPos2]][[2]][[1]]}},edgeData[[tempPos2]] = {EuclideanDistance[vertexList[[edgeData[[tempPos2]][[2]][[2]]]],vertexList[[Complement[edgeData[[edgesToDelete[[2]][[1]]]][[2]],{p1}][[1]]]]],{Complement[edgeData[[edgesToDelete[[2]][[1]]]][[2]],{p1}][[1]], edgeData[[tempPos2]][[2]][[2]]}}];

cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]] = DeleteCases[Flatten[Insert[cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]],Table[cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]][[i]],{i,Join[Range[Position[cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]],p1][[1]][[1]],Length[cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]]]],Range[Position[cellularVoronoiDiagram[[2]][[cNums[[2]]]][[2]],p1][[1]][[1]]-1]]}],Position[cellularVoronoiDiagram[[2]][[cNums[[1]]]][[2]],p1][[1]][[1]]]],x_/;x==verticesToDelete[[1]]||x==verticesToDelete[[2]]];

A0s[[cNums[[1]]]] = A0s[[cNums[[1]]]]+A0s[[cNums[[2]]]];
A0s = Delete[A0s,cNums[[2]]];
P0s[[cNums[[1]]]] = Sqrt[A0s[[cNums[[1]]]]]*shP;
affectedEdges = Complement[Flatten[Position[edgeData[[All,2]],x_/;Length[Intersection[x,allVerticesWRepeats]]==2,1,Heads->False]],Flatten[edgesToDelete]];
newTension = (contractileTension[[cNums[[1]]]]/originalScale)*(-P0s[[cNums[[1]]]]) ;
newTensionEdges = Table[newTension+If[Length[getCellsAdjacentToEdge[edgeData[[i]][[2]],cellularVoronoiDiagram]]>0,(contractileTension[[Complement[getCellsAdjacentToEdge[edgeData[[i]][[2]],cellularVoronoiDiagram],{cNums[[1]]}][[1]]]]/originalScale)*(-P0s[[Complement[getCellsAdjacentToEdge[edgeData[[i]][[2]],cellularVoronoiDiagram],{cNums[[1]]}][[1]]]])],{i,affectedEdges}];
newLineTensions = Table[
proc = ItoProcess[\[DifferentialD]lineTension[t]== -1/decayRate*(lineTension[t]-newTensionEdges[[i]])\[DifferentialD]t+\[DifferentialD]w[t],lineTension[t],{lineTension,newTensionEdges[[i]]},t,w\[Distributed]WienerProcess[0,lineTensionVariability]];
(Interpolation[MovingAverage[RandomFunction[proc,{0,tMax,0.5}],5]][t+2]+ If[Complement[getCellsAdjacentToEdge[edgeData[[affectedEdges[[i]]]][[2]],cellularVoronoiDiagram],{cNums[[1]]}][[1]]==woundNum,increasingLineTension,0]),{i,Length[affectedEdges]}
];
Do[lineTensions[[affectedEdges[[i]]]] = newLineTensions[[i]],{i,Length[affectedEdges]}];
P0s = Delete[P0s,cNums[[2]]];

cellularVoronoiDiagram[[2]] = Delete[cellularVoronoiDiagram[[2]],cNums[[2]]];
cellularVoronoiDiagram[[1]] = Delete[cellularVoronoiDiagram[[1]],List/@verticesToDelete];
edgeData = Delete[edgeData,edgesToDelete];
lineTensions = Delete[lineTensions,edgesToDelete];
healths = Delete[healths,cNums[[2]]];
contractileTension = Delete[contractileTension,cNums[[2]]];
larger = If[p1>p2,p1,p2];
smaller = Complement[{p1,p2},{larger}][[1]];
Do[counter = 1;While[counter<=Length[cellularVoronoiDiagram[[2]][[i]][[2]]],If[larger>cellularVoronoiDiagram[[2]][[i]][[2]][[counter]]>smaller,cellularVoronoiDiagram[[2]][[i]][[2]][[counter]] = cellularVoronoiDiagram[[2]][[i]][[2]][[counter]]-1,If[cellularVoronoiDiagram[[2]][[i]][[2]][[counter]]>larger,cellularVoronoiDiagram[[2]][[i]][[2]][[counter]] = cellularVoronoiDiagram[[2]][[i]][[2]][[counter]]-2,If[cellularVoronoiDiagram[[2]][[i]][[2]][[counter]]==larger||cellularVoronoiDiagram[[2]][[i]][[2]][[counter]]==smaller,cellularVoronoiDiagram[[2]][[i]][[2]]=Delete[cellularVoronoiDiagram[[2]][[i]][[2]],counter];counter +=-1]]];counter +=1;],{i,Length[cellularVoronoiDiagram[[2]]]}];

Do[cellularVoronoiDiagram[[2]][[i]][[1]] =cellularVoronoiDiagram[[2]][[i]][[1]]-1,{i,cNums[[2]],Length[cellularVoronoiDiagram[[2]]]}];
vertexList = Array[{Subscript[x, #][t], Subscript[y, #][t]}&,Length[cellularVoronoiDiagram[[1]]]];
Do[edgeData[[i]] = {EuclideanDistance[If[larger>(tPos1 = edgeData[[i]][[2]][[1]])>smaller,vertexList[[tPos1-1]],If[larger<(tPos1 = edgeData[[i]][[2]][[1]]),vertexList[[tPos1-2]],vertexList[[tPos1]]]],If[larger>(tPos2 = edgeData[[i]][[2]][[2]])>smaller,vertexList[[tPos2-1]],If[larger<(tPos2 = edgeData[[i]][[2]][[2]]),vertexList[[tPos2-2]],vertexList[[tPos2]]]]],{If[larger>tPos1>smaller,tPos1-1,If[larger<tPos1,tPos1-2,tPos1]],If[larger>tPos2>smaller,tPos2-1,If[larger<tPos2,tPos2-2,tPos2]]}},{i,Length[edgeData]}];

)

fuse[cNums_]:=
((*need to modify A0s*)
positions = Position[fusionList,x_/;Length[Intersection[x,cNums]]>0 ,1,Heads->False];
If[Length[positions]>0,
totalA0 = Sum[A0s[[i]],{i,fusionList[[Flatten[positions]]][[All,1]]}];
fusionList[[positions[[1]][[1]]]] = Flatten[fusionList[[Flatten[positions]]]];
If[Length[listTemp =Complement[cNums,fusionList[[positions[[1]][[1]]]]]]>0,totalA0+=Sum[A0s[[listTemp[[j]]]],{j,Length[listTemp]}];AppendTo[fusionList[[positions[[1]][[1]]]],listTemp[[1]]]];
fusionList = Delete[fusionList,positions[[2;;]]];
Do[A0s[[i]] = totalA0,{i,fusionList[[positions[[1]][[1]]]]}],
AppendTo[fusionList,cNums];
totalA0 = Sum[A0s[[i]],{i,fusionList[[-1]]}];
Do[A0s[[i]] = totalA0,{i,fusionList[[-1]]}]
];
)
(*can we use this to analyze shrinking cells vs border breakdown--negative pressure in LE cells from A bigger than A0 pulling in material?*)
(*link Areas of cells to model fusion before border breakdown or shrinkage, see if shrinkage happens naturally*)

