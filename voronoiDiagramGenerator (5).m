(* ::Package:: *)

(* GENERATE A RANDOM LIST OF POINTS LYING WITHIN THE ELLIPSE DEFINED ABOVE *)
pdf = ProbabilityDistribution[(2b*Sqrt[1-(x/a)^2])/(Pi*a*b),{x,-a,a}];
xList = RandomVariate[pdf,nCells];
centroidList = Table[{x,RandomReal[{-b*Sqrt[1-(x/a)^2],b*Sqrt[1-(x/a)^2]}]},{x,xList}];
(* FUNCTIONS NEEDED FOR THE MONTE CARLO MINIMIZATION OF ENERGY *)
A=.;
dc =.;
T =.;
stepSize=.;
nMCS =.;
isInsideEllipse[point_]:=(point[[1]]/a)^2 + (point[[2]]/b)^2 <1 -0.001
isOutsideEllipse[point_]:=(point[[1]]/a)^2 + (point[[2]]/b)^2 >1
isOutsideEllipse[point_,tolerance_]:=(point[[1]]/a)^2 + (point[[2]]/b)^2 >1+tolerance
(* above two functions are slightly different in how they treat points ON the ellipse *)

distanceSquaredToEllipse[point_,\[Theta]_]:=(a*Cos[\[Theta]]-point[[1]])^2 + (b*Sin[\[Theta]]-point[[2]])^2

minDistanceSquaredToEllipse[point_]:=Minimize[distanceSquaredToEllipse[point,\[Theta]],\[Theta]][[1]]

getClosestPointOnEllipse[point_]:={a*Cos[\[Theta]],b*Sin[\[Theta]]}/.(Minimize[distanceSquaredToEllipse[point,\[Theta]],\[Theta]][[2]])

minDistanceSquaredToEllipsePointList[point_, ellipsePointList_]:=Min[SquaredEuclideanDistance[point,#]&/@ellipsePointList]

sumInverseSquaredDistances[i_,displacement_]:=Total[1/SquaredEuclideanDistance[centroidList[[i]]+displacement,#]&/@ Drop[centroidList,{i}]]

deltaEnergy1[i_,displacement_]:=(A^2)(sumInverseSquaredDistances[i,displacement]-sumInverseSquaredDistances[i,{0,0} ])

deltaEnergy2[i_, displacement_,useEllipsePointList_]:=If[useEllipsePointList,
(dc^2)(1/minDistanceSquaredToEllipsePointList[centroidList[[i]]+displacement, ellipsePointList]-1/minDistanceSquaredToEllipsePointList[centroidList[[i]], ellipsePointList]),
(dc^2)(1/minDistanceSquaredToEllipse[centroidList[[i]]+displacement]-1/minDistanceSquaredToEllipse[centroidList[[i]]])
]

deltaEnergy[i_, displacement_,useEllipsePointList_]:= 
deltaEnergy1[i, displacement]+deltaEnergy2[i, displacement,useEllipsePointList]
(* RUN MONTE CARLO MINIMIZATION OF ENERGY TO SPREAD POINTS NEARLY EVENLY WITHIN ELLIPSE *)
verbose = False;
useEllipsePointList = True;

A=100;
dc = 50;
T = 0.001;
stepSize=1;
nMCS = 200;

originalCentroidList=centroidList;

If[verbose,Dynamic["cell #"<>ToString[i]<>", \[CapitalDelta]r= "<>ToString[\[CapitalDelta]r]<>", \[CapitalDelta]E= "<>ToString[N[\[CapitalDelta]E]]<>", accept="<>ToString[accept]] ]
mcs=0;
nAccepted =0;
Quiet[Do[
(* A single Monte Carlo step *)
mcs+=1;
iter = 0;
Do[
iter+=1;
(* Choose a cell center at random *)
i = RandomInteger[{1,nCells}];
(* Get a random unit displacement *)
\[Phi]=RandomReal[{0,2\[Pi]}];
\[CapitalDelta]r = {stepSize*Cos[\[Phi]],stepSize*Sin[\[Phi]]};
(* Would this step remain within the ellipse? If not, skip to next step. *)
If[!isInsideEllipse[centroidList[[i]]+\[CapitalDelta]r],Continue[]];
(* Calculate the energy change associated with moving cell i by \[CapitalDelta]r *)
(* \[CapitalDelta]E=If[useDistToEllipseEnergy,deltaEnergy[i,\[CapitalDelta]r],deltaEnergy1[i,\[CapitalDelta]r]]; *)
\[CapitalDelta]E=deltaEnergy[i,\[CapitalDelta]r,useEllipsePointList];
accept = If[\[CapitalDelta]E<0,True,RandomReal[]<Exp[-\[CapitalDelta]E/T]];
If[!accept,Continue[] ];
nAccepted +=1;
acceptanceRate = nAccepted/((mcs-1)*nCells + iter);
centroidList[[i]]=centroidList[[i]]+\[CapitalDelta]r
,{nCells}]
,{nMCS}]]

(* BEGIN WITH AN UNBOUNDED VORONOI DIAGRAM BASED ON EQUILIBRATED CELL CENTERS ABOVE *)
Needs["ComputationalGeometry`"]
unboundedVoronoiDiagram=VoronoiDiagram[centroidList];
(* REPLACE RAYS WITH POINTS AT WHICH EACH RAY INTERSECTS THE BOUNDING ELLIPSE *)(* Some rays are completely outside ellipse. For now, these are just set equal to their final point *)rayIndexList = Flatten[Position[unboundedVoronoiDiagram[[1]],Ray[_,__]]];
rayIntersectionList={};
Do[
{x1,y1} = unboundedVoronoiDiagram[[1,index,1]];
{x2,y2} = unboundedVoronoiDiagram[[1,index,2]];
yLine = y1 + (x-x1)*(y2-y1)/(x2-x1);
intersections={x,yLine}/.Solve[yLine^2 == (b^2)(1-(x/a)^2),x];
internalIntersections=Select[intersections,Min[x1,x2]<=#[[1]]<=Max[x1,x2]&,1];
primaryIntersection = If[Length[internalIntersections]>=1, internalIntersections[[1]],{x2,y2}];
If[verbose,Print[index,": " ,{{x1,y1},{x2,y2}}]];
If[verbose,Print["       " ,intersections,", " ,primaryIntersection]];
AppendTo[rayIntersectionList,primaryIntersection]
,{index,rayIndexList}]

replacementList=Table[{1,rayIndexList[[i]]}-> rayIntersectionList[[i]],{i,Length[rayIndexList]}];
boundedVoronoiDiagram=ReplacePart[unboundedVoronoiDiagram,replacementList];
(* FUNCTIONS NEEDED FOR TRIMMING THE VORONOI DIAGRAM *)

(* function for finding vertices connecting to any vertex in a cell *)
getAdjVertices[cellAdjList_,index_]:=(indexInAdjList=Flatten[Position[cellAdjList,index]][[1]];
prevVertex=If[indexInAdjList>1,cellAdjList[[indexInAdjList-1]],cellAdjList[[-1]] ];
nextVertex=If[indexInAdjList<Length[cellAdjList],cellAdjList[[indexInAdjList+1]],cellAdjList[[1]] ];
{prevVertex, nextVertex})

(* function for finding intersection of a line segment with the ellipse boundary *)
intersectionWithBoundingEllipse[{x1_,y1_},{x2_,y2_}]:=( 
yLine = y1 + (x-x1)*(y2-y1)/(x2-x1);
intersections={x,yLine}/.Solve[yLine^2 == (b^2)(1-(x/a)^2),x];
internalIntersections=Select[intersections,Min[x1,x2]<=#[[1]]<=Max[x1,x2]&,1];
primaryIntersection = If[Length[internalIntersections]>=1, internalIntersections[[1]],{x2,y2}])

(* function for inserting a new vertex between two existing vertices in an adjacency list *)
(* returns original adjacency list if existing vertices not found in  list *)
adjListAfterInsertVertex[cellAdjList_,i1_,i2_,iNew_]:=( 
posList1=Flatten[Position[cellAdjList,i1]];
posList2=Flatten[Position[cellAdjList,i2]];
If[Length[posList1]>0 && Length[posList2]>0,
(indexInAdjList1=posList1[[1]];
indexInAdjList2=posList2[[1]];
If[Abs[indexInAdjList1-indexInAdjList2]>1 && Min[indexInAdjList1,indexInAdjList2]==1,
Insert[cellAdjList,iNew,1],
Insert[cellAdjList,iNew,Max[indexInAdjList1,indexInAdjList2]]
]),
cellAdjList
]
)

(* function to delete and decrement indices in the adjacency lists *)
deleteIndexAndDecrementThoseAbove[indexList_,indexToDelete_]:=If[#>indexToDelete,#-1,#]&/@Select[indexList,#!=indexToDelete&]
(* ADD VERTICES AT ELLIPSE BOUNDARY FOR ALL SEGMENTS THAT CROSS THE BOUNDARY *)
(*I think there's a bug here where it adds a duplicate occaisionally that ruins the intercalation code. I jury-rigged a fix at the bottom (if there is more than there should be, it deletes the extra), but it would probably be worthwhile fixing it in a better way by finding the actual bug*)
verbose = False;
outsideVertices=Select[boundedVoronoiDiagram[[1]],isOutsideEllipse];
intersectionList={};
newBoundedVoronoiDiagram = boundedVoronoiDiagram;
Do[
index=Flatten[Position[newBoundedVoronoiDiagram[[1]],outsideVertex]][[1]];
relevantCells=Select[newBoundedVoronoiDiagram[[2]],MemberQ[#[[2]],index]&];
connectingVertexIndices =Union[Flatten[getAdjVertices[#[[2]], index]&/@relevantCells]]; 
If[verbose,Print[index,": ",relevantCells]];
If[verbose,Print["    adjacent vertices= ",connectingVertexIndices]];
Do[
adjVertex = newBoundedVoronoiDiagram[[1,vIndex]];
If[!isOutsideEllipse[adjVertex],(If[verbose,Print["         edge ",index,"-",vIndex," crosses ellipse boundary"]];
If[verbose,Print["         end points = ",{outsideVertex,adjVertex}]];
intersection=intersectionWithBoundingEllipse[adjVertex,outsideVertex];
If[intersection != adjVertex,
(If[verbose,Print["          intersection point = ",intersection]];
AppendTo[intersectionList,intersection];
AppendTo[newBoundedVoronoiDiagram[[1]],intersection];
newIndex=Flatten[Position[newBoundedVoronoiDiagram[[1]],intersection]][[1]];
If[verbose,Print["        adding new vertex # ",newIndex]];
Do[
newAdjList = adjListAfterInsertVertex[cell[[2]],index,vIndex,newIndex];
newBoundedVoronoiDiagram[[2,cell[[1]],2 ]]=newAdjList; 
,{cell,relevantCells}];
relevantCells=Select[newBoundedVoronoiDiagram[[2]],MemberQ[#[[2]],index]&];
If[verbose,Print["    new AdjLists: ",relevantCells]];)
]
)
];
,{vIndex,connectingVertexIndices}]
,{outsideVertex, outsideVertices[[1;;]]}]
(* FINALLY REMOVE VERTICES THAT LIE OUTSIDE THE ELLIPSE *)
newBoundedVoronoiDiagram2=newBoundedVoronoiDiagram;

(* Remove outside vertices from the adjacency lists *)
Do[
cellAdjList =newBoundedVoronoiDiagram[[2,i]];
outsideList=isOutsideEllipse[#,0.001]&/@newBoundedVoronoiDiagram[[1,cellAdjList[[2]] ]];
newBoundedVoronoiDiagram[[2,i,2]]=Pick[cellAdjList[[2]],!#&/@outsideList];
,{i,Length[newBoundedVoronoiDiagram2[[2]]]}]


(* get a list of all vertices that are now not in the adjacency lists *)
unconnectedVertices=Complement[Range[Length[newBoundedVoronoiDiagram[[1]]] ],Flatten[newBoundedVoronoiDiagram[[2,All,2]]] ];
(* the above should return an ordered list unconnected vertices, order is important for functions below *)

(* Delete the unconnected vertices themselves *)
numDeleted=0;
Do[
newVertList = Delete[newBoundedVoronoiDiagram[[1]],indexToDelete-numDeleted];
newAdjList = {#[[1]],deleteIndexAndDecrementThoseAbove[#[[2]],indexToDelete-numDeleted]}&/@newBoundedVoronoiDiagram[[2]];
newBoundedVoronoiDiagram = {newVertList,newAdjList};
numDeleted +=1;
,{indexToDelete,unconnectedVertices}]

(* check that no unconnected vertices remain *)
unconnectedVertices=Complement[Range[Length[newBoundedVoronoiDiagram[[1]]] ],Flatten[newBoundedVoronoiDiagram[[2,All,2]]] ];

cellularVoronoiDiagram=newBoundedVoronoiDiagram;

numVertices = 2originalNCells-2;
numVertices = 2originalNCells-2;
While[Length[cellularVoronoiDiagram[[1]]]>numVertices, cellularVoronoiDiagram[[1]] = Delete[cellularVoronoiDiagram[[1]],-1]];
Quiet[
Do[Do[
If[cellularVoronoiDiagram[[2]][[i]][[2]][[j]]>numVertices,cellularVoronoiDiagram[[2]][[i]][[2]]=Delete[cellularVoronoiDiagram[[2]][[i]][[2]],j],If[cellularVoronoiDiagram[[2]][[i]][[2]][[j]]==cellularVoronoiDiagram[[2]][[i]][[2]][[If[j-1==0,-1,j-1]]], cellularVoronoiDiagram[[2]][[i]][[2]]=Delete[cellularVoronoiDiagram[[2]][[i]][[2]],j];Print["Deletion"]]],
{j, Length[cellularVoronoiDiagram[[2]][[i]][[2]]]}],{i, Length[cellularVoronoiDiagram[[2]]]}]];
