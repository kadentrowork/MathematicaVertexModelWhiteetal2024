Contained in this repo is the mathematica code for the vertex model simulations in White et al 2024, "Wound-Induced Syncytia Outpace Mononucleate Neighbors during Drosophila Wound Repair"
Explanations of the parameter choices and other aspects of the model are contained there. Developed by Kaden J Tro and M. Shane Hutson. 

The main code is titled "syncytiacodefinal", which references two helper wolfram packages, "voronoiDiagramFunctions" and "voronoiDiagramGenerator (5)". 
To run the code one can download syncytiacodefinal and voronoidiagram functions along with a premade wounded cell initialization (wound1,2,3,4),
or make a new random intialization by downloading "voronoiDiagramGenerator (5)" and running the "create new.." block in the syncytiaCodeFinal notebook. The import statements 
should be updated to point to your computers downloaded copy of voronoiDiagramGenerator and voronoiDiagramFunctions

Additionally included are movies of the model running on the included example wounds, as well as the data used in the paper noted above. 
This code is highly non-optimal, and runs rather slowly. It is not gaurunteed to be bug free, although it should run correctly for the cases used in the paper.
The author recommends usage of other, better mantained and open-source code bases for your vertex modeling neads,
such as tissue forge (https://www.nature.com/articles/s41598-023-45127-x). 

Questions may be directed to kaden_tro@berkeley.edu

How to site: 
James S. White, Jasmine J. Su, Elizabeth M. Ruark, Junmin HuaM. Shane Hutson, Andrea Page-McCaw, 2024, Wound-Induced Syncytia Outpace Mononucleate Neighbors during Drosophila Wound Repair, eLife13:RP92593
