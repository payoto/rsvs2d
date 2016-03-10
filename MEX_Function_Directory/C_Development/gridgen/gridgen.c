#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <stdbool.h>

#include "gridgen.h"
#include "mex.h" /* Always include this */
// Constant declaration
//int dim() = 2;

// Global Variables declaration
int nLevels, nCells,nEdges,nVerts,nCellGrid,nEdgeGrid,nVertGrid;
// File data Arrays
int *levelSize=NULL,*cells=NULL;
cellTemplate *celldatstruct=NULL; // Array containing data from File
// Array for active template
cellTemplate *cellCurrentTemplate=NULL;
vertexTemplate *vertCurrentTemplate=NULL;
edgeTemplate *edgeCurrentTemplate=NULL;
// Arrays for final grid storage
cellTemplate *cellstruct=NULL;
vertexTemplate *vertstruct=NULL;
edgeTemplate *edgestruct=NULL;
// Array for template creation
cellTemplate *cellstructTemp=NULL;
vertexTemplate *vertstructTemp=NULL;
edgeTemplate *edgestructTemp=NULL;

// Function Declarations

// Main text body

void mexFunction(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[]){
	InitialiseGridFromFile();
	
	
	
	return;
}

/*
int main(){
InitialiseGridFromFile();
return(0);
}
*/
void InitialiseGridFromFile(){
	 //Arrays
	//printf("%i\n",sizeof(edgeTemplate));
	DataIn();
	Allocatecelldatstruct();
	DataArrayToStruct();
	GenerateGrid();
	//Checks
	/*
	int ii, jj;
	printf("Check Data: nLevels = %i; nCells = %i\n",nLevels,nCells);
	
	for (ii=0;ii<=nLevels-1;ii++){
		printf("%i, %i\n",levelSize[2*ii],levelSize[2*ii+1]);
	}
	
	printf("Check Data: Indices\n");
	for (ii=0;ii<=nCells-1;ii++){
		for (jj=0;jj<=cells[ii*(nLevels+2)]+1;jj++){
			printf("%i ",cells[(nLevels+2)*ii+jj]);
		}
		printf("\n");
	}
	*/
	
}

// Allocation functions
void Allocatecelldatstruct(){
	int ii;
	
	celldatstruct=(cellTemplate*)malloc(nCells * sizeof(cellTemplate));
	for(ii=0;ii<nCells;ii++){
		celldatstruct[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
	} 
}

void DeAllocateTemplate(int domSize[dim()], cellTemplate *cellAct, edgeTemplate *edgeAct, vertexTemplate *vertAct){
	
	
	int nCellTemplate,nEdgeTemplate,nVertTemplate,ii;
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	
	for(ii=0;ii<nCellTemplate;ii++){
		free(cellAct[ii].refineVec);
	} 
	free(cellAct);
	free(vertAct);
	free(edgeAct);

}

void DeAllocate(){
	
	free(celldatstruct);
	free(cells);
	free(levelSize);
	free(cellstruct);
	free(vertstruct);
	free(edgestruct);

}

void AllocateGridStruct(int domSize[dim()],
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut){
	
	int nCellCurr,nEdgeCurr,nVertCurr, ii;
	nCellCurr=domSize[0]*domSize[1];
	nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	nVertCurr=(domSize[0]+1)*(domSize[1]+1);

	*vertstructTempOut=(vertexTemplate*)malloc(nVertCurr * sizeof(vertexTemplate));
	*edgestructTempOut=(edgeTemplate*)malloc(nEdgeCurr * sizeof(edgeTemplate));
	*cellstructTempOut=(cellTemplate*)malloc(nCellCurr * sizeof(cellTemplate));
	
	
	if (NULL == vertstructTempOut) {
		printf( "malloc failed for vertstructTemp\n");
	
	}
	if (NULL == edgestructTempOut) {
		printf( "malloc failed for edgestructTemp\n");
	
	}
	if (NULL == cellstructTempOut) {
		printf( "malloc failed for cellstructTemp\n");
	
	}
	for(ii=0;ii<nCellCurr;ii++){
		(*cellstructTempOut)[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
	} 
}

// EXECUTION FUNCTIONS

// Import data
void DataIn(){
	// Reads in data from external files
	
	FILE *cellgridFID;
	int ii,jj;
	
	
	cellgridFID=fopen("cellgrid.dat","r");
	if((cellgridFID!=NULL)){
		fscanf(cellgridFID,"%i",&nLevels); // Reads in the number of levels
		levelSize=(int*)calloc(dim()*nLevels, sizeof(int));
		for (ii=0;ii<=nLevels-1;ii++){ // Reads in the size of those levels
			for (jj=0; jj<dim();jj++){
				fscanf(cellgridFID,"%i ",&levelSize[dim()*ii+jj]);
				// printf("%i %i\n",levelSize[2*ii],levelSize[2*ii+1]);
			}
		}
		
		fscanf(cellgridFID,"%i",&nCells); // Reads in the number of cells
		//printf("	Reading cells in\n");
		cells=(int*)calloc((2+nLevels)*nCells, sizeof(int));
		
		// for (ii=0;ii<(2+nLevels)*nCells;ii++){cells[ii]=0;}
		
		for (ii=0;ii<=nCells-1;ii++){ // reads in the content of those cells
			fscanf(cellgridFID,"%i %i",&cells[(2+nLevels)*ii],&cells[(2+nLevels)*ii+1]);
			//printf("%i %i ",cells[(2+nLevels)*ii],cells[(2+nLevels)*ii+1]);
			for (jj=1;jj<=cells[(2+nLevels)*ii];jj++){
				fscanf(cellgridFID,"%i ",&cells[(2+nLevels)*ii+jj+1]);
				//printf("%i ",cells[(2+nLevels)*ii+jj+1]);
			}
			//printf("\n");
		}
		//printf("Cell Grid Structure succesfully read in\n\n");
		} else {
		perror("Data file failed to open!");
		//exit(0);
	}
}

void DataArrayToStruct(){
	int ii,jj, rootSub;
	
	for(ii=0;ii<nCells;ii++){
		rootSub=ii*(nLevels+2);
		celldatstruct[ii].index=cells[rootSub+1];
		celldatstruct[ii].refineLvl=cells[rootSub];
		
		for (jj=0;jj<celldatstruct[ii].refineLvl;jj++){
			celldatstruct[ii].refineVec[jj]=cells[rootSub+2+jj];
		}
	}
	free(cells);
	// Checks
	/*
	for(ii=0;ii<nCells;ii++){
		printf("cell(%i): index %i :: refineLvl %i :: vec ",ii,celldatstruct[ii].index,celldatstruct[ii].refineLvl);
		for (jj=0;jj<celldatstruct[ii].refineLvl;jj++){
			printf("%i ",celldatstruct[ii].refineVec[jj]);
		}
		printf("\n");
	} 
	*/
}

// Text output for Template grid

void OutputTemplateGrid(int domSize[dim()], int lvlGrid){

	int nCellCurr,nEdgeCurr,nVertCurr, ii;
	char *fileName;
	FILE *checkFID;
	
	int *nCellP, *nEdgeP, *nVertP;
	nCellP=&nCellCurr;
	nEdgeP=&nEdgeCurr;
	nVertP=&nVertCurr;
	CalculateNumElements(domSize,nCellP,nEdgeP,nVertP);
	
	if (lvlGrid>=10){
		fileName=(char*)malloc((15+(int)(log10(lvlGrid))+2)*sizeof(char));
	} else if ((lvlGrid>=1) & (lvlGrid<10)) {
		fileName=(char*)malloc((15+(int)(log10(lvlGrid))+2)*sizeof(char));
	} else {
		fileName=(char*)malloc((15)*sizeof(char));
	}
	
	sprintf(fileName,"checkTemplate%i.m",lvlGrid);
	checkFID=fopen(fileName,"w");
	free(fileName);
	printf("\n    OUTPUT: checkTemplate%i.m... ",lvlGrid);
	for (ii=0;ii<=nEdgeCurr-1;ii++){
		fprintf(checkFID,"templateGrid.edge(%i).index=%i;\n",ii+1,edgeCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(1)=%i;\n",ii+1,edgeCurrentTemplate[ii].cellind[0]);
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(2)=%i;\n",ii+1,edgeCurrentTemplate[ii].cellind[1]);
		
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(1)=%i;\n",ii+1,edgeCurrentTemplate[ii].vertex[0]);
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(2)=%i;\n",ii+1,edgeCurrentTemplate[ii].vertex[1]);
	}
	//printf(" . ");
	for (ii=0;ii<=nCellCurr-1;ii++){
		fprintf(checkFID,"templateGrid.cell(%i).index=%i;\n",ii+1,cellCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.cell(%i).fill=%lf;\n",ii+1,cellCurrentTemplate[ii].fill);
	}
	//printf(" . ");
	for (ii=0;ii<=nVertCurr-1;ii++){
		fprintf(checkFID,"templateGrid.vertex(%i).index=%i;\n",ii+1,vertCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.vertex(%i).coord(1)=%lf;\n",ii+1,vertCurrentTemplate[ii].coord[0]);
		fprintf(checkFID,"templateGrid.vertex(%i).coord(2)=%lf;\n",ii+1,vertCurrentTemplate[ii].coord[1]);
		
	} 
	printf(" done!");
	fclose(checkFID);
	
}

void OutputGrid(int jj){

	int ii;
	char fileName[20];
	FILE *checkFID;
	sprintf(fileName,"checkGrid%i.m",jj);
	checkFID=fopen(fileName,"w");
	printf("\n    OUTPUT: Checkgrid%i.m ", jj);
	for (ii=0;ii<nEdgeGrid;ii++){
		fprintf(checkFID,"templateGrid.edge(%i).index=%i;\n",ii+1,edgestruct[ii].index);
		
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(1)=%i;\n",ii+1,edgestruct[ii].cellind[0]);
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(2)=%i;\n",ii+1,edgestruct[ii].cellind[1]);
		
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(1)=%i;\n",ii+1,edgestruct[ii].vertex[0]);
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(2)=%i;\n",ii+1,edgestruct[ii].vertex[1]);
	}
	printf(" . ");
	for (ii=0;ii<nCellGrid;ii++){
		fprintf(checkFID,"templateGrid.cell(%i).index=%i;\n",ii+1,cellstruct[ii].index);
		
		fprintf(checkFID,"templateGrid.cell(%i).fill=%lf;\n",ii+1,cellstruct[ii].fill);
	}
	printf(" . ");
	for (ii=0;ii<nVertGrid;ii++){
		fprintf(checkFID,"templateGrid.vertex(%i).index=%i;\n",ii+1,vertstruct[ii].index);
		
		fprintf(checkFID,"templateGrid.vertex(%i).coord(1)=%lf;\n",ii+1,vertstruct[ii].coord[0]);
		fprintf(checkFID,"templateGrid.vertex(%i).coord(2)=%lf;\n",ii+1,vertstruct[ii].coord[1]);
		
	} 
	printf(" done!");
	fclose(checkFID);
	
}

// Template grid main process
void GenerateTemplateGrid(int lvlGenerate){
	
	int domSize[dim()],baseRefineLvl,ii;
	ii=lvlGenerate-1;
	domSize[0]=levelSize[2*ii];
	
	domSize[1]=levelSize[2*ii+1];
	//printf("%i\n",domSize[1]);
	baseRefineLvl=lvlGenerate;
	BuildLvlTemplate(domSize, baseRefineLvl,nLevels,&cellCurrentTemplate,&vertCurrentTemplate,&edgeCurrentTemplate);
	/*
	printf("\n*****  &(cellCurrentTemplate[0].index) in function is: %d\n",&(cellCurrentTemplate[0].index));
	printf("*****  (cellCurrentTemplate[0].index) in function is: %d\n",(cellCurrentTemplate[0].index));
	printf("*****  &(*cellCurrentTemplate) in function is: %d\n",&(*cellCurrentTemplate));
	printf("*****  ADRESS OF &cellCurrentTemplate in function is: %d\n",&cellCurrentTemplate);
	printf("*****  ADRESS OF cellCurrentTemplate in function is: %d\n",cellCurrentTemplate);
	printf("*****  ADRESS OF *cellCurrentTemplate in function is: %d\n",*cellCurrentTemplate);
	*/
	OutputTemplateGrid(domSize,lvlGenerate);
	
}

void GridInitialisation(){
	
	int domSize[dim()];	
	int *nCellP, *nEdgeP, *nVertP;
	int nCellCurr, nEdgeCurr, nVertCurr;
	
	// First step initialises grid
	domSize[0]=levelSize[0];
	domSize[1]=levelSize[1];
	//printf("%i\n",domSize[1]);
	AllocateGridStruct(domSize,&cellstruct,&vertstruct,&edgestruct);
	GenerateTemplateGrid(1);
	
	nCellP=&nCellCurr;
	nEdgeP=&nEdgeCurr;
	nVertP=&nVertCurr;
	CalculateNumElements(domSize,nCellP,nEdgeP,nVertP);
	
	memcpy(vertstruct,vertCurrentTemplate,sizeof(vertexTemplate)*(nVertCurr));
	memcpy(edgestruct,edgeCurrentTemplate,sizeof(edgeTemplate)*(nEdgeCurr));
	MemCopyCellStruct(cellCurrentTemplate, cellstruct, nCellCurr);
	
	/*
	printf("\n Template in Main \n");
	for(jj=0;jj<nCellCurr;jj++){
	
		printf("%i ",cellCurrentTemplate[jj].index);
		for(kk=0;kk<nLevels;kk++){
			printf("%i ",cellCurrentTemplate[jj].refineVec[kk]);
		}
		printf("\n");
	}*/
	
	nCellGrid=nCellCurr;
	nEdgeGrid=nEdgeCurr;
	nVertGrid=nVertCurr;
	
	//printf("memory was copied\n");
	DeAllocateTemplate(domSize, cellCurrentTemplate, edgeCurrentTemplate, vertCurrentTemplate);
	//printf("Template deallocated\n");
	
	
}

void GenerateGrid(){
	int ii;//jj,kk;
	int *posCellRefine=NULL,*indCellRefine=NULL,*posEdgeRefine=NULL,*indEdgeRefine=NULL;
	int domSize[dim()];
	int nCellRefine=0, nEdgeRefine=0;
	printf("\n***** START GRID INITIALISATION *****\n");
	GridInitialisation();
	printf("\n    ACTION: Starting Grid Generated");
	OutputGrid(1);
	printf("\n\n");
	
	/*
	printf("\n Edge Orientation\n");
	for(jj=0;jj<nEdgeGrid;jj++){
		printf("%i %i",edgestruct[jj].index,edgestruct[jj].orientation);
		printf("\n");
	}*/
	for (ii=2;ii<=nLevels;ii++){
		printf("\n***** START GRID REFINEMENT %i of %i *****\n",ii-1,nLevels-1);
		GenerateTemplateGrid(ii);
		printf("\n    ACTION: Template Generated");
		IdentifyRefineCell(ii,&posCellRefine,&indCellRefine,&nCellRefine);
		IdentifyRefineEdge(posCellRefine, indCellRefine,nCellRefine,
				&posEdgeRefine,&indEdgeRefine,&nEdgeRefine,edgestruct,nEdgeGrid);
		
		printf("\n    ACTION: Refinement Targets Identified");
		domSize[0]=levelSize[2*(ii-1)];
		domSize[1]=levelSize[2*(ii-1)+1];
		/*
		printf("\n Edge Orientation start\n");
		for(jj=0;jj<nEdgeGrid;jj++){
			printf("%i %i",edgestruct[jj].index,edgestruct[jj].orientation);
			printf("\n");
		}*/
		RefineSelectedEdges(domSize,posEdgeRefine,indEdgeRefine,nEdgeRefine);
		printf("\n    ACTION: Edges Refined");
		/*
		printf("\n Edge Orientation post edge refine\n");
		for(jj=0;jj<nEdgeGrid;jj++){
			printf("%i %i",edgestruct[jj].index,edgestruct[jj].orientation);
			printf("\n");
		}*/
		RefineSelectedCells(domSize,posCellRefine,indCellRefine,nCellRefine);
		printf("\n    ACTION: Cells Refined");
		/*
		printf("\n Cells refine Vecs\n");
		for(jj=0;jj<nCellGrid;jj++){
		
			printf("%i ",cellstruct[jj].index);
			for(kk=0;kk<nLevels;kk++){
				printf("%i ",cellstruct[jj].refineVec[kk]);
			}
			printf("\n");
		}
		*/
		/*
		printf("\n Edge Orientationpost cell refine\n");
		for(jj=0;jj<nEdgeGrid;jj++){
			printf("%i %i",edgestruct[jj].index,edgestruct[jj].orientation);
			printf("\n");
		}*/
		free(posCellRefine);
		free(indCellRefine);
		free(posEdgeRefine);
		free(indEdgeRefine);
		DeAllocateTemplate(domSize, cellCurrentTemplate, edgeCurrentTemplate, vertCurrentTemplate);
		OutputGrid(ii);
		printf("\n\n");
	}
	
	
}


// Base Template grid Generation operations
void EdgeIJtoGrid(int domSize[dim()],int IJK[dim()], int l){
	// l is either 0 or 1 and is used as a mathematical switch
	// Allows the selection of the right position to assign the data
	
	int edgSub;
	
	edgSub=edgsub(IJK[0], IJK[1], l,domSize);
	
	edgestructTemp[edgSub].index=
			(1-l)*(IJK[0]+(IJK[1]-1)*domSize[0]) // case where l=0
			+(l)*((IJK[0]-1)*(domSize[1])+(IJK[1]) // case where l=1
				+((domSize[1]+1)*domSize[0]));
				
	edgestructTemp[edgSub].orientation=l;	
	
	edgestructTemp[edgSub].vertex[0]=(IJK[0]+(IJK[1]-1)*(domSize[0]+1))*(1-l)+ // l=1
			(IJK[0]+(IJK[1]-1)*(domSize[0]+1))*l;// l=2
	edgestructTemp[edgSub].vertex[1]=(1+IJK[0]+(IJK[1]-1)*(domSize[0]+1))*(1-l)+ // l=1
		(IJK[0]+(IJK[1]-1+1)*(domSize[0]+1))*l;// l=2
	
	edgestructTemp[edgSub].cellind[0]=(IJK[0]+(IJK[1]-2)*domSize[0])*(1-l)+// l=1
		(-1+IJK[0]+(IJK[1]-1)*domSize[0])*l;// l=2
	edgestructTemp[edgSub].cellind[1]=(IJK[0]+(IJK[1]-1)*domSize[0])*(1-l)+//ll=1
		(IJK[0]+(IJK[1]-1)*domSize[0])*l; // l=2
	
}

void EdgeIJtoGridLowBound(int domSize[dim()],int IJK[dim()], int l){
	// l is either 0 or 1 and is used as a mathematical switch
	// Allows the selection of the right position to assign the data
	
	int edgSub;
	
	edgSub=edgsub(IJK[0], IJK[1], l,domSize);
	
	EdgeIJtoGrid(domSize,IJK,l);
	
	edgestructTemp[edgSub].cellind[0]=(-1)*(1-l)+// l=1
		(-2)*l;// l=2
		
	
}

void EdgeIJtoGridHighBound(int domSize[dim()],int IJK[dim()], int l){
	// l is either 0 or 1 and is used as a mathematical switch
	// Allows the selection of the right position to assign the data
	
	int edgSub;
	
	edgSub=edgsub(IJK[0], IJK[1], l,domSize);
	
	EdgeIJtoGrid(domSize,IJK,l);
	
	edgestructTemp[edgSub].cellind[1]=(-3)*(1-l)+//ll=1
		(-4)*l; // l=2
	
}

void VertexIJGrid(int domSize[dim()],int IJK[dim()]){
	// For vertices IJK goes to dim()+1
	// Assigns vertex data to template structure
	int vertSub;
	
	vertSub=vertsub(IJK[0],IJK[1],domSize);
	vertstructTemp[vertSub].index=IJK[0]+(IJK[1]-1)*(domSize[0]+1);
	//printf("Accessed Index\n");
	vertstructTemp[vertSub].coord[0]=((double)(IJK[0]-1))
				/((double)(domSize[0]));
	vertstructTemp[vertSub].coord[1]=((double)(IJK[1]-1))
				/((double)(domSize[1]));
	//printf("Accessed coord\n");
	
}

void CellIJGrid(int domSize[dim()],int IJK[dim()], int baseRefineLvl){
	// For vertices IJK goes to dim()+1
	// baseRefineLvl starts from 1
	int cellSub;
	
	cellSub=cellsub(IJK[0],IJK[1],domSize);
	cellstructTemp[cellSub].index=IJK[0]+(IJK[1]-1)*domSize[0];
	cellstructTemp[cellSub].fill=0;
	cellstructTemp[cellSub].refineLvl=baseRefineLvl;
	cellstructTemp[cellSub].refineVec[baseRefineLvl-1]=IJK[0]+(IJK[1]-1)*domSize[0];

}

// Upper level grid template generation
void AssignEdgestructContent(int domSize[dim()]){
	int ii,jj,IJK[dim()];
	
	// Assign all data
	for (ii=1;ii<=domSize[0];ii++){
		for (jj=2;jj<=domSize[1];jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			EdgeIJtoGrid(domSize,IJK,0);
		}
	}
	for (ii=2;ii<=domSize[0];ii++){
		for (jj=1;jj<=domSize[1];jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			EdgeIJtoGrid(domSize,IJK,1);
		}
	}
	
	// Boundary conditions
	ii=1; 
	for (jj=1;jj<=domSize[1];jj++){
		IJK[0]=ii;
		IJK[1]=jj;
		EdgeIJtoGridLowBound(domSize,IJK,1);
	}
	ii=domSize[0]+1;
	for (jj=1;jj<=domSize[1];jj++){
		
		IJK[0]=ii;
		IJK[1]=jj;
		EdgeIJtoGridHighBound(domSize,IJK,1);
	}
	jj=1;
	for (ii=1;ii<=domSize[0];ii++){

		IJK[0]=ii;
		IJK[1]=jj;
		EdgeIJtoGridLowBound(domSize,IJK,0);
	}
	jj=domSize[1]+1;
	for (ii=1;ii<=domSize[0];ii++){

		IJK[0]=ii;
		IJK[1]=jj;
		EdgeIJtoGridHighBound(domSize,IJK,0);
	}
	
	
}

void AssignVertextructContent(int domSize[dim()]){
	
	int ii, jj, IJK[dim()];
	
	// Assign all data
	for (ii=1;ii<=domSize[0]+1;ii++){
		for (jj=1;jj<=domSize[1]+1;jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			VertexIJGrid(domSize,IJK);
		}
	}
	
	
}

void AssignCelltructContent(int domSize[dim()], int baseRefineLvl){
	
	int ii, jj, IJK[dim()];
	
	// Assign all data
	for (ii=1;ii<=domSize[0];ii++){
		for (jj=1;jj<=domSize[1];jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			CellIJGrid(domSize,IJK,baseRefineLvl);
		}
	}
	
	
}

void CalculateNumElements(int domSize[dim()],int *nCellCurr,int *nEdgeCurr,int *nVertCurr){
	
	*nCellCurr=domSize[0]*domSize[1];
	*nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	*nVertCurr=(domSize[0]+1)*(domSize[1]+1);
}

void MemCopyCellStruct(cellTemplate *original, cellTemplate *destination, int nElm){
	int ii;
	for(ii=0;ii<nElm;ii++){
		memcpy(&(destination[ii].index),&(original[ii].index),(sizeof(cellTemplate)-sizeof(int*)));
		memcpy(&(destination[ii].refineVec[0]),&(original[ii].refineVec[0]),sizeof(int)*nLevels);
	}
	
}

void BuildLvlTemplate(int domSize[dim()], int baseRefineLvl, int nLevelsInput,
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut){
	

	//int ii;
	//int jj;
	//int kk;

	int *nCellP;int *nEdgeP;int *nVertP;
	int nCellCurr; int nEdgeCurr;int nVertCurr;
	
	nCellP=&nCellCurr;
	nEdgeP=&nEdgeCurr;
	nVertP=&nVertCurr;
	
	CalculateNumElements(domSize,nCellP,nEdgeP,nVertP);
	nLevels=nLevelsInput;
	AllocateGridStruct(domSize,cellstructTempOut,vertstructTempOut,edgestructTempOut);
	AllocateGridStruct(domSize,&cellstructTemp,&vertstructTemp,&edgestructTemp);
	//printf("Allocated the things\n");
	AssignVertextructContent(domSize);
	AssignEdgestructContent(domSize);
	AssignCelltructContent(domSize,baseRefineLvl);

	
	
	//printf("Size of Copied Cell in function %i\n",sizeof(*cellstructTemp)*(nCellCurr));
	
	memcpy(*vertstructTempOut,vertstructTemp,sizeof(vertexTemplate)*(nVertCurr));
	memcpy(*edgestructTempOut,edgestructTemp,sizeof(edgeTemplate)*(nEdgeCurr));
	//memcpy(*cellstructTempOut,cellstructTemp,sizeof(*cellstructTemp)*(nCellCurr));
	MemCopyCellStruct(cellstructTemp,*cellstructTempOut,nCellCurr);
	
	DeAllocateTemplate(domSize,cellstructTemp,edgestructTemp,vertstructTemp);
	//free(cellstructTemp);
	//free(vertstructTemp);
	//free(edgestructTemp);
}


// Grid Refinement Processes

void IdentifyRefineCell(int refinLvl,int** posGridRef,int **indGridRef, int *nRefineCell){
	// Identify cells which still need refinement
	//	 This is done by going through the data and finding the 
	//	 positions (specified in the vector) where refinement is needed
	//	 For all these operations the index of cells CANNOT BE USED
		
	// celldatstruct[].refineLvl>=refinLvl -> Needs to be refined
	// celldatstruct[].refinevec[0:refinLvl-2] 
	// 				-> indicate the current position of the cell of interest
	//              => This is the universal indicator that keeps the cell consistant
	// the corresponding cells in the current grid will be:
	// cellstruct[].refinevec[0:refinLvl-2] 
	
	// The output that is needed here is the position and indices of the cells 
	// that need refinement 
	// cellstruct[kk].index
	
	int ii,kk, jj,ll, nPosRefine;
	int *posRefV=NULL;
	int *posCellRef=NULL;
	//int *posGridRef=NULL, *indGridRef=NULL;
	// first count the number of refinements needed

	posRefV= (int *)calloc(nCells,sizeof(int));
	posCellRef= (int *)calloc(nCells*(refinLvl-1),sizeof(int));
	//printf("%i \n",refinLvl-2);
	kk=0;
	for (ii=0;ii<nCells;ii++){
	
		//printf("%i %i\n",ii,kk);
		if (celldatstruct[ii].refineLvl>=refinLvl){
			posRefV[kk]=ii;
			for (ll=0;ll<refinLvl-1;ll++){
				posCellRef[kk*(refinLvl-1)+ll]=celldatstruct[ii].refineVec[ll];
			}
			kk++;
		}
	}
	//printf("%i",kk);
	posRefV=(int*)realloc(posRefV,kk*sizeof(int));
	posCellRef=(int*)realloc(posCellRef,kk*sizeof(int)*(refinLvl-1));
	
	nPosRefine=kk;
	(*posGridRef)=(int*)calloc(nCellGrid,sizeof(int));
	kk=0;
	
	// Match vectors to cell location 
	// -> Note will be more efficient to remove 
	//		vectors which are the same and then do this step
	//printf("\n");
	for (ii=0;ii<nPosRefine;ii++){
		for (jj=0;jj<nCellGrid;jj++){
			ll=0;
			while((cellstruct[jj].refineVec[ll]==posCellRef[ii*(refinLvl-1)+ll])
				& (ll<(refinLvl-1))){
				ll++;
				//printf("%i ",ll);
			}
			if ((ll)==(refinLvl-1)){
				(*posGridRef)[kk]=jj;
				kk++;
				//printf("\n");
				break;
			}
		}
	}
	
	if (kk==0){perror("ERROR: The specified refineVec does not exist");exit(EXIT_FAILURE);}
	
	(*posGridRef)=(int*)realloc((*posGridRef),kk*sizeof(int));
	/*
	// Checks - Seem to work
	printf("\n");
	for (ii=0;ii<kk;ii++) {
		printf("%i - ",(*posGridRef)[ii]);
		for (jj=0;jj<(refinLvl-1);jj++){
			printf("%i ",cellstruct[(*posGridRef)[ii]].refineVec[jj]);
		}
		printf("\n");
	}
	
	
	//printf("\n");
	//for (ii=0;ii<kk;ii++) {printf("%i ",posCellRef[ii]);}
	//printf("\n");
	*/
	/*
	qsort((*posGridRef),nPosRefine,sizeof(int),SortOrder);
	jj=1;
	printf("\n");
	for (ii=1;ii<nPosRefine;ii++){
		if ((*posGridRef)[ii]>(*posGridRef)[jj-1]){
			(*posGridRef)[jj]=(*posGridRef)[ii];
			jj++;
			//printf("%i ", jj);
		}
	}
	(*posGridRef)=(int*)realloc((*posGridRef),jj*sizeof(int));
	*/
	
	RemoveIdenticalEntries_int(posGridRef, nPosRefine, &jj);
	*indGridRef=(int*)calloc(jj,sizeof(int));
	for(ii=0;ii<jj;ii++){
		(*indGridRef)[ii]=cellstruct[(*posGridRef)[ii]].index;
	}
	/*
	printf("\n");
	for (ii=0;ii<jj;ii++){
		printf("%i - ",(*posGridRef)[ii]);
		for (ll=0;ll<(refinLvl-1);ll++){
			printf("%i ",cellstruct[(*posGridRef)[ii]].refineVec[ll]);
		}
		printf("\n");
	}
	*/
	*nRefineCell=jj;
	free(posRefV);
	free(posCellRef);
	//free(posGridRef);
	//free(indGridRef)
}

void IdentifyRefineEdge(int *posCellRefine, int *indCellRefine,int nCellRefine,
		int **posEdgeRefine, int **indEdgeRefine,int *nEdgeRefine, edgeTemplate *edgestructAct, int nEdge){
	
	int ii,jj,ll,kk,kkStart,flagEdge;
	
	*posEdgeRefine=(int*)calloc(nEdge,sizeof(int));
	kk=0;
	for(ii=0;ii<nEdge;ii++){
		for(jj=0;jj<nCellRefine;jj++){
			ll=0;
			kkStart=kk;
			do{
				flagEdge=((edgestructAct[ii].cellind[ll]==indCellRefine[jj]));
				(*posEdgeRefine)[kk]=ii*flagEdge;
				kk=kk+flagEdge;
				ll++;
			} while((kkStart==kk)&(ll<2));
		}
	}
		
	
	*posEdgeRefine=(int*)realloc(*posEdgeRefine,kk*sizeof(int));
	RemoveIdenticalEntries_int(posEdgeRefine, kk, &kk);
	*nEdgeRefine=kk;
	*indEdgeRefine=(int*)calloc(kk,sizeof(int));
	
	for(ii=0;ii<kk;ii++){
		(*indEdgeRefine)[ii]=edgestructAct[(*posEdgeRefine)[ii]].index;
	}
	/*
	for(ii=0;ii<kk;ii++){
		printf("%i ",(*indEdgeRefine)[ii]);
	}*/
}
	
void RefineSelectedEdges(int domSize[dim()],int *posEdgeRefine,int *indEdgeRefine, int nEdgeRefine){
	// This function needs to find out how many mre edges and vertices are needed.
	// Reallocate the gridstruct arrays to match that
	// Then actually do the modification
	
	int ii,jj;
	int nAddEdge,nAddVertex, maxEdgeIndex,maxVertexIndex;
	int nEdgeSplit,startEdgeSub,currEdgeSub,currVertSub,
		startVertInd, endVertInd,flipSwitch,startVertSub,endVertSub;
	double coordCooeff[dim()],coordStart[dim()],coordEnd[dim()];
	int *listIndEdge=NULL, *listIndVert=NULL;
	
	nAddEdge=0;
	
	for (ii=0;ii<nEdgeRefine;ii++){
		nAddEdge=nAddEdge
					+(1-edgestruct[posEdgeRefine[ii]].orientation)*(domSize[0]-1)
					+(edgestruct[posEdgeRefine[ii]].orientation)*(domSize[1]-1);
	}
	nAddVertex=nAddEdge;
	
	listIndEdge=(int*)ArrayFromStruct(edgestruct, nEdgeGrid, sizeof(edgeTemplate),sizeof(int));
	maxEdgeIndex=max_array(listIndEdge, nEdgeGrid);
	listIndVert=(int*)ArrayFromStruct(vertstruct, nVertGrid, sizeof(vertexTemplate),sizeof(int));
	maxVertexIndex=max_array(listIndVert, nVertGrid);
	
	edgestruct=(edgeTemplate*)realloc(edgestruct,(nEdgeGrid+nAddEdge)*sizeof(edgeTemplate));
	vertstruct=(vertexTemplate*)realloc(vertstruct,(nVertGrid+nAddVertex)*sizeof(vertexTemplate));
	
	currEdgeSub=nEdgeGrid;
	currVertSub=nVertGrid;
	
	for(ii=0;ii<nEdgeRefine;ii++){
		nEdgeSplit=(1-edgestruct[posEdgeRefine[ii]].orientation)*(domSize[0]-1)
					+(edgestruct[posEdgeRefine[ii]].orientation)*(domSize[1]-1);
		//printf("nEdgeSplit: %i \n",nEdgeSplit);

		//printf("Start Finding ObjNum ... ");
		startVertInd=edgestruct[posEdgeRefine[ii]].vertex[0];
		endVertInd=edgestruct[posEdgeRefine[ii]].vertex[1];
		FindObjNum(&startVertInd, listIndVert, &startVertSub, 1,nVertGrid);
		FindObjNum(&endVertInd, listIndVert, &endVertSub, 1,nVertGrid);
		//printf(" done!\n");
		//printf("Start Coordinate Calc ... ");
		coordStart[0]=vertstruct[startVertSub].coord[0];
		coordStart[1]=vertstruct[startVertSub].coord[1];
		coordEnd[0]=vertstruct[endVertSub].coord[0];
		coordEnd[1]=vertstruct[endVertSub].coord[1];
		coordCooeff[0]=1.0/((double)(nEdgeSplit+1))*(coordEnd[0]-coordStart[0]);
		coordCooeff[1]=1.0/((double)(nEdgeSplit+1))*(coordEnd[1]-coordStart[1]);
		//printf(" done!\n");
		
		flipSwitch=1;
		startEdgeSub=posEdgeRefine[ii];
		edgestruct[startEdgeSub].vertex[flipSwitch]=maxVertexIndex+1;
		
		for(jj=0;jj<nEdgeSplit;jj++){
		
			memcpy(&(edgestruct[currEdgeSub].index),
				&(edgestruct[startEdgeSub].index),sizeof(edgeTemplate));
				
			flipSwitch=abs(flipSwitch-1);
			maxVertexIndex++;
			maxEdgeIndex++;
			
			edgestruct[currEdgeSub].index=maxEdgeIndex;
			edgestruct[currEdgeSub].vertex[flipSwitch]=maxVertexIndex+1;
			
			vertstruct[currVertSub].index=maxVertexIndex;
			vertstruct[currVertSub].coord[0]=coordCooeff[0]*((double)(jj+1))+coordStart[0];
			vertstruct[currVertSub].coord[1]=coordCooeff[1]*((double)(jj+1))+coordStart[1];
			
			startEdgeSub=currEdgeSub;
			currEdgeSub++;
			currVertSub++;
		}
		edgestruct[startEdgeSub].vertex[flipSwitch]=endVertInd;
	}
	nEdgeGrid=(nEdgeGrid+nAddEdge);
	nVertGrid=(nVertGrid+nAddVertex);
	
	free(listIndEdge);
	free(listIndVert);
}

void ExtractEdgeChain(int *posEdge, int *indEdge,int nEdge ,int *ordPosEdge,
	int *ordIndEdge, int *ordIndVert,edgeTemplate *edgestructAct){
	
	int ii,jj,kk,ll;
	int actPos, actInd, remActPos;
	int *listIndVert=NULL;
	
	listIndVert=(int*)calloc(2*nEdge,sizeof(int));
	
	for(ii=0;ii<nEdge;ii++){
		listIndVert[2*ii]=edgestructAct[posEdge[ii]].vertex[0];
		listIndVert[2*ii+1]=edgestructAct[posEdge[ii]].vertex[1];
	}
	
	actPos=0;
	for(ii=0;ii<nEdge;ii++){
		actInd=listIndVert[actPos];
		ordIndVert[ii]=actInd;
		if (actInd==0){
			printf("\n **** active Vertex Index is 0 ii= %i ; nEdge= %i \n",ii,nEdge);
			/*
			for(ii=0;ii<nEdge;ii++){
				listIndVert[2*ii]=edgestructAct[posEdge[ii]].vertex[0];
				listIndVert[2*ii+1]=edgestructAct[posEdge[ii]].vertex[1];
				printf("%i   %i\n",listIndVert[2*ii],listIndVert[2*ii+1]);
			}
			
			
			*/
			exit(EXIT_FAILURE);
		}
		
		listIndVert[actPos]=0;
		FindObjNum(&actInd, listIndVert, &actPos, 1, 2*nEdge);
		listIndVert[actPos]=0;
		remActPos=(actPos % 2);
		ordPosEdge[ii]=posEdge[(actPos-remActPos)/2];
		ordIndEdge[ii]=indEdge[(actPos-remActPos)/2];
		actPos=actPos+(1-2*remActPos);
	}
	
	free(listIndVert);
}

int LowerLeftVertex(double *coordList, int nEdge){
	int ii, jj;
	int minPos;
	double minCoord;
	double *distList=NULL;
	
	distList=(double*)calloc(nEdge,sizeof(double));
	minCoord=fabs(fmin_array(coordList,dim()*nEdge));
	for(ii=0;ii<nEdge;ii++){
		for(jj=0;jj<dim();jj++){
			distList[ii]=distList[ii]+pow((coordList[ii*dim()+jj]+minCoord),2.0);
		}
	}
	minPos=pos_fmin_array(distList,nEdge);
	/*
	printf("\n minPos=%i ; ",minPos);
	for(ii=0;ii<nEdge;ii++){
		printf("%lf ",distList[ii]);
	}
	printf("\n");
	*/
	free(distList);
	return(minPos);
}

double* ExtractCoord(int nEdge, int *ordIndVert,int * ordPosVert, vertexTemplate *vertstructAct, int nVertAct){
	int ii,jj;
	double *coordList=NULL;
	int *listIndVert=NULL;
	
	listIndVert=(int*)calloc(nVertGrid,sizeof(int));
	coordList=(double*)calloc(dim()*nEdge,sizeof(double));
	
	
	for(ii=0;ii<nVertAct;ii++){
		listIndVert[ii]=vertstructAct[ii].index;
	}
	FindObjNum(ordIndVert,listIndVert,ordPosVert,nEdge,nVertAct);
	
	for(ii=0;ii<nEdge;ii++){
		for(jj=0;jj<dim();jj++){
			coordList[ii*dim()+jj]=vertstructAct[ordPosVert[ii]].coord[jj];
		}
	}
	free(listIndVert);
	return(coordList);
}

double* ExtractAngles(double *vecTest, int nVec, int nDim){
	
	int ii,jj,signSin;
	double *angles=NULL, *dist=NULL;
	
	if(nDim!=2){perror("Angles do not support more than 2 dim()ensions");exit(EXIT_FAILURE);}
	
	angles=(double*)calloc(nVec,sizeof(double));
	dist=(double*)calloc(nVec,sizeof(double));
	//printf("\n angles: presign Postsign ");
	for (ii=0; ii<nVec;ii++){
		for(jj=0;jj<nDim;jj++){
			dist[ii]=dist[ii]+pow(vecTest[nDim*ii+jj],2);
		}
		dist[ii]=sqrt(dist[ii]);
		angles[ii]=acos(vecTest[nDim*ii]/dist[ii]);
		//printf("     %lf",angles[ii]);
		signSin=CompareDouble(vecTest[nDim*ii+1],0.0);
		signSin=(1-abs(signSin))+signSin; 
		angles[ii]=((double)signSin)*angles[ii];
		//printf("   %lf  \n",angles[ii]);
	}
	
	free(dist);
	return(angles);
}

int IsClockWiseChain(int *minVertPosRet,int nEdge, int *ordIndVert, int *ordPosVert, vertexTemplate *vertstructAct, int nVertAct){
	int ii,jj;
	int minVertPos,cwPos;
	double *coordList=NULL, *angleTest=NULL, vecTest[2*dim()];
	int vertPosTest[2];
	
	coordList=(double*)ExtractCoord(nEdge, ordIndVert,ordPosVert, vertstructAct,nVertAct);
	minVertPos=LowerLeftVertex(coordList,nEdge);
	
	vertPosTest[0]=minVertPos-1;
	vertPosTest[1]=minVertPos+1;
	if(minVertPos==0){vertPosTest[0]=nEdge-1;}
	if(minVertPos==(nEdge-1)){vertPosTest[1]=0;}
	
	for(ii=0;ii<2;ii++){
		for(jj=0;jj<dim();jj++){
			vecTest[dim()*ii+jj]=coordList[vertPosTest[ii]*dim()+jj]-coordList[minVertPos*dim()+jj];
		}
	}
	angleTest=(double*)ExtractAngles( vecTest,2,dim());
	cwPos=pos_fmax_array(angleTest,2);
	//printf("\ncw? %i angles: %lf %lf \n",cwPos,angleTest[0],angleTest[1]);
	//printf("vertprec %i vert next %i \n",ordIndVert[vertPosTest[0]],ordIndVert[vertPosTest[1]]);
	
	free(coordList);
	free(angleTest);
	*minVertPosRet=minVertPos;
	return(cwPos);
}

void ReorderArray(void *oldArray, int *newOrder ,int nArray, int sizArray){
	
	int ii;
	void *tempArray=NULL;
	tempArray=(void*)calloc(nArray,sizArray);
	
	for(ii=0;ii<nArray;ii++){
		memcpy(((char*)tempArray+ii*sizArray),((char*)oldArray+(newOrder[ii])*sizArray),sizArray);
	}
	memcpy(oldArray,tempArray,sizArray*nArray);
	
	free(tempArray);
}

void OrderEdgeChain(int nEdge, int *ordPosEdge, int *ordIndEdge, int *ordIndVert,int *ordPosVert,
	vertexTemplate *vertstructAct, int nVertAct){
	
	int ii,jj;
	int minVertPos,cwPos;
	int *newOrder;
	
	newOrder=(int*)calloc(nEdge,sizeof(int));
	
	cwPos=IsClockWiseChain(&minVertPos,nEdge, ordIndVert, ordPosVert, vertstructAct, nVertAct);
	//printf("\n lower left corner detected as %i isCW %i\n",ordIndVert[minVertPos], cwPos);
	//printf("\n  ");
	/*
	for(ii=0;ii<nEdge;ii++){
		
		printf("%i ",ordIndVert[ii]);
	}
	printf("\n");
	for(ii=0;ii<nEdge;ii++){
		
		printf("%i ",ordIndEdge[ii]);
	}
	printf("\n");
	*/
	for(ii=0;ii<nEdge;ii++){
		newOrder[ii]=(nEdge+minVertPos+cwPos*ii-(1-cwPos)*ii)%nEdge;
		//printf("neworder: %i orderedPos: %i \n",newOrder[ii],ordIndVert[ii]);
	}
	ReorderArray(ordIndVert, newOrder ,nEdge, sizeof(int));
	ReorderArray(ordPosVert, newOrder ,nEdge, sizeof(int));
	/*
	for(ii=0;ii<nEdge;ii++){
		newOrder[ii]=(nEdge+minVertPos+cwPos*ii-(1-cwPos)*(ii+1))%nEdge;
		printf(" %i ",newOrder[ii]);
	}*/
	ReorderArray(ordPosEdge, newOrder ,nEdge, sizeof(int));
	ReorderArray(ordIndEdge, newOrder ,nEdge, sizeof(int));
	/*
	printf("\n  ");
	for(ii=0;ii<nEdge;ii++){
		printf("neworder: %i orderedPos: %i %i\n",newOrder[ii],ordIndVert[ii],ordIndEdge[ii]);
	}
	*/
	free(newOrder);
}

void GenerateIndMatch(int domSize[dim()], int posCellRefine ,int *posEdgeSideTemp,int *posVertSideTemp,int *posCellAddTemp,
	int *posEdgeAddTemp, int *posVertAddTemp, int *ordPosEdge, int *ordPosVert, int nCurrEdge,
	int **convCellList, int **convEdgeList, int **convVertList){
	
	//int *convCellList=NULL,*convEdgeList=NULL,*convVertList=NULL;
	int *tempCellList=NULL,*tempEdgeList=NULL,*tempVertList=NULL;
	int *gridCellList=NULL,*gridEdgeList=NULL,*gridVertList=NULL;
	
	int ii;
	int maxCellTemplate,maxEdgeTemplate,maxVertTemplate;
	int maxCellGrid,maxEdgeGrid,maxVertGrid;
	int nCellTemplate,nEdgeTemplate,nVertTemplate;
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	
	tempCellList=(int*)calloc(nCellTemplate,sizeof(int));
	tempEdgeList=(int*)calloc(nEdgeTemplate,sizeof(int));
	tempVertList=(int*)calloc(nVertTemplate,sizeof(int));
	for (ii=0;ii<nCellTemplate;ii++){tempCellList[ii]=cellCurrentTemplate[ii].index;}
	for (ii=0;ii<nEdgeTemplate;ii++){tempEdgeList[ii]=edgeCurrentTemplate[ii].index;}
	for (ii=0;ii<nVertTemplate;ii++){tempVertList[ii]=vertCurrentTemplate[ii].index;}
	maxCellTemplate=max_array(tempCellList,nCellTemplate);
	maxEdgeTemplate=max_array(tempEdgeList,nEdgeTemplate);
	maxVertTemplate=max_array(tempVertList,nVertTemplate);
	//printf("\n cell: %i edge: %i vert: %i",maxCellTemplate,maxEdgeTemplate,maxVertTemplate);
	
	gridCellList=(int*)calloc(nCellGrid,sizeof(int));
	gridEdgeList=(int*)calloc(nEdgeGrid,sizeof(int));
	gridVertList=(int*)calloc(nVertGrid,sizeof(int));
	for (ii=0;ii<nCellGrid;ii++){gridCellList[ii]=cellstruct[ii].index;}
	for (ii=0;ii<nEdgeGrid;ii++){gridEdgeList[ii]=edgestruct[ii].index;}
	for (ii=0;ii<nVertGrid;ii++){gridVertList[ii]=vertstruct[ii].index;}
	maxCellGrid=max_array(gridCellList,nCellGrid);
	maxEdgeGrid=max_array(gridEdgeList,nEdgeGrid);
	maxVertGrid=max_array(gridVertList,nVertGrid);
	
	// Generates matching lists
	*convCellList=(int*)calloc((maxCellTemplate+1),sizeof(int));
	*convEdgeList=(int*)calloc((maxEdgeTemplate+1),sizeof(int));
	*convVertList=(int*)calloc((maxVertTemplate+1),sizeof(int));
	for(ii=0;ii<nCurrEdge;ii++){
		(*convEdgeList)[edgeCurrentTemplate[posEdgeSideTemp[ii]].index]=edgestruct[ordPosEdge[ii]].index;
		(*convVertList)[vertCurrentTemplate[posVertSideTemp[ii]].index]=vertstruct[ordPosVert[ii]].index;
	}
	(*convCellList)[cellCurrentTemplate[0].index]=cellstruct[posCellRefine].index;
	for(ii=0;ii<(nCellTemplate-1);ii++){(*convCellList)[cellCurrentTemplate[ii+1].index]=maxCellGrid+1+ii;}
	for(ii=0;ii<(nEdgeTemplate-nCurrEdge);ii++){
		(*convEdgeList)[edgeCurrentTemplate[posEdgeAddTemp[ii]].index]=maxEdgeGrid+1+ii;
	}
	for(ii=0;ii<(nVertTemplate-nCurrEdge);ii++){
		(*convVertList)[vertCurrentTemplate[posVertAddTemp[ii]].index]=maxVertGrid+1+ii;
	}
	
	free(tempCellList);
	free(tempEdgeList);
	free(tempVertList);
	free(gridCellList);
	free(gridEdgeList);
	free(gridVertList);
	
}

void ModifBorderEdges(int nCurrEdge, int posCellRefine, int *ordPosEdge, int *posEdgeSideTemp, int *convCellList){
	
	int ii,jj,kk;
	int posCellInd;
	posCellInd=cellstruct[posCellRefine].index;
	
	for(ii=0;ii<nCurrEdge;ii++){
		jj=0;
		while((jj<2)&(edgestruct[ordPosEdge[ii]].cellind[jj]!=posCellInd)){
			jj++;
		}
		kk=0;
		while((kk<2)&(edgeCurrentTemplate[posEdgeSideTemp[ii]].cellind[kk]<0)){
			kk++;
		}
		if((kk==2) | (jj==2)){perror("Cell Index Matching Failed");exit(EXIT_FAILURE);}
		edgestruct[ordPosEdge[ii]].cellind[jj]=
			convCellList[edgeCurrentTemplate[posEdgeSideTemp[ii]].cellind[kk]];
	}
	
}

void ExtendGridStructures(int nCurrEdge,int nCellTemplate,int nEdgeTemplate,int nVertTemplate){
	
	int ii;
	
	vertstruct=(vertexTemplate*)realloc(vertstruct,(nVertGrid+nVertTemplate-nCurrEdge)*sizeof(vertexTemplate));
	edgestruct=(edgeTemplate*)realloc(edgestruct,(nEdgeGrid+nEdgeTemplate-nCurrEdge) * sizeof(edgeTemplate));
	cellstruct=(cellTemplate*)realloc(cellstruct,(nCellGrid+nCellTemplate-1) * sizeof(cellTemplate));
	
	for(ii=nCellGrid;ii<(nCellGrid+nCellTemplate-1);ii++){
		cellstruct[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
	}
	
}

void CopyEdgeStruct(int nCurrEdge,int domSize[dim()],int *posEdgeAddTemp,
	int *convEdgeList,int *convCellList,int *convVertList){
	
	int nCellTemplate,nEdgeTemplate,nVertTemplate;
	int ii,jj,kk,ll;
	int currPosGrid;
	
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	
	//printf("\n");
	for(ii=0;ii<(nEdgeTemplate-nCurrEdge);ii++){
		currPosGrid=nEdgeGrid+ii;
		memcpy(&(edgestruct[currPosGrid].index),&(edgeCurrentTemplate[posEdgeAddTemp[ii]].index),
			sizeof(edgeTemplate));
		
		edgestruct[currPosGrid].index=convEdgeList[edgestruct[currPosGrid].index];
		
		edgestruct[currPosGrid].cellind[0]=convCellList[edgestruct[currPosGrid].cellind[0]];
		edgestruct[currPosGrid].cellind[1]=convCellList[edgestruct[currPosGrid].cellind[1]];
		
		edgestruct[currPosGrid].vertex[0]=convVertList[edgestruct[currPosGrid].vertex[0]];
		edgestruct[currPosGrid].vertex[1]=convVertList[edgestruct[currPosGrid].vertex[1]];
		//printf("%i %i %i\n",edgestruct[currPosGrid].index,edgeCurrentTemplate[posEdgeAddTemp[ii]].orientation,edgestruct[currPosGrid].orientation);
	}
	
	
}

void CopyVertStruct(int nCurrEdge,int domSize[dim()],int *posVertAddTemp,
	int *convVertList, int *ordPosVert){
	
	int nCellTemplate,nEdgeTemplate,nVertTemplate;
	int ii,jj,kk,ll;
	int currPosGrid;
	double minDim[dim()], maxDim[dim()];
	
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	
	
	//printf("\n min Max");
	for(jj=0;jj<dim();jj++){
		minDim[jj]=vertstruct[ordPosVert[0]].coord[jj];
		maxDim[jj]=vertstruct[ordPosVert[0]].coord[jj];
		for(ii=1;ii<nCurrEdge;ii++){
			minDim[jj]=min(minDim[jj],vertstruct[ordPosVert[ii]].coord[jj]);
			maxDim[jj]=max(maxDim[jj],vertstruct[ordPosVert[ii]].coord[jj]);
		}
		//printf("\n %lf %lf",minDim[jj],maxDim[jj]);
	}
	
	for(ii=0;ii<(nVertTemplate-nCurrEdge);ii++){
		currPosGrid=nVertGrid+ii;
		vertstruct[currPosGrid].index=convVertList[vertCurrentTemplate[posVertAddTemp[ii]].index];
		for(jj=0;jj<dim();jj++){
			vertstruct[currPosGrid].coord[jj]=minDim[jj]+
				vertCurrentTemplate[posVertAddTemp[ii]].coord[jj]*(maxDim[jj]-minDim[jj]);
		}
	}
	
	
}

void CopyCellStruct(int domSize[dim()],int posCellRefine,int *convCellList){
	
	int nCellTemplate,nEdgeTemplate,nVertTemplate;
	
	int ii,jj,kk,ll;
	int *currPosGrid;
	int *refineVecOrig;
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	
	
	currPosGrid=(int*)calloc(nCellTemplate,sizeof(int));
	refineVecOrig=(int*)calloc(nLevels,sizeof(int));
	
	memcpy(refineVecOrig,cellstruct[posCellRefine].refineVec,sizeof(int)*nLevels);
	
	currPosGrid[0]=posCellRefine;
	for (ii=1;ii<nCellTemplate;ii++){currPosGrid[ii]=nCellGrid+ii-1;}
	//printf("\n cell assign");
	for(ii=0;ii<nCellTemplate;ii++){
		cellstruct[currPosGrid[ii]].index=convCellList[cellCurrentTemplate[ii].index];
		//printf("\n%i %i ",cellCurrentTemplate[ii].index,convCellList[cellCurrentTemplate[ii].index]);
		cellstruct[currPosGrid[ii]].fill=0;
		
		cellstruct[currPosGrid[ii]].refineLvl=cellCurrentTemplate[ii].refineLvl;
		for(jj=0;jj<nLevels;jj++){
			cellstruct[currPosGrid[ii]].refineVec[jj]=cellCurrentTemplate[ii].refineVec[jj]+refineVecOrig[jj];
		}
	}
	
	free(currPosGrid);
	free(refineVecOrig);
}

void MergeNewCell(int domSize[dim()], int posCellRefine ,int *posEdgeSideTemp,int *posVertSideTemp,int *posCellAddTemp,
	int *posEdgeAddTemp, int *posVertAddTemp, int *ordPosEdge, int *ordPosVert, int nCurrEdge){
	
	int ii;
	int *convCellList=NULL,*convEdgeList=NULL,*convVertList=NULL;
	int nCellTemplate,nEdgeTemplate,nVertTemplate;
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	
	
	GenerateIndMatch(domSize,posCellRefine, posEdgeSideTemp, posVertSideTemp, posCellAddTemp,
		posEdgeAddTemp,  posVertAddTemp,  ordPosEdge,  ordPosVert,nCurrEdge,
		&convCellList,&convEdgeList,&convVertList);
	/*
	printf("\n convCellList: ");
	for(ii=0;ii<(nCellTemplate+1);ii++){
		printf(" %i ",convCellList[ii]);
	}
	printf("\n convEdgeList: ");
	for(ii=0;ii<(nEdgeTemplate+1);ii++){
		printf(" %i ",convEdgeList[ii]);
	}
	printf("\n convVertList: ");
	for(ii=0;ii<(nVertTemplate+1);ii++){
		printf(" %i ",convVertList[ii]);
	}
	*/
	
	ModifBorderEdges(nCurrEdge,posCellRefine,ordPosEdge,posEdgeSideTemp, convCellList);
	ExtendGridStructures(nCurrEdge,nCellTemplate,nEdgeTemplate,nVertTemplate);
	CopyEdgeStruct(nCurrEdge,domSize,posEdgeAddTemp,convEdgeList,convCellList,convVertList);
	CopyCellStruct(domSize,posCellRefine,convCellList);
	CopyVertStruct(nCurrEdge,domSize,posVertAddTemp,convVertList,ordPosVert);
	
	nEdgeGrid=nEdgeGrid+(nEdgeTemplate-nCurrEdge);
	nCellGrid=nCellGrid+(nCellTemplate-1);
	nVertGrid=nVertGrid+(nVertTemplate-nCurrEdge);
	
	free(convCellList);
	free(convEdgeList);
	free(convVertList);
}

void RefineCell(int domSize[dim()],int posCellRefine,int indCellRefine,int *posEdgeSideTemp,
	int *posVertSideTemp,int *posCellAddTemp,int *posEdgeAddTemp, int *posVertAddTemp){
	
	int ii;
	int *posCurrEdge=NULL,*indCurrEdge=NULL;
	int *ordPosEdge=NULL, *ordIndEdge=NULL, *ordIndVert=NULL,*ordPosVert=NULL;
	int nCurrEdge;
	
	IdentifyRefineEdge(&posCellRefine, &indCellRefine,1,
			&posCurrEdge, &indCurrEdge,&nCurrEdge,edgestruct,nEdgeGrid);
			
	//printf("\n***Number of Edge of Cell: %i",nCurrEdge);
	ordPosEdge=(int*)calloc(nCurrEdge,sizeof(int));
	ordIndEdge=(int*)calloc(nCurrEdge,sizeof(int));
	ordIndVert=(int*)calloc(nCurrEdge,sizeof(int));
	ordPosVert=(int*)calloc(nCurrEdge,sizeof(int));		
			
	ExtractEdgeChain(posCurrEdge, indCurrEdge,nCurrEdge ,ordPosEdge, ordIndEdge, ordIndVert,edgestruct);
	OrderEdgeChain(nCurrEdge, ordPosEdge, ordIndEdge, ordIndVert,ordPosVert, vertstruct,nVertGrid);
	/*printf("\nRefine Cell edge index list:\n");
	for (ii=0;ii<nCurrEdge;ii++){printf("%i ",ordIndEdge[ii]);}
	printf("\n");*/
	
	MergeNewCell(domSize,posCellRefine, posEdgeSideTemp, posVertSideTemp, posCellAddTemp,
		posEdgeAddTemp,  posVertAddTemp,  ordPosEdge,  ordPosVert,nCurrEdge);
	//printf("\n");
	//for (ii=0;ii<nCurrEdge;ii++)
	//	printf("nOrdered: %i  ; Ordered Vert: %i ; Ordered Edge: %i \n",indCurrEdge[ii],ordIndEdge[ii],ordIndVert[ii]);
	
	free(posCurrEdge);
	free(indCurrEdge);
	free(ordPosEdge);
	free(ordIndEdge);
	free(ordIndVert);
	free(ordPosVert);
}

int* OppositeList(int *fullList, int *partial, int nFull, int nPart){
	int ii,jj,kk,test;
	int *oppList=NULL;
	
	oppList=(int*)calloc(nFull-nPart+1,sizeof(int));
	kk=0;
	for(ii=0;ii<nFull;ii++){
		jj=0;
		do{
			test=fullList[ii]!=partial[jj];
			jj++;
		} while (test & (jj<nPart));
		oppList[kk]=fullList[ii];
		kk=kk+((jj==nPart) & test);
	}
	oppList=(int*)realloc(oppList,(nFull-nPart)*sizeof(int));
	
	if (kk!=(nFull-nPart)){
		printf("OppositeList Not Working as expected\n");
		printf("%i %i %i",kk, nFull, nPart);
		exit(EXIT_FAILURE);
	}
	
	return(oppList);
}

void PrepareTemplateInfo(int domSize[dim()],int **posEdgeSide,int **posVertSide,int **posCellAdd,
	int **posEdgeAdd,int **posVertAdd){
	
	int ii;
	int *posCurrEdge=NULL,*indCurrEdge=NULL, indCellRefine[4]={-1,-2,-3,-4};
	int *ordIndEdge=NULL, *ordIndVert=NULL, *listTempEdge=NULL,*listTempVert=NULL;
	int nCurrEdge;
	
	int nCellTemplate,nEdgeTemplate,nVertTemplate;
	CalculateNumElements(domSize,&nCellTemplate,&nEdgeTemplate,&nVertTemplate);
	// Cell
	*posCellAdd=(int*)calloc(nCellTemplate,sizeof(int));
	for (ii=0;ii<nCellTemplate;ii++){ (*posCellAdd)[ii] = ii; }
	// Edge and Vertex Sides
	IdentifyRefineEdge(indCellRefine, indCellRefine,4,
			&posCurrEdge, &indCurrEdge,&nCurrEdge,edgeCurrentTemplate,nEdgeTemplate);
			
	//printf("\n***Number of Border Edges for template %i",nCurrEdge);
	*posEdgeSide=(int*)calloc(nCurrEdge,sizeof(int));
	ordIndEdge=(int*)calloc(nCurrEdge,sizeof(int));
	ordIndVert=(int*)calloc(nCurrEdge,sizeof(int));
	*posVertSide=(int*)calloc(nCurrEdge,sizeof(int));		
	
	ExtractEdgeChain(posCurrEdge, indCurrEdge,nCurrEdge ,*posEdgeSide, ordIndEdge, ordIndVert,edgeCurrentTemplate);
	OrderEdgeChain(nCurrEdge, *posEdgeSide, ordIndEdge, ordIndVert,*posVertSide, vertCurrentTemplate,nVertTemplate);
	/*printf("\nPrepare template order edge index list:\n");
	for (ii=0;ii<nCurrEdge;ii++){printf("%i ",ordIndEdge[ii]);}
	printf("\n");*/
	// Inside 
	listTempEdge=(int*)calloc(nEdgeTemplate,sizeof(int));
	listTempVert=(int*)calloc(nVertTemplate,sizeof(int));
	for (ii=0;ii<nEdgeTemplate;ii++){listTempEdge[ii]=ii;}
	for (ii=0;ii<nVertTemplate;ii++){listTempVert[ii]=ii;}
	*posEdgeAdd=(int*)OppositeList(listTempEdge, *posEdgeSide, nEdgeTemplate, nCurrEdge);
	*posVertAdd=(int*)OppositeList(listTempVert, *posVertSide, nVertTemplate, nCurrEdge);
	
	free(posCurrEdge);
	free(indCurrEdge);
	//free(ordPosEdge);
	free(ordIndEdge);
	free(ordIndVert);
	free(listTempEdge);
	free(listTempVert);
}

void RefineSelectedCells(int domSize[dim()],int *posCellRefine,int *indCellRefine,int nCellRefine){
	
	int ii,jj,kk,ll;
	int *posEdgeSideTemp=NULL, *posVertSideTemp=NULL, *posCellAddTemp=NULL;
	int *posEdgeAddTemp=NULL, *posVertAddTemp=NULL;
	PrepareTemplateInfo(domSize,&posEdgeSideTemp,&posVertSideTemp,&posCellAddTemp,
			&posEdgeAddTemp, &posVertAddTemp);
	
	for(ii=0;ii<nCellRefine;ii++){
		 RefineCell(domSize,posCellRefine[ii],indCellRefine[ii],posEdgeSideTemp,
		 posVertSideTemp,posCellAddTemp,posEdgeAddTemp, posVertAddTemp);
	}
	
	
}


// UTILITIES
	//Various

void* ArrayFromStruct(void* ptr, int nStruct, int sizMember, int sizType){
	int ii;
	void *dest;
	dest=calloc(nStruct,sizType);
	
	for(ii=0;ii<nStruct;ii++){
		memcpy(((char*)dest+ii*sizType),((char*)ptr+(ii*sizMember)),sizType);
	}
	return(dest);
}

void RemoveIdenticalEntries_int(int **array, int nArray, int *nNewArray){
	
	int jj,ii;
	
	qsort((*array),nArray,sizeof(int),SortOrder);
	jj=1;
	//printf("\n");
	for (ii=1;ii<nArray;ii++){
		if ((*array)[ii]>(*array)[jj-1]){
			(*array)[jj]=(*array)[ii];
			jj++;
			//printf("%i ", jj);
		}
	}
	(*array)=(int*)realloc((*array),jj*sizeof(int));
	*nNewArray=jj;
}

void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList){
	// Look up is the Indices we need the position of
	// listInd is the list of all indices in the right order
	// dest is where the position will be stored (-1 means the index wasn't found)
	// While not the most efficient yet this should be pretty good
	int ii,jj;
	
	for(ii=0;ii<nLook;ii++){
		dest[ii]=-1;
		jj=0;
		do{
			
			dest[ii]=-abs(listInd[jj]-lookup[ii])*nList+jj;
			jj++;
		} while((dest[ii]!=(jj-1)) & (jj<nList));
	}

}

int SortOrder(const void * elem1, const void * elem2) {
    int f = *((int*)elem1);
    int s = *((int*)elem2);
	
    return (f > s) - (f < s);
}

int CompareDouble(double elem1, double elem2) {
    double f = elem1;
    double s = elem2;
	
    return (f > s) - (f < s);
}

	// Index functions
int edgsub(int I, int J, int l, int domSize[dim()]){
	/*
	// l is either 0 for dim() 1 or 1 for dim() 2
	// domSize is the actual requested domain size rather than the edge domain size
	// I and J vary normallly from 1:domSize
	// J can go to domSize[1]+1 with l=0
	// I can go to domSize[0]+1 with l=1 
	*/
	
	int index;
	index=-1 // Use of a mathematical switch
			+(1-l)*(I+(J-1)*domSize[0]) // case where l=0
			+(l)*((I-1)*(domSize[1])+(J) // case where l=1
				+((domSize[1]+1)*domSize[0]));
	return(index);
}

int vertsub(int I, int J, int domSize[dim()]){
	/*
	// domSize is the actual requested domain size rather than the edge domain size
	// I and J vary normallly from 1:domSize+1
	*/
	
	int index;
	index=-1+(I+(J-1)*(domSize[0]+1));
	//printf("domsize %i %i , I %i J %i index %i",domSize[0],domSize[1],I,J,index);
	return(index);
}

int cellsub(int I, int J, int domSize[dim()]){
	/*
	// domSize is the actual requested domain size rather than the edge domain size
	// I and J vary normallly from 1:domSize
	*/
	
	int index;
	index=-1+(I+(J-1)*(domSize[0]));
	return(index);
}

	// Min Max functions 
	// (see Macro)
int min_array(int *array, int nArray){
	int ii,minNum;
	minNum=array[0];
	for (ii=1;ii<nArray;ii++){
		minNum=min(minNum,array[ii]);
	}
	return(minNum);
}

int pos_min_array(int *array, int nArray){
	int ii,minNum,minPrev,minPos;
	minNum=array[0];
	minPos=0;
	for (ii=1;ii<nArray;ii++){
		minPrev=minNum;
		minNum=min(minNum,array[ii]);
		minPos=minPos*(minNum==minPrev)+ii*(minPrev!=minNum);
	}
	return(minPos);
}

int max_array(int *array, int nArray){
	int ii,minNum;
	minNum=array[0];
	for (ii=1;ii<nArray;ii++){
		minNum=max(minNum,array[ii]);
	}
	return(minNum);
}

int pos_max_array(int *array, int nArray){
	int ii,minNum,minPrev,minPos;
	minNum=array[0];
	minPos=0;
	for (ii=1;ii<nArray;ii++){
		minPrev=minNum;
		minNum=max(minNum,array[ii]);
		minPos=minPos*(minNum==minPrev)+ii*(minPrev!=minNum);
	}
	return(minPos);
}

double fmin_array(double *array, int nArray){
	int ii;
	double minNum;
	minNum=array[0];
	for (ii=1;ii<nArray;ii++){
		minNum=min(minNum,array[ii]);
	}
	return(minNum);
}

int pos_fmin_array(double *array, int nArray){
	int ii,minPos;
	double minNum,minPrev;
	minNum=array[0];
	minPos=0;
	for (ii=1;ii<nArray;ii++){
		minPrev=minNum;
		minNum=min(minNum,array[ii]);
		minPos=minPos*(minNum==minPrev)+ii*(minPrev!=minNum);
	}
	return(minPos);
}

double fmax_array(double *array, int nArray){
	int ii;
	double minNum;
	minNum=array[0];
	for (ii=1;ii<nArray;ii++){
		minNum=max(minNum,array[ii]);
	}
	return(minNum);
}

int pos_fmax_array(double *array, int nArray){
	int ii,minPos;
	double minNum,minPrev;
	minNum=array[0];
	minPos=0;
	for (ii=1;ii<nArray;ii++){
		minPrev=minNum;
		minNum=max(minNum,array[ii]);
		minPos=minPos*(minNum==minPrev)+ii*(minPrev!=minNum);
	}
	return(minPos);
}

int imin(int a, int b){return (((a>b)*b)+((a<b)*a));}
int imax(int a, int b){return (((a<b)*b)+((a>b)*a));}
double f_min(double a, double b){return (((a>b)*b)+((a<b)*a));}
double f_max(double a, double b){return (((double)(a<b)*b)+((double)(a>b)*a));}


