#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\Pre-Param-Tools\templategrid\templategrid.h"

// Constant declaration
// dim=2;

// Global Variables declaration
int nLevels, nCells,nEdges,nVerts,nCellGrid,nEdgeGrid,nVertGrid;
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

// Function Declarations
void Allocatecelldatstruct();
void DataIn();
void DataArrayToStruct();
void OutputTemplateGrid(int domSize[dim], int lvlGrid);
void GenerateTemplateGrid(int lvlGenerate);
void GenerateGrid();
void DeAllocate();
void DeAllocateTemplate();
void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList);
void GridInitialisation();
void OutputGrid(int domSize[dim]);
int SortOrder(const void * elem1, const void * elem2);

// Main text body
int main(){
	 //Arrays
	 printf("%i\n",sizeof(double));
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
	
	 //DeAllocate();
	return(0);
}

// Allocation functions
void Allocatecelldatstruct(){
	int ii;
	
	celldatstruct=(cellTemplate*)malloc(nCells * sizeof(cellTemplate));
	for(ii=0;ii<nCells;ii++){
		celldatstruct[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
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
		levelSize=(int*)calloc(dim*nLevels, sizeof(int));
		for (ii=0;ii<=nLevels-1;ii++){ // Reads in the size of those levels
			for (jj=0; jj<dim;jj++){
				fscanf(cellgridFID,"%i ",&levelSize[dim*ii+jj]);
				// printf("%i %i\n",levelSize[2*ii],levelSize[2*ii+1]);
			}
		}
		
		fscanf(cellgridFID,"%i",&nCells); // Reads in the number of cells
		printf("	Reading cells in\n");
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
		printf("Cell Grid Structure succesfully read in\n\n");
		} else {
		perror("Data file failed to open!");
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

void OutputTemplateGrid(int domSize[dim], int lvlGrid){

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
	printf("checkTemplate%i.m... ",lvlGrid);
	for (ii=0;ii<=nEdgeCurr-1;ii++){
		fprintf(checkFID,"templateGrid.edge(%i).index=%i;\n",ii+1,edgeCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(1)=%i;\n",ii+1,edgeCurrentTemplate[ii].cellind[0]);
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(2)=%i;\n",ii+1,edgeCurrentTemplate[ii].cellind[1]);
		
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(1)=%i;\n",ii+1,edgeCurrentTemplate[ii].vertex[0]);
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(2)=%i;\n",ii+1,edgeCurrentTemplate[ii].vertex[1]);
	}
	printf(" . ");
	for (ii=0;ii<=nCellCurr-1;ii++){
		fprintf(checkFID,"templateGrid.cell(%i).index=%i;\n",ii+1,cellCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.cell(%i).fill=%lf;\n",ii+1,cellCurrentTemplate[ii].fill);
	}
	printf(" . ");
	for (ii=0;ii<=nVertCurr-1;ii++){
		fprintf(checkFID,"templateGrid.vertex(%i).index=%i;\n",ii+1,vertCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.vertex(%i).coord(1)=%lf;\n",ii+1,vertCurrentTemplate[ii].coord[0]);
		fprintf(checkFID,"templateGrid.vertex(%i).coord(2)=%lf;\n",ii+1,vertCurrentTemplate[ii].coord[1]);
		
	} 
	printf(" done!\n");
	fclose(checkFID);
	
}


void OutputGrid(int domSize[dim]){

	int nCellCurr,nEdgeCurr,nVertCurr, ii;
	FILE *checkFID;
	
	int *nCellP, *nEdgeP, *nVertP;
	nCellP=&nCellCurr;
	nEdgeP=&nEdgeCurr;
	nVertP=&nVertCurr;
	CalculateNumElements(domSize,nCellP,nEdgeP,nVertP);
	
	checkFID=fopen("checkGrid.m","w");
	
	for (ii=0;ii<=nEdgeCurr-1;ii++){
		fprintf(checkFID,"templateGrid.edge(%i).index=%i;\n",ii+1,edgestruct[ii].index);
		
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(1)=%i;\n",ii+1,edgestruct[ii].cellind[0]);
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(2)=%i;\n",ii+1,edgestruct[ii].cellind[1]);
		
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(1)=%i;\n",ii+1,edgestruct[ii].vertex[0]);
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(2)=%i;\n",ii+1,edgestruct[ii].vertex[1]);
	}
	printf(" . ");
	for (ii=0;ii<=nCellCurr-1;ii++){
		fprintf(checkFID,"templateGrid.cell(%i).index=%i;\n",ii+1,cellstruct[ii].index);
		
		fprintf(checkFID,"templateGrid.cell(%i).fill=%lf;\n",ii+1,cellstruct[ii].fill);
	}
	printf(" . ");
	for (ii=0;ii<=nVertCurr-1;ii++){
		fprintf(checkFID,"templateGrid.vertex(%i).index=%i;\n",ii+1,vertstruct[ii].index);
		
		fprintf(checkFID,"templateGrid.vertex(%i).coord(1)=%lf;\n",ii+1,vertstruct[ii].coord[0]);
		fprintf(checkFID,"templateGrid.vertex(%i).coord(2)=%lf;\n",ii+1,vertstruct[ii].coord[1]);
		
	} 
	printf(" done!\n");
	fclose(checkFID);
	
}


// Template grid main process
void GenerateTemplateGrid(int lvlGenerate){
	
	int domSize[dim],baseRefineLvl,ii;
	ii=lvlGenerate-1;
	domSize[0]=levelSize[2*ii];
	
	domSize[1]=levelSize[2*ii+1];
	printf("%i\n",domSize[1]);
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
	
	int domSize[dim];	
	int *nCellP, *nEdgeP, *nVertP;
	int nCellCurr, nEdgeCurr, nVertCurr;
	
	// First step initialises grid
	domSize[0]=levelSize[0];
	domSize[1]=levelSize[1];
	printf("%i\n",domSize[1]);
	AllocateGridStruct(domSize,&cellstruct,&vertstruct,&edgestruct);
	GenerateTemplateGrid(1);
	
	nCellP=&nCellCurr;
	nEdgeP=&nEdgeCurr;
	nVertP=&nVertCurr;
	CalculateNumElements(domSize,nCellP,nEdgeP,nVertP);
	
	memcpy(cellstruct,cellCurrentTemplate,sizeof(*cellCurrentTemplate)*(nCellCurr));
	memcpy(vertstruct,vertCurrentTemplate,sizeof(*vertCurrentTemplate)*(nVertCurr));
	memcpy(edgestruct,edgeCurrentTemplate,sizeof(*edgeCurrentTemplate)*(nEdgeCurr));
	nCellGrid=nCellCurr;
	nEdgeGrid=nEdgeCurr;
	nVertGrid=nVertCurr;
	
	printf("memory was copied\n");
	DeAllocateTemplate();
	printf("Template deallocated\n");
	OutputGrid(domSize);
	
}

void IdentifyRefineCell(int refinLvl){
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
	int *posGridRef=NULL, *indGridRef=NULL;
	// first count the number of refinements needed

	posRefV= (int *)calloc(nCells,sizeof(int));
	posCellRef= (int *)calloc(nCells*(refinLvl-1),sizeof(int));
	printf("%i \n",refinLvl-2);
	kk=0;
	for (ii=0;ii<nCells;ii++){
	
		printf("%i %i\n",ii,kk);
		if (celldatstruct[ii].refineLvl>=refinLvl){
			posRefV[kk]=ii;
			for (ll=0;ll<refinLvl-1;ll++){
				posCellRef[kk*(refinLvl-1)+ll]=celldatstruct[ii].refineVec[ll];
			}
			kk++;
		}
	}
	printf("%i",kk);
	posRefV=(int*)realloc(posRefV,kk*sizeof(int));
	posCellRef=(int*)realloc(posCellRef,kk*sizeof(int)*(refinLvl-1));
	
	nPosRefine=kk;
	posGridRef=(int*)calloc(nCellGrid,sizeof(int));
	kk=0;
	
	// Match vectors to cell location 
	// -> Note will be more efficient to remove 
	//		vectors which are the same and then do this step
	printf("\n");
	for (ii=0;ii<nPosRefine;ii++){
		for (jj=0;jj<nCellGrid;jj++){
			ll=0;
			while((cellstruct[jj].refineVec[ll]==posCellRef[ii*(refinLvl-1)+ll])
				& (ll<(refinLvl-1))){
				ll++;
				printf("%i ",ll);
			}
			if ((ll)==(refinLvl-1)){
				posGridRef[kk]=jj;
				kk++;
				printf("\n");
				break;
			}
		}
	}
	
	if (kk==0){printf("ERROR: The specified refineVec does not exist");}
	
	posGridRef=(int*)realloc(posGridRef,kk*sizeof(int));
	/*
	// Checks - Seem to work
	printf("\n");
	for (ii=0;ii<kk;ii++) {
		printf("%i - ",posGridRef[ii]);
		for (jj=0;jj<(refinLvl-1);jj++){
			printf("%i ",cellstruct[posGridRef[ii]].refineVec[jj]);
		}
		printf("\n");
	}
	
	
	//printf("\n");
	//for (ii=0;ii<kk;ii++) {printf("%i ",posCellRef[ii]);}
	//printf("\n");
	*/
	
	qsort(posGridRef,nPosRefine,sizeof(int),SortOrder);
	jj=1;
	printf("\n");
	for (ii=1;ii<nPosRefine;ii++){
		if (posGridRef[ii]>posGridRef[jj-1]){
			posGridRef[jj]=posGridRef[ii];
			jj++;
			//printf("%i ", jj);
		}
	}
	posGridRef=(int*)realloc(posGridRef,jj*sizeof(int));
	indGridRef=(int*)calloc(jj,sizeof(int));
	for(ii=0;ii<jj;ii++){
		indGridRef[ii]=cellstruct[posGridRef[ii]].index;
	}
	/*
	printf("\n");
	for (ii=0;ii<jj;ii++){
		printf("%i - ",posGridRef[ii]);
		for (ll=0;ll<(refinLvl-1);ll++){
			printf("%i ",cellstruct[posGridRef[ii]].refineVec[ll]);
		}
		printf("\n");
	}
	*/
	
	free(posRefV);
	free(posCellRef);
	free(posGridRef);
}

int SortOrder(const void * elem1, const void * elem2) {
    int f = *((int*)elem1);
    int s = *((int*)elem2);
	
    return (f > s) - (f < s);
}

void GenerateGrid(){
	int ii;
	int *posRefine;
	GridInitialisation();
	
	for (ii=2;ii<=nLevels;ii++){
		GenerateTemplateGrid(ii);
		
		IdentifyRefineCell(ii);
		
	}
}

void DeAllocateTemplate(){
	
	free(cellCurrentTemplate);
	free(vertCurrentTemplate);
	free(edgeCurrentTemplate);

}

void DeAllocate(){
	
	free(celldatstruct);
	free(cells);
	free(levelSize);
	free(cellstruct);
	free(vertstruct);
	free(edgestruct);

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
			jj+1;
			dest[ii]=-abs(listInd[jj]-lookup[ii])*nList+jj;
		} while((dest[ii]!=jj) & (jj<nList));
	}

}
