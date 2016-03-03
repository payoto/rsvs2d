#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "C:\Users\ap1949\Local Documents\PhD\Development Work\Snakes Volume Parametrisation\source\Pre-Param-Tools\templategrid\templategrid.h"

// Constant declaration
// dim=2;

// Global Variables declaration
int nLevels, nCells,nEdges,nVerts;
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
void GenerateTemplates();
void DeAllocate();

// Main text body
int main(){
	 //Arrays
	 printf("%i\n",sizeof(double));
	DataIn();
	Allocatecelldatstruct();
	DataArrayToStruct();
	GenerateTemplates();
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
	nCellCurr=domSize[0]*domSize[1];
	nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	nVertCurr=(domSize[0]+1)*(domSize[1]+1);
	if (lvlGrid>=10){
		fileName=(char*)malloc((15+(int)(log10(lvlGrid))+2)*sizeof(char));
	} else if ((lvlGrid<=1) & (lvlGrid<10)) {
		fileName=(char*)malloc((15+(int)(log10(lvlGrid))+2)*sizeof(char));
	} else {
		fileName=(char*)malloc((15)*sizeof(char));
	}
	printf("checkTemplate%i.m\n",lvlGrid);
	sprintf(fileName,"checkTemplate%i.m",lvlGrid);
	checkFID=fopen(fileName,"w");
	
	for (ii=0;ii<=nEdgeCurr-1;ii++){
		fprintf(checkFID,"templateGrid.edge(%i).index=%i;\n",ii+1,edgeCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(1)=%i;\n",ii+1,edgeCurrentTemplate[ii].cellind[0]);
		fprintf(checkFID,"templateGrid.edge(%i).cellindex(2)=%i;\n",ii+1,edgeCurrentTemplate[ii].cellind[1]);
		
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(1)=%i;\n",ii+1,edgeCurrentTemplate[ii].vertex[0]);
		fprintf(checkFID,"templateGrid.edge(%i).vertexindex(2)=%i;\n",ii+1,edgeCurrentTemplate[ii].vertex[1]);
	}
	
	for (ii=0;ii<=nCellCurr-1;ii++){
		fprintf(checkFID,"templateGrid.cell(%i).index=%i;\n",ii+1,cellCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.cell(%i).fill=%lf;\n",ii+1,cellCurrentTemplate[ii].fill);
	}
	
	for (ii=0;ii<=nVertCurr-1;ii++){
		fprintf(checkFID,"templateGrid.vertex(%i).index=%i;\n",ii+1,vertCurrentTemplate[ii].index);
		
		fprintf(checkFID,"templateGrid.vertex(%i).coord(1)=%lf;\n",ii+1,vertCurrentTemplate[ii].coord[0]);
		fprintf(checkFID,"templateGrid.vertex(%i).coord(2)=%lf;\n",ii+1,vertCurrentTemplate[ii].coord[1]);
		
	} 
	
	fclose(checkFID);
	
}

// Template grid main process
void GenerateTemplateGrid(int lvlGenerate){
	
	int domSize[dim],baseRefineLvl,ii;
	ii=lvlGenerate-1;
	domSize[0]=levelSize[2*ii];
	domSize[1]=levelSize[2*ii+1];
	baseRefineLvl=lvlGenerate;
	BuildLvlTemplate(domSize, baseRefineLvl,nLevels,&cellCurrentTemplate,&vertCurrentTemplate,&edgeCurrentTemplate);
	
		printf("\n*****  ADRESS OF &cellCurrentTemplate in function is: %d\n",&(cellCurrentTemplate[0].index));
	printf("\n*****  ADRESS OF cellCurrentTemplate[1] in function is: %d\n",cellCurrentTemplate);
	printf("\n*****  ADRESS OF *cellCurrentTemplate in function is: %d\n",*cellCurrentTemplate);
	OutputTemplateGrid(domSize,lvlGenerate);
	
}

void GenerateTemplates(){
	
	int domSize[dim],baseRefineLvl,ii;
	
	for (ii=1;ii<=nLevels;ii++){
		GenerateTemplateGrid(ii);
	}
}

void DeAllocate(){
	
	free(celldatstruct);
	free(cells);
	free(levelSize);
	free(cellCurrentTemplate);
	free(vertCurrentTemplate);
	free(edgeCurrentTemplate);
	free(cellstruct);
	free(vertstruct);
	free(edgestruct);

}
