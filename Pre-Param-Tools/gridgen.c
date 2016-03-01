#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
// Type Definitions
typedef struct {
	int index;
	double fill;
	int refineLvl;
	int *refineVec; // This has to be allocated 
	
} cellTemplate;

typedef struct {
	int index;
	int cellind[2];
	int vertex[2];
	
} edgeTemplate;

typedef struct {
	int index;
	double coord[dim];
	
} vertexTemplate;


// Constant declaration
int dim=2;


// Global Variables declaration
int nLevels, nCells,nEdges,nVerts;
int *levelSize,*cells;
cellTemplate *celldatstruct; // Array containing data from File
cellTemplate *cellstruct;
vertTemplate *vertstruct;
edgeTemplate *edgestruct;


// Function Declarations



// Allocation functions
void Allocatecelldatstruct(){
	int ii;
	
	celldatstruct=(cellTemplate*)malloc(nCells * sizeof(cellTemplate));
	for(ii=0;ii<nCells;ii++){
		celldatstruct[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
	} 
}

void AllocateGridStruct(int domSize[dim]){
	
	int nCellCurr,nEdgeCurr,nVertCurr;
	nCellCurr=domSize[0]*domSize[1];
	nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	nVertCurr=(domSize[0]+1)*(domSize[1]+1);
	
	
	vertstruct=(vertTemplate*)malloc(nVertCurr * sizeof(vertTemplate));
	edgestruct=(edgeTemplate*)malloc(nEdgeCurr * sizeof(edgeTemplate));
	
	cellstruct=(cellTemplate*)malloc(nCellCurr * sizeof(cellTemplate));
	for(ii=0;ii<nCellCurr;ii++){
		celldatstruct[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
	} 
}

// Actions

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
				fscanf(cellgridFID,"%i ",&levelSize[dim*ii]);
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

void EdgeIJtoGrid(int domSize[dim],int IJK[dim], int l){
	
	index=IJK[0]+(IJK[1]-1)*domSize[1]; // l=1
	index=domSize[0]*(domSize[1]+1)+IJK[1]+(IJK[0]-1)*domSize[0]; // ll=2
	
	vertexind[0]=IJK[0]+(IJK[1]-1)*(domSize[1]+1);
	vertexind[1]=1+IJK[0]+(IJK[1]-1)*(domSize[1]+1);// l=1
	
	vertexind[0]=IJK[0]+(IJK[1]-1)*(domSize[1]+1);
	vertexind[1]=IJK[0]+(IJK[1]-1+1)*(domSize[1]+1);// ll=2
	
	cellind[0]=IJK[0]+(IJK[1]-2)*domSize[1];
	cellind[1]=IJK[0]+(IJK[1]-1)*domSize[1]; //ll=1

	cellind[0]=-1+IJK[0]+(IJK[1]-1)*domSize[1];
	cellind[1]=IJK[0]+(IJK[1]-1)*domSize[1]; //ll=2
	
}

void AssignEdgeStructContent(int domSize[dim]){
	int ii,jj,IJK[dim];
	
	
	for (ii=1;ii<=domSize[0];ii++){
		for (jj=1;jj<=domSize[0];jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			EdgeIJtoGrid(domSize[dim],IJK[dim],1);
			EdgeIJtoGrid(domSize[dim],IJK[dim],2);
		}
	}
	
}


// Main text body
int main(){
	 //Arrays
	DataIn();
	Allocatecelldatstruct();
	DataArrayToStruct();
	
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
	return(0);
}

