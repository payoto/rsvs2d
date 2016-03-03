#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "templategrid.h"


// Constant declaration
int dim=2;

// Global Variables declaration
int nLevels, nCells,nEdges,nVerts;
int *levelSize,*cells;
cellTemplate *celldatstruct; // Array containing data from File
// Array for active template
cellTemplate *cellstructTemp;
vertexTemplate *vertstructTemp;
edgeTemplate *edgestructTemp;


// Function Declarations



// Allocation functions

void AllocateGridStruct(int domSize[dim]){
	
	int nCellCurr,nEdgeCurr,nVertCurr, ii;
	nCellCurr=domSize[0]*domSize[1];
	nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	nVertCurr=(domSize[0]+1)*(domSize[1]+1);
	
	
	vertstructTemp=(vertexTemplate*)malloc(nVertCurr * sizeof(vertexTemplate));
	edgestructTemp=(edgeTemplate*)malloc(nEdgeCurr * sizeof(edgeTemplate));
	cellstructTemp=(cellTemplate*)malloc(nCellCurr * sizeof(cellTemplate));
	
	
	if (NULL == vertstructTemp) {
		printf( "malloc failed for vertstructTemp\n");
	
	}
	if (NULL == edgestructTemp) {
		printf( "malloc failed for edgestructTemp\n");
	
	}
	if (NULL == cellstructTemp) {
		printf( "malloc failed for cellstructTemp\n");
	
	}
	for(ii=0;ii<nCellCurr;ii++){
		cellstructTemp[ii].refineVec=(int*)calloc(nLevels,sizeof(int));
	} 
}

// EXECUTION FUNCTIONS
// Index functions

int edgsub(int I, int J, int l, int domSize[dim]){
	/*
	// l is either 0 for dim 1 or 1 for dim 2
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

int vertsub(int I, int J, int domSize[dim]){
	/*
	// domSize is the actual requested domain size rather than the edge domain size
	// I and J vary normallly from 1:domSize+1
	*/
	
	int index;
	index=-1+(I+(J-1)*(domSize[0]+1));
	//printf("domsize %i %i , I %i J %i index %i",domSize[0],domSize[1],I,J,index);
	return(index);
}

int cellsub(int I, int J, int domSize[dim]){
	/*
	// domSize is the actual requested domain size rather than the edge domain size
	// I and J vary normallly from 1:domSize
	*/
	
	int index;
	index=-1+(I+(J-1)*(domSize[0]));
	return(index);
}

// Base Template grid Generation operations
void EdgeIJtoGrid(int domSize[dim],int IJK[dim], int l){
	// l is either 0 or 1 and is used as a mathematical switch
	// Allows the selection of the right position to assign the data
	
	int edgSub;
	
	edgSub=edgsub(IJK[0], IJK[1], l,domSize);
	
	edgestructTemp[edgSub].index=
			(1-l)*(IJK[0]+(IJK[1]-1)*domSize[0]) // case where l=0
			+(l)*((IJK[0]-1)*(domSize[1])+(IJK[1]) // case where l=1
				+((domSize[1]+1)*domSize[0]));
				
	edgestructTemp[edgSub].vertex[0]=(IJK[0]+(IJK[1]-1)*(domSize[0]+1))*(1-l)+ // l=1
			(IJK[0]+(IJK[1]-1)*(domSize[0]+1))*l;// l=2
	edgestructTemp[edgSub].vertex[1]=(1+IJK[0]+(IJK[1]-1)*(domSize[0]+1))*(1-l)+ // l=1
		(IJK[0]+(IJK[1]-1+1)*(domSize[0]+1))*l;// l=2
	
	edgestructTemp[edgSub].cellind[0]=(IJK[0]+(IJK[1]-2)*domSize[0])*(1-l)+// l=1
		(-1+IJK[0]+(IJK[1]-1)*domSize[0])*l;// l=2
	edgestructTemp[edgSub].cellind[1]=(IJK[0]+(IJK[1]-1)*domSize[0])*(1-l)+//ll=1
		(IJK[0]+(IJK[1]-1)*domSize[0])*l; // l=2
	
}

void EdgeIJtoGridLowBound(int domSize[dim],int IJK[dim], int l){
	// l is either 0 or 1 and is used as a mathematical switch
	// Allows the selection of the right position to assign the data
	
	int edgSub;
	
	edgSub=edgsub(IJK[0], IJK[1], l,domSize);
	
	EdgeIJtoGrid(domSize,IJK,l);
	
	edgestructTemp[edgSub].cellind[0]=(-1)*(1-l)+// l=1
		(-2)*l;// l=2
		
	
}

void EdgeIJtoGridHighBound(int domSize[dim],int IJK[dim], int l){
	// l is either 0 or 1 and is used as a mathematical switch
	// Allows the selection of the right position to assign the data
	
	int edgSub;
	
	edgSub=edgsub(IJK[0], IJK[1], l,domSize);
	
	EdgeIJtoGrid(domSize,IJK,l);
	
	edgestructTemp[edgSub].cellind[1]=(-3)*(1-l)+//ll=1
		(-4)*l; // l=2
	
}

void VertexIJGrid(int domSize[dim],int IJK[dim]){
	// For vertices IJK goes to dim+1
	// Assigns vertex data to template structure
	int vertSub;
	
	vertSub=vertsub(IJK[0],IJK[1],domSize);
	vertstructTemp[vertSub].index=IJK[0]+(IJK[1]-1)*(domSize[0]+1);
	
	vertstructTemp[vertSub].coord[0]=((double)(IJK[0]-1))
				/((double)(domSize[0]));
	vertstructTemp[vertSub].coord[1]=((double)(IJK[1]-1))
				/((double)(domSize[1]));
	
}

void CellIJGrid(int domSize[dim],int IJK[dim], int baseRefineLvl){
	// For vertices IJK goes to dim+1
	// baseRefineLvl starts from 1
	int cellSub;
	
	cellSub=cellsub(IJK[0],IJK[1],domSize);
	cellstructTemp[cellSub].index=IJK[0]+(IJK[1]-1)*domSize[0];
	cellstructTemp[cellSub].fill=0;
	cellstructTemp[cellSub].refineLvl=baseRefineLvl;
	cellstructTemp[cellSub].refineVec[baseRefineLvl-1]=IJK[0]+(IJK[1]-1)*domSize[1];

}

// Upper level grid template generation
void AssignEdgestructContent(int domSize[dim]){
	int ii,jj,IJK[dim];
	
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

void AssignVertextructContent(int domSize[dim]){
	
	int ii, jj, IJK[dim];
	
	// Assign all data
	for (ii=1;ii<=domSize[0]+1;ii++){
		for (jj=1;jj<=domSize[1]+1;jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			VertexIJGrid(domSize,IJK);
		}
	}
	
	
}

void AssignCelltructContent(int domSize[dim], int baseRefineLvl){
	
	int ii, jj, IJK[dim];
	
	// Assign all data
	for (ii=1;ii<=domSize[0];ii++){
		for (jj=1;jj<=domSize[1];jj++){
			
			IJK[0]=ii;
			IJK[1]=jj;
			CellIJGrid(domSize,IJK,baseRefineLvl);
		}
	}
	
	
}

void BuildLvlTemplate(int domSize[dim], int baseRefineLvl, int nLevelsInput){

	nLevels=nLevelsInput;
	AllocateGridStruct(domSize);
	AssignVertextructContent(domSize);
	AssignEdgestructContent(domSize);
	AssignCelltructContent(domSize,baseRefineLvl);
	
}

