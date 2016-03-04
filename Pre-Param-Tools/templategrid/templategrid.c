#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "templategrid.h"


// Constant declaration
int dim = 2;

// Global Variables declaration
int nLevels, nCells,nEdges,nVerts;

// Array for active template
cellTemplate *cellstructTemp=NULL;
vertexTemplate *vertstructTemp=NULL;
edgeTemplate *edgestructTemp=NULL;


// Allocation functions

void AllocateGridStruct(int domSize[dim],
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut){
	
	int nCellCurr,nEdgeCurr,nVertCurr, ii;
	nCellCurr=domSize[0]*domSize[1];
	nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	nVertCurr=(domSize[0]+1)*(domSize[1]+1);
	
	//printf("size vertstructTempOut %i\n",nVertCurr * sizeof(vertexTemplate));
	//printf("size edgestructTempOut %i\n",nEdgeCurr * sizeof(edgeTemplate));
	//printf("size cellstructTempOut %i\n",nCellCurr * (sizeof(cellTemplate)+nLevels*sizeof(int)));
	
	
	//printf("Size of Allocated Cell %i\n",nCellCurr * sizeof(cellTemplate));
	
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
	//printf("Accessed Index\n");
	vertstructTemp[vertSub].coord[0]=((double)(IJK[0]-1))
				/((double)(domSize[0]));
	vertstructTemp[vertSub].coord[1]=((double)(IJK[1]-1))
				/((double)(domSize[1]));
	//printf("Accessed coord\n");
	
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

void CalculateNumElements(int domSize[dim],int *nCellCurr,int *nEdgeCurr,int *nVertCurr){
	
	*nCellCurr=domSize[0]*domSize[1];
	*nEdgeCurr=domSize[0]*(1+domSize[1])+domSize[1]*(1+domSize[0]);
	*nVertCurr=(domSize[0]+1)*(domSize[1]+1);
}

void BuildLvlTemplate(int domSize[dim], int baseRefineLvl, int nLevelsInput,
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut){
	
	//cellstructTemp=cellstructTempOut;
	//vertstructTemp=vertstructTempOut;
	//edgestructTemp=edgestructTempOut;
	
	nLevels=nLevelsInput;
	AllocateGridStruct(domSize,cellstructTempOut,vertstructTempOut,edgestructTempOut);
	AllocateGridStruct(domSize,&cellstructTemp,&vertstructTemp,&edgestructTemp);
	printf("Allocated the things\n");
	AssignVertextructContent(domSize);
	AssignEdgestructContent(domSize);
	AssignCelltructContent(domSize,baseRefineLvl);
	printf("Generated the template\n");
	//*cellstructTempOut=cellstructTemp;
	//*vertstructTempOut=vertstructTemp;
	//*edgestructTempOut=edgestructTemp;
	/*
	printf("\n*****  &(cellstructTemp[0].index) in function is: %d\n",&(cellstructTemp[0].index));
	printf("*****  (cellstructTemp[0].index) in function is: %d\n",(cellstructTemp[0].index));
	printf("*****  ADRESS OF &cellStructTemp in function is: %d\n",&cellstructTemp);
	printf("*****  ADRESS OF cellStructTemp in function is: %d\n",cellstructTemp);
	printf("*****  ADRESS OF *cellStructTemp in function is: %d\n",*cellstructTemp);
	*/
	
	int *nCellP, *nEdgeP, *nVertP;
	int nCellCurr, nEdgeCurr, nVertCurr;
	
	nCellP=&nCellCurr;
	nEdgeP=&nEdgeCurr;
	nVertP=&nVertCurr;
	
	CalculateNumElements(domSize,nCellP,nEdgeP,nVertP);
	
	//printf("Size of Copied Cell in function %i\n",sizeof(*cellstructTemp)*(nCellCurr));
	memcpy(*cellstructTempOut,cellstructTemp,sizeof(*cellstructTemp)*(nCellCurr));
	memcpy(*vertstructTempOut,vertstructTemp,sizeof(*vertstructTemp)*(nVertCurr));
	memcpy(*edgestructTempOut,edgestructTemp,sizeof(*edgestructTemp)*(nEdgeCurr));
	/*
	//printf("\n*****  &(cellstructTemp[0].index) in function is: %d\n",&(cellstructTempOut[0].index));
	//printf("*****  (cellstructTemp[0].index) in function is: %d\n",(cellstructTempOut[0].index));
	printf("\n*****  ADRESS OF &cellstructTempOut in function is: %d\n",&cellstructTempOut);
	printf("*****  ADRESS OF cellstructTempOut in function is: %d\n",cellstructTempOut);
	printf("*****  ADRESS OF *cellstructTempOut in function is: %d\n",*cellstructTempOut);
	//*vertstructTempOut=vertstructTemp;
	//*edgestructTempOut=edgestructTemp;
	printf("Returned the template\n");
	printf("cellstructTemp is %i Bytes\n",sizeof(*cellstructTemp));
	printf("cellTemplate is %i Bytes\n",sizeof(cellTemplate));
	printf("cellstructTempOut is %i Bytes\n",sizeof(**cellstructTempOut));
*/
	
	//for (ii=0;ii<sizeof(cellstructTemp);ii++){
	//	((char*)cellstructTemp)[ii]=1;}
		
	free(cellstructTemp);
	free(vertstructTemp);
	free(edgestructTemp);
}

