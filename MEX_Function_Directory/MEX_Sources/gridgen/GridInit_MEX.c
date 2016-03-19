#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <stdbool.h>


#include "mex.h" /* Always include this */
#include "matrix.h"
// Constant declaration
//int dim() = 2;
#define MEX_COMPILE

#include "gridgen.h"
#include "gridgen.c"

#ifndef GRIDGEN_VAR_INCLUDED
#define GRIDGEN_VAR_INCLUDED
int plotFlag=0;
int outCount=0;
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
// Dependancy Arrays
int **newEdgesInd=NULL, *nNewEdges=NULL, *splitEdgesInd=NULL, nSplitEdges=0;
int **newCellsInd=NULL, *nNewCells=NULL, *splitCellsInd=NULL, nSplitCells=0;
#endif
// Function Declarations

// Main text body
mxArray *OutputEdgestruct(mxArray *outLHS);
mxArray *OutputCellstruct(mxArray *outLHS);
mxArray *OutputVertstruct(mxArray *outLHS);
void OutputGridStruct(mxArray **outLHS);
void InitialiseGridFromFile_MEX(mxArray **plhs);

void mexFunction(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[]){
	
	
	printf("\nCalculations Started . . .");
	InitialiseGridFromFile_MEX(plhs);
	printf("done !\n");
	
	ClearWorkSpace();
	return;
}

void InitialiseGridFromFile_MEX(mxArray **plhs){
	 //Arrays
	//printf("%i\n",sizeof(edgeTemplate));
	DataIn();
	Allocatecelldatstruct();
	DataArrayToStruct();
	
	printf("\n***** START GRID INITIALISATION *****\n");
	GridInitialisation();
	printf("\n    ACTION: Starting Grid Generated");
	OutputGrid(1);
	printf("\nMex Output . . .");
	OutputGridStruct(&plhs[0]);
	printf("done !\n");
	RefineGrid();
	printf("\nMex Output . . .");
	OutputGridStruct(&plhs[1]);
	printf("done !\n");
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

void OutputGridStruct(mxArray **outLHS){

	const char **fields;
	mxArray *plhs2[3];
	int ii,jj;
	int nFields=3;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(6*sizeof(char));
	sprintf(fields[0],"edge");
	fields[1]=(char*)malloc(7*sizeof(char));
	sprintf(fields[1],"vertex");
	fields[2]=(char*)malloc(7*sizeof(char));
	sprintf(fields[2],"cell");
	
	//plhs2=(mxArray**)malloc(3*sizeof(mxArray*));
	
	(*outLHS) = mxCreateStructMatrix(1, 1, nFields, (const char **)fields);
	//plhs2[0]=mxGetField((*outLHS),1,fields[0]);
	//plhs2[1]=mxGetField((*outLHS),1,fields[1]);
	//plhs2[2]=mxGetField((*outLHS),1,fields[2]);
	//OutputEdgestruct(&plhs2);
	//printf("is Struct? %s\n",mxGetFieldNameByNumber(*outLHS,1));
	printf("Grid . . .");
	plhs2[0]=OutputEdgestruct(plhs2[0]);
	printf("Edge . . .");
	plhs2[1]=OutputVertstruct(plhs2[1]);
	printf("Vertex . . .");
	plhs2[2]=OutputCellstruct(plhs2[2]);
	printf("Cell . . .");
	
	mxSetField((*outLHS),0,fields[0],plhs2[0]);
	mxSetField((*outLHS),0,fields[1],plhs2[1]);
	mxSetField((*outLHS),0,fields[2],plhs2[2]);
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	//free(plhs2);
}

mxArray *OutputEdgestruct(mxArray *outLHS){
	
	mxArray **edgeOut1;
	mxArray **edgeOut2;
	double *outArray1;
	double *outArray2;
	const char **fields;
	int nFields=4;
	int ii,jj;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(6*sizeof(char));
	sprintf(fields[0],"index");
	fields[1]=(char*)malloc(10*sizeof(char));
	sprintf(fields[1],"cellindex");
	fields[2]=(char*)malloc(12*sizeof(char));
	sprintf(fields[2],"vertexindex");
	fields[3]=(char*)malloc(12*sizeof(char));
	sprintf(fields[3],"orientation");
	
	edgeOut1=(mxArray**)malloc(nEdgeGrid*sizeof(mxArray*));
	edgeOut2=(mxArray**)malloc(nEdgeGrid*sizeof(mxArray*));
	(outLHS) = mxCreateStructMatrix(nEdgeGrid, 1, nFields, (const char **)fields);
	//edgeOut=mxGetPr(plhs[0]);
	//memcpy(&(edgeOut[0].index),&(edgestruct[0].index),nEdgeGrid*sizeof(edgeTemplate));
	//edgeOut=(mxCreateDoubleMatrix(2,1,mxREAL));
	//outArray=mxGetPr(edgeOut);
	for(ii=0;ii<nEdgeGrid;ii++){
	
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,2,mxREAL));
		edgeOut2[ii]=(mxCreateDoubleMatrix(1,2,mxREAL));
		
		mxSetField((outLHS),ii,fields[0],mxCreateDoubleScalar(edgestruct[ii].index));
		
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<2;jj++){
			outArray1[jj]=(double)edgestruct[ii].cellind[jj];
		}
		mxSetField((outLHS),ii,fields[1],edgeOut1[ii]);
		
		outArray2=mxGetPr(edgeOut2[ii]);
		for(jj=0;jj<2;jj++){
			outArray2[jj]=(double)edgestruct[ii].vertex[jj];
		}
		mxSetField((outLHS),ii,fields[2],edgeOut2[ii]);
		mxSetField((outLHS),ii,fields[3],mxCreateDoubleScalar(edgestruct[ii].orientation));
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	free(edgeOut1);
	free(edgeOut2);
	return(outLHS);
}

mxArray *OutputCellstruct(mxArray *outLHS){
	
	mxArray **edgeOut1;
	double *outArray1;
	
	const char **fields;
	int ii,jj;
	int nFields=4;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(6*sizeof(char));
	sprintf(fields[0],"index");
	fields[1]=(char*)malloc(5*sizeof(char));
	sprintf(fields[1],"fill");
	fields[2]=(char*)malloc(10*sizeof(char));
	sprintf(fields[2],"refineLvl");
	fields[3]=(char*)malloc(10*sizeof(char));
	sprintf(fields[3],"refineVec");
	
	edgeOut1=(mxArray**)malloc(nEdgeGrid*sizeof(mxArray*));
	(outLHS) = mxCreateStructMatrix(nCellGrid, 1, nFields, (const char **)fields);
	//edgeOut=mxGetPr(plhs[0]);
	//memcpy(&(edgeOut[0].index),&(edgestruct[0].index),nEdgeGrid*sizeof(edgeTemplate));
	//edgeOut=(mxCreateDoubleMatrix(2,1,mxREAL));
	//outArray=mxGetPr(edgeOut);
	for(ii=0;ii<nCellGrid;ii++){
	
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,nLevels,mxREAL));
		
		mxSetField((outLHS),ii,fields[0],mxCreateDoubleScalar(cellstruct[ii].index));
		mxSetField((outLHS),ii,fields[1],mxCreateDoubleScalar(cellstruct[ii].fill));
		mxSetField((outLHS),ii,fields[2],mxCreateDoubleScalar(cellstruct[ii].refineLvl));
		
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<nLevels;jj++){
			outArray1[jj]=(double)cellstruct[ii].refineVec[jj];
		}
		mxSetField((outLHS),ii,fields[3],edgeOut1[ii]);
		
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	free(edgeOut1);
	return(outLHS);
}

mxArray *OutputVertstruct(mxArray *outLHS){
	
	mxArray **edgeOut1;
	double *outArray1;
	const char **fields;
	int ii,jj;
	int nFields=2;
	
	
	fields=(const char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(6*sizeof(char));
	sprintf(fields[0],"index");
	fields[1]=(char*)malloc(6*sizeof(char));
	sprintf(fields[1],"coord");
	
	edgeOut1=(mxArray**)malloc(nVertGrid*sizeof(mxArray*));
	(outLHS) = mxCreateStructMatrix(nVertGrid, 1, nFields, (const char **)fields);

	
	for(ii=0;ii<nVertGrid;ii++){
	
		mxSetField((outLHS),ii,fields[0],mxCreateDoubleScalar(vertstruct[ii].index));
		
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,dim(),mxREAL));
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<dim();jj++){
			outArray1[jj]=(double)vertstruct[ii].coord[jj];
		}
		mxSetField((outLHS),ii,fields[1],edgeOut1[ii]);
	
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	return(outLHS);
	//free(edgeOut1);
}
