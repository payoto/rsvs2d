#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <stdbool.h>

#include "gridgen.h"
#include "gridgen.c"
#include "mex.h" /* Always include this */
#include "matrix.h"
// Constant declaration
//int dim() = 2;

#ifndef GRIDGEN_VAR_INCLUDED
#define GRIDGEN_VAR_INCLUDED
int plotFlag=0;
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
#endif
// Function Declarations

// Main text body
void OutputEdgestruct(mxArray **outLHS);
void OutputCellstruct(mxArray **outLHS);
void OutputVertstruct(mxArray **outLHS);

void mexFunction(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[]){
	
	//*cellTemplate cellOut;
	
	//#define Edge_out	plhs[0]
	//#define Vert_out	plhs[1]
	//#define Cell_out	plhs[2]
	printf("\nCalculations Started . . .");
	InitialiseGridFromFile();
	printf("done !\n");
	
	printf("\nMex Output . . .");
	OutputEdgestruct(&plhs[0]);
	OutputVertstruct(&plhs[1]);
	OutputCellstruct(&plhs[2]);
	printf("done !\n");
	
	return;
}

void OutputEdgestruct(mxArray **outLHS){
	
	mxArray **edgeOut1;
	mxArray **edgeOut2;
	double *outArray1;
	double *outArray2;
	const char **fields;
	int nFields=3;
	int ii,jj;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(6*sizeof(char));
	sprintf(fields[0],"index");
	fields[1]=(char*)malloc(10*sizeof(char));
	sprintf(fields[1],"cellindex");
	fields[2]=(char*)malloc(12*sizeof(char));
	sprintf(fields[2],"vertexindex");
	//fields[3]=(char*)malloc(12*sizeof(char));
	//sprintf(fields[3],"orientation");
	
	edgeOut1=(mxArray**)malloc(nEdgeGrid*sizeof(mxArray*));
	edgeOut2=(mxArray**)malloc(nEdgeGrid*sizeof(mxArray*));
	(*outLHS) = mxCreateStructMatrix(nEdgeGrid, 1, nFields, (const char **)fields);
	//edgeOut=mxGetPr(plhs[0]);
	//memcpy(&(edgeOut[0].index),&(edgestruct[0].index),nEdgeGrid*sizeof(edgeTemplate));
	//edgeOut=(mxCreateDoubleMatrix(2,1,mxREAL));
	//outArray=mxGetPr(edgeOut);
	for(ii=0;ii<nEdgeGrid;ii++){
	
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,2,mxREAL));
		edgeOut2[ii]=(mxCreateDoubleMatrix(1,2,mxREAL));
		
		mxSetField((*outLHS),ii,fields[0],mxCreateDoubleScalar(edgestruct[ii].index));
		
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<2;jj++){
			outArray1[jj]=(double)edgestruct[ii].cellind[jj];
		}
		mxSetField((*outLHS),ii,fields[1],edgeOut1[ii]);
		
		outArray2=mxGetPr(edgeOut2[ii]);
		for(jj=0;jj<2;jj++){
			outArray2[jj]=(double)edgestruct[ii].vertex[jj];
		}
		mxSetField((*outLHS),ii,fields[2],edgeOut2[ii]);
		//mxSetField((*outLHS),ii,fields[3],mxCreateDoubleScalar(edgestruct[ii].orientation));
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	free(edgeOut1);
	free(edgeOut2);
}

void OutputCellstruct(mxArray **outLHS){
	
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
	(*outLHS) = mxCreateStructMatrix(nCellGrid, 1, nFields, (const char **)fields);
	//edgeOut=mxGetPr(plhs[0]);
	//memcpy(&(edgeOut[0].index),&(edgestruct[0].index),nEdgeGrid*sizeof(edgeTemplate));
	//edgeOut=(mxCreateDoubleMatrix(2,1,mxREAL));
	//outArray=mxGetPr(edgeOut);
	for(ii=0;ii<nCellGrid;ii++){
	
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,nLevels,mxREAL));
		
		mxSetField((*outLHS),ii,fields[0],mxCreateDoubleScalar(cellstruct[ii].index));
		mxSetField((*outLHS),ii,fields[1],mxCreateDoubleScalar(cellstruct[ii].fill));
		mxSetField((*outLHS),ii,fields[2],mxCreateDoubleScalar(cellstruct[ii].refineLvl));
		
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<nLevels;jj++){
			outArray1[jj]=(double)cellstruct[ii].refineVec[jj];
		}
		mxSetField((*outLHS),ii,fields[3],edgeOut1[ii]);
		
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	free(edgeOut1);
}

void OutputVertstruct(mxArray **outLHS){
	
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
	(*outLHS) = mxCreateStructMatrix(nVertGrid, 1, nFields, (const char **)fields);

	
	for(ii=0;ii<nVertGrid;ii++){
	
		mxSetField((*outLHS),ii,fields[0],mxCreateDoubleScalar(vertstruct[ii].index));
		
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,dim(),mxREAL));
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<dim();jj++){
			outArray1[jj]=(double)vertstruct[ii].coord[jj];
		}
		mxSetField((*outLHS),ii,fields[1],edgeOut1[ii]);
	
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	//free(edgeOut1);
}
