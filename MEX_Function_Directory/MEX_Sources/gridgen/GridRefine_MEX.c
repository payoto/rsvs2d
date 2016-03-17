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
int **newEdgesInd=NULL; int *nNewEdges=NULL, *splitEdgesInd=NULL; int nSplitEdges=0;
int **newCellsInd=NULL, *nNewCells=NULL, *splitCellsInd=NULL, nSplitCells=0;
#endif

// Function Declarations
mxArray *OutputEdgestruct(mxArray *outLHS);
mxArray *OutputCellstruct(mxArray *outLHS);
mxArray *OutputVertstruct(mxArray *outLHS);
void OutputGridStruct(mxArray **outLHS);
void InputGridstruct(mxArray *inLHS);
void InitialiseGridFromFile_MEX(mxArray **plhs);
void RefineGrid_MEX(int nCellRefine, int *indCellRefine,int *posCellRefine);
void OutputConnectStruct(mxArray **outLHS);
mxArray *OutputEdgeConnect();
mxArray *OutputCellConnect();
// Main text body


void mexFunction(int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[]){
	
	int *cellrefineInd=NULL, *cellrefinePos=NULL;
	int nRefine,ii,jj;
	double *inputArrayPtr;
	/*Read in Scalars*/
	nLevels=(int)mxGetScalar(prhs[0]);
	nCellGrid=(int)mxGetScalar(prhs[2]);
	nEdgeGrid=(int)mxGetScalar(prhs[3]);
	nVertGrid=(int)mxGetScalar(prhs[4]);
	nRefine=(int)mxGetScalar(prhs[6]);
	
	/*Allocate Input arrays*/
	AllocateSolutionStruct(nCellGrid, nEdgeGrid, nVertGrid, nLevels,&cellstruct,&vertstruct,&edgestruct);
	levelSize=(int*)calloc(nLevels*dim(),sizeof(int));
	cellrefineInd=(int*)calloc(nRefine,sizeof(int));
	cellrefinePos=(int*)calloc(nRefine,sizeof(int));
	
	/*Add data to input arrays*/
	inputArrayPtr=mxGetPr(prhs[1]);
	for (ii=0;ii<nLevels;ii++){
		for (jj=0;jj<dim();jj++){
			levelSize[dim()*ii+jj]=(int)(*(inputArrayPtr+(dim()*ii+jj)));
			//printf("%i ", levelSize[dim()*ii+jj]);
		}
	}
	inputArrayPtr=mxGetPr(prhs[7]);
	for (ii=0;ii<nRefine;ii++){
		cellrefineInd[ii]=(int)(*(inputArrayPtr+(ii)));
		//printf("%i ", cellrefineInd[ii]);
	}
	inputArrayPtr=mxGetPr(prhs[8]);
	for (ii=0;ii<nRefine;ii++){
		cellrefinePos[ii]=(int)(*(inputArrayPtr+(ii)));
		//printf("%i ", cellrefineInd[ii]);
	}
	InputGridstruct(prhs[5]);
	
	/*Do stuff*/
	printf("\nCalculations Started . . .");
	RefineGrid_MEX(nRefine,cellrefineInd,cellrefinePos);
	printf("done !\n");
	OutputGridStruct(&plhs[0]);
	OutputConnectStruct(&plhs[1]);
	//free(levelSize);
	 ClearWorkSpace();
	free(cellrefineInd);
	free(cellrefinePos);
	/*
	printf("\nMex Output . . .");
	OutputGridStruct(&plhs[0]);
	printf("done !\n");
	*/

	//return;
}

void InputGridstruct(mxArray *inLHS){
	
	mxArray *edgeSubStruct;
	mxArray *cellSubStruct;
	mxArray *vertSubStruct;
	
	double *readPtr;
	int ii,jj;
	
	
	edgeSubStruct=mxGetField(inLHS,0,"edge");
	cellSubStruct=mxGetField(inLHS,0,"cell");
	vertSubStruct=mxGetField(inLHS,0,"vertex");
	
	
	for(ii=0;ii<nEdgeGrid;ii++){
		readPtr=mxGetPr(mxGetField(edgeSubStruct,ii,"index"));
		edgestruct[ii].index=(int)*readPtr;
		
		readPtr=mxGetPr(mxGetField(edgeSubStruct,ii,"cellindex"));
		for(jj=0;jj<2;jj++){ edgestruct[ii].cellind[jj]=(int)*(readPtr+jj);}
		
		readPtr=mxGetPr(mxGetField(edgeSubStruct,ii,"vertexindex"));
		for(jj=0;jj<2;jj++){edgestruct[ii].vertex[jj]=(int)*(readPtr+jj);}
		
		readPtr=mxGetPr(mxGetField(edgeSubStruct,ii,"orientation"));
		edgestruct[ii].orientation=(int)*readPtr;
		
	}
	
	for(ii=0;ii<nCellGrid;ii++){
		readPtr=mxGetPr(mxGetField(cellSubStruct,ii,"index"));
		cellstruct[ii].index=(int)(*readPtr);
		
		readPtr=mxGetPr(mxGetField(cellSubStruct,ii,"fill"));
		cellstruct[ii].fill=(*readPtr);
		
		readPtr=mxGetPr(mxGetField(cellSubStruct,ii,"refineLvl"));
		cellstruct[ii].refineLvl=(int)(*readPtr);
		
		readPtr=mxGetPr(mxGetField(cellSubStruct,ii,"refineVec"));
		for(jj=0;jj<nLevels-1;jj++){ cellstruct[ii].refineVec[jj]=(int)*(readPtr+jj);}
	}
	
	for(ii=0;ii<nVertGrid;ii++){
		readPtr=mxGetPr(mxGetField(vertSubStruct,ii,"index"));
		vertstruct[ii].index=(int)*readPtr;
		
		readPtr=mxGetPr(mxGetField(vertSubStruct,ii,"coord"));
		for(jj=0;jj<dim();jj++){ vertstruct[ii].coord[jj]=*(readPtr+jj);}
		
	}
	OutputGrid(5);
	
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
	printf("Output Grid . . .");
	plhs2[0]=OutputEdgestruct(plhs2[0]);
	printf("Edge . . .");
	plhs2[1]=OutputVertstruct(plhs2[1]);
	printf("Vertex . . .");
	plhs2[2]=OutputCellstruct(plhs2[2]);
	printf("Cell . . .");
	
	mxSetField((*outLHS),0,fields[0],plhs2[0]);
	mxSetField((*outLHS),0,fields[1],plhs2[1]);
	mxSetField((*outLHS),0,fields[2],plhs2[2]);
	printf("done !\n");
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	//free(plhs2);
}

void OutputConnectStruct(mxArray **outLHS){

	const char **fields;
	mxArray *plhs2[3];
	int ii,jj;
	int nFields=2;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(6*sizeof(char));
	sprintf(fields[0],"edge");
	fields[1]=(char*)malloc(5*sizeof(char));
	sprintf(fields[1],"cell");
	
	//plhs2=(mxArray**)malloc(3*sizeof(mxArray*));
	
	(*outLHS) = mxCreateStructMatrix(1, 1, nFields, (const char **)fields);

	printf("Connect . . .");
	//plhs2[0]=OutputEdgestruct(plhs2[0]);
	plhs2[0]=OutputEdgeConnect();
	printf("Edge . . .");
	//plhs2[1]=OutputCellConnect(plhs2[1]);
	plhs2[1]=OutputCellConnect();
	printf("Cell . . .");
	
	mxSetField((*outLHS),0,fields[0],plhs2[0]);
	mxSetField((*outLHS),0,fields[1],plhs2[1]);
	
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

void RefineGrid_MEX(int nCellRefine, int *indCellRefine,int *posCellRefine){

	int ii;//jj,kk;
	int *posEdgeRefine=NULL,*indEdgeRefine=NULL;
	int domSize[dim()];
	int nEdgeRefine=0;
	
	ii=nLevels;
	printf("\n\n");
	printf("\n***** START GRID REFINEMENT %i of %i *****\n",ii-1,nLevels-1);
	GenerateTemplateGrid(ii);
	
	printf("\n    ACTION: Template Generated");
	IdentifyRefineEdge(posCellRefine, indCellRefine,nCellRefine,
			&posEdgeRefine,&indEdgeRefine,&nEdgeRefine,edgestruct,nEdgeGrid);
	
	printf("\n    ACTION: Refinement Targets Identified");
	domSize[0]=levelSize[2*(ii-1)];
	domSize[1]=levelSize[2*(ii-1)+1];
	
	RefineSelectedEdges(domSize,posEdgeRefine,indEdgeRefine,nEdgeRefine);
	printf("\n    ACTION: Edges Refined");

	RefineSelectedCells(domSize,posCellRefine,indCellRefine,nCellRefine);
	printf("\n    ACTION: Cells Refined");

	free(posEdgeRefine);
	free(indEdgeRefine);
	DeAllocateTemplate(domSize, cellCurrentTemplate, edgeCurrentTemplate, vertCurrentTemplate);
	
	OutputGrid(ii);
	printf("\n\n");
	
	
	
}

mxArray *OutputEdgeConnect(){
	mxArray *outLHS;
	mxArray **edgeOut1;
	double *outArray1;
	const char **fields;
	int nFields=2;
	int ii,jj;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(4*sizeof(char));
	sprintf(fields[0],"old");
	fields[1]=(char*)malloc(4*sizeof(char));
	sprintf(fields[1],"new");

	
	edgeOut1=(mxArray**)malloc(nSplitEdges*sizeof(mxArray*));
	//printf("\n\n nSplitEdges=%i \n",nSplitEdges);
	(outLHS) = mxCreateStructMatrix(nSplitEdges, 1, nFields, (const char **)fields);
	//edgeOut=mxGetPr(plhs[0]);
	//memcpy(&(edgeOut[0].index),&(edgestruct[0].index),nEdgeGrid*sizeof(edgeTemplate));
	//edgeOut=(mxCreateDoubleMatrix(2,1,mxREAL));
	//outArray=mxGetPr(edgeOut);
	
	for(ii=0;ii<nSplitEdges;ii++){
	
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,nNewEdges[ii],mxREAL));
		
		mxSetField((outLHS),ii,fields[0],mxCreateDoubleScalar(splitEdgesInd[ii]));
		
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<nNewEdges[ii];jj++){
			outArray1[jj]=(double)(*((*(newEdgesInd+ii))+jj));
		}
		mxSetField((outLHS),ii,fields[1],edgeOut1[ii]);
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	free(edgeOut1);
	
	return(outLHS);
}

mxArray *OutputCellConnect(){
	mxArray *outLHS;
	mxArray **edgeOut1;
	double *outArray1;
	const char **fields;
	int nFields=2;
	int ii,jj;
	
	fields=(char**)malloc(nFields*sizeof(char*));
	fields[0]=(char*)malloc(4*sizeof(char));
	sprintf(fields[0],"old");
	fields[1]=(char*)malloc(4*sizeof(char));
	sprintf(fields[1],"new");

	
	edgeOut1=(mxArray**)malloc(nSplitCells*sizeof(mxArray*));
	//printf("\n\n nSplitCells=%i \n",nSplitCells);
	(outLHS) = mxCreateStructMatrix(nSplitCells, 1, nFields, (const char **)fields);
	//edgeOut=mxGetPr(plhs[0]);
	//memcpy(&(edgeOut[0].index),&(edgestruct[0].index),nEdgeGrid*sizeof(edgeTemplate));
	//edgeOut=(mxCreateDoubleMatrix(2,1,mxREAL));
	//outArray=mxGetPr(edgeOut);
	
	for(ii=0;ii<nSplitCells;ii++){
	
		edgeOut1[ii]=(mxCreateDoubleMatrix(1,nNewCells[ii],mxREAL));
		
		mxSetField((outLHS),ii,fields[0],mxCreateDoubleScalar(splitCellsInd[ii]));
		
		outArray1=mxGetPr(edgeOut1[ii]);
		for(jj=0;jj<nNewCells[ii];jj++){
			outArray1[jj]=(double)(*((*(newCellsInd+ii))+jj));
		}
		mxSetField((outLHS),ii,fields[1],edgeOut1[ii]);
	}
	
	for (ii=0;ii<nFields;ii++){
		free(fields[ii]);
	}
	free(fields);
	free(edgeOut1);
	
	return(outLHS);
}
