#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "meshsym.h"

#ifndef MESHSYM_VAR_INCLUDED
#define MESHSYM_VAR_INCLUDED
int nCell, nEdge, nVert,nCellKeep, nEdgeKeep, nVertKeep,nBound;
int cutDir,cutSide;
double cutLoc,cutEps;
int *edgeArray=NULL,*vertInd=NULL,*keepVertPos=NULL,*keepCellPos=NULL,*keepEdgePos=NULL;
int *cellMatchList=NULL,*vertMatchList=NULL,*keepCellLog=NULL,*vertDelLog=NULL;
double *cellVol=NULL,*vertCoord=NULL;

#endif



int main(){
	
	int nVertNeg,nEdgeDel;
	int *posVertNeg=NULL,*posEdgeDel=NULL;
	
	 DataIn();
	posVertNeg=FindNegVertex(&nVertNeg);
	posEdgeDel=FindDelEdge(&nEdgeDel,posVertNeg, nVertNeg);
	printf("%i negative vertices of %i \n",nVertNeg,nVert);
	printf("%i negative Edges of %i \n",nEdgeDel,nEdge);
	printf("%i Cells to keep if %i \n",nCellKeep,nCell);
	BuildCellMatchList();
	BuildVertMatchList();
	WriteOutGrid();
	WriteOutSymPlane();
	free(posVertNeg);
	free(posEdgeDel);
	FreeAllMeshSym();
	return(0);
}

void FreeAllMeshSym(){
free(edgeArray);
free(vertInd);
free(keepVertPos);
free(keepCellPos);
free(keepEdgePos);
free(cellMatchList);
free(vertMatchList);
free(keepCellLog);
free(vertDelLog);
free(cellVol);
free(vertCoord);

}

void DataIn(){
	/* Reads in data from external files */
	
	FILE *cellgridFID;
	int ii,jj;
	
	cellgridFID=fopen("griduns","r");
	
	if((cellgridFID!=NULL)){
	
		fscanf(cellgridFID,"%i",&nCell);
		fscanf(cellgridFID,"%i",&nEdge);
		fscanf(cellgridFID,"%i",&nVert);
		
		edgeArray=(int*)calloc(nEdge*4, sizeof(int));
		for(ii=0;ii<nEdge;ii++){
			for(jj=0;jj<4;jj++){
				fscanf(cellgridFID,"%i",&(edgeArray[ii*4+jj]));
			}
		}
		
		vertInd=(int*)calloc(nVert, sizeof(int));
		vertCoord=(double*)calloc(nVert*2, sizeof(double));
		for(ii=0;ii<nVert;ii++){
			fscanf(cellgridFID,"%i",&(vertInd[ii]));
			for(jj=0;jj<2;jj++){
				fscanf(cellgridFID,"%lf",&(vertCoord[ii*2+jj]));
			}
		}
		
		cellVol=(double*)calloc(nCell, sizeof(double));
		for(ii=0;ii<nCell;ii++){
			fscanf(cellgridFID,"%lf",&(cellVol[ii]));
		}
		
		
		/*printf("Cell Grid Structure succesfully read in\n\n"); */
		} else {
		perror("Data file failed to open!");
		exit(-1);
	}
	
	fclose(cellgridFID);
	cellgridFID=fopen("symsettings","r");
	
	if((cellgridFID!=NULL)){
		fscanf(cellgridFID,"%i",&cutDir);
		fscanf(cellgridFID,"%i",&cutSide);
		fscanf(cellgridFID,"%lf",&cutLoc);
		fscanf(cellgridFID,"%lf",&cutEps);
	} else {
		perror("Settings file failed to open!");
		exit(-1);
	}
	fclose(cellgridFID);
}

int *FindNegVertex(int *nVertNeg){
	int ii,jj,jjKeep,vertIndMax,cond;
	int *vertNegPos;
	
	
	vertIndMax=max_array(vertInd, nVert);
	vertDelLog=(int*)calloc(vertIndMax+1,sizeof(int));
	vertNegPos=(int*)calloc(nVert+1,sizeof(int));
	keepVertPos=(int*)calloc(nVert+1,sizeof(int));
	jj=0;
	jjKeep=0;
	for(ii=0;ii<nVert;ii++){
		
		vertNegPos[jj]=ii;
		keepVertPos[jjKeep]=ii;
		cond=((cutSide*(vertCoord[ii*2+cutDir]-cutLoc))<(-cutEps)); 
		jj=jj+cond;/*(vertCoord[ii*2+1]<pow(10,-10))*/
		jjKeep=jjKeep+!(cond);
		vertDelLog[vertInd[ii]]=cond;
	}
	
	*nVertNeg=jj;
	vertNegPos=realloc(vertNegPos,(*nVertNeg)*sizeof(int));
	nVertKeep=jjKeep;
	keepVertPos=realloc(keepVertPos,nVertKeep*sizeof(int));
	
	return(vertNegPos);
}

int *FindDelEdge(int *nEdgeDel,int *posVertNeg,int nVertNeg){
	int ii,jj,kk,jjKeep,cond,ind;
	int *edgeDelPos=NULL;
	int dump[2];
	
	
	edgeDelPos=(int*)calloc(nEdge,sizeof(int));
	keepEdgePos=(int*)calloc(nEdge,sizeof(int));
	keepCellLog=(int*)calloc(nCell+1,sizeof(int));
	keepCellPos=(int*)calloc(2*nEdge,sizeof(int));
	jj=0;
	jjKeep=0;
	dump[0]=0;
	dump[1]=0;
	for(ii=0;ii<nCell;ii++) keepCellLog[ii]=1;
	for(ii=0;ii<nEdge;ii++){
		
		edgeDelPos[jj]=ii;
		keepEdgePos[jjKeep]=ii;
		keepCellPos[jjKeep*2]=edgeArray[ii*4+2];
		keepCellPos[jjKeep*2+1]=edgeArray[ii*4+3];
		/*
		cond=0;
		kk=0;
		while ((cond==0) & (kk<nVertNeg)){
			cond=(edgeArray[ii*4]==vertInd[posVertNeg[kk]]) | (edgeArray[ii*4+1]==vertInd[posVertNeg[kk]]); 
			
			kk++;
		}
		*/
		cond=vertDelLog[edgeArray[ii*4]] | vertDelLog[edgeArray[ii*4+1]]; 
		
		ind=(edgeArray[ii*4+2]-1)*(edgeArray[ii*4+2]>0)+nCell*(!(edgeArray[ii*4+2]>0));
		//  def cur ind = cellPos(fortran) convert to C, math switch for positive index + case when negative index. 
		keepCellLog[ind]=keepCellLog[ind]*(!cond);
		ind=(edgeArray[ii*4+3]-1)*(edgeArray[ii*4+3]>0)+nCell*(!(edgeArray[ii*4+3]>0));
		keepCellLog[ind]=keepCellLog[ind]*(!cond);
		jj=jj+cond; /*(vertCoord[ii*2+1]<pow(10,-10))*/
		jjKeep=jjKeep+!cond;
	}
	
	
	
	*nEdgeDel=jj;
	edgeDelPos=realloc(edgeDelPos,(*nEdgeDel)*sizeof(int));
	nEdgeKeep=jjKeep;
	keepEdgePos=realloc(keepEdgePos,(nEdgeKeep)*sizeof(int));
	
	jj=0;
	for (ii=0;ii<3;ii++){keepCellPos[ii]=ii-3;}
	for (ii=0;ii<nCell;ii++){
		keepCellPos[jj+3]=ii+1;
		jj=jj+keepCellLog[ii];
	}
	keepCellPos=realloc(keepCellPos,(jj+3)*sizeof(int));
	RemoveIdenticalEntries_int(&keepCellPos, (jj+3), &nCellKeep);

	return(edgeDelPos);
}


void BuildCellMatchList(){
	int ii,kk=0;
	int cellPosMin,cellPosMax;
	
	cellPosMin=min_array(keepCellPos, nCellKeep);
	cellPosMax=max_array(keepCellPos, nCellKeep);
	
	cellMatchList=(int*)calloc(nCell-cellPosMin+1,sizeof(int));
	for (ii=0;ii<nCell-cellPosMin+1;ii++){cellMatchList[ii]=-3;}
	while (keepCellPos[kk]<=0) kk++;
	nBound=kk;
	printf("%i %i %i %i %i\n",cellPosMin,kk, cellPosMax,nCellKeep, nCell);
	for (ii=0;ii<nCellKeep;ii++){
		cellMatchList[keepCellPos[ii]-cellPosMin]=(keepCellPos[ii])*(ii<kk)+(ii-kk+1)*(ii>=kk);
		}
	/*cellPosMin=min_array(cellMatchList, nCell-cellPosMin);
	cellPosMax=max_array(cellMatchList, nCell-cellPosMin);
	printf("%i %i %i %i %i\n",cellPosMin,kk, cellPosMax,nCellKeep, nCell);*/
	//for (ii=0;ii<nCellKeep+2;ii++){printf("%6i ",keepCellPos[ii]);}
}


void BuildCellMatchList2(){
	int ii,kk=0,jj=0;
	int cellPosMin,cellPosMax;
	
	cellPosMin=min_array(keepCellPos, nCellKeep);
	cellPosMax=max_array(keepCellPos, nCellKeep);
	
	cellMatchList=(int*)calloc(nCell-cellPosMin+1,sizeof(int));
	for (ii=0;ii<nCell-cellPosMin+1;ii++){cellMatchList[ii]=-3;}
	while (keepCellPos[kk]<=0) kk++;
	nBound=kk;
	printf("%i %i %i %i %i\n",cellPosMin,kk, cellPosMax,nCellKeep, nCell);
	for (ii=0;ii<kk;ii++){
		cellMatchList[keepCellPos[ii]-cellPosMin]=(keepCellPos[ii])*(ii<kk);
		}
	for (ii=0;ii<nCell;ii++){
		cellMatchList[1-cellPosMin+ii]=cellMatchList[1-cellPosMin+ii]*(!keepCellLog[ii])+jj*(keepCellLog[ii]);
		jj=jj+(keepCellLog[ii]);
		}
	/*cellPosMin=min_array(cellMatchList, nCell-cellPosMin);
	cellPosMax=max_array(cellMatchList, nCell-cellPosMin);
	printf("%i %i %i %i %i\n",cellPosMin,kk, cellPosMax,nCellKeep, nCell);*/
	//for (ii=0;ii<nCellKeep+2;ii++){printf("%6i ",keepCellPos[ii]);}
}

void BuildVertMatchList(){
	int ii;
	vertMatchList=(int*)calloc(nVert+1,sizeof(int));
	for (ii=0;ii<nVertKeep;ii++){
		vertMatchList[vertInd[keepVertPos[ii]]]=ii+1;
		}
	
}

void WriteOutGrid(){
	int ii,jj,kk;
	FILE *cellgridFID;
	int cellPosMin,cellPosMax;
	
	cellPosMin=min_array(keepCellPos, nCellKeep);

	cellgridFID=fopen("griduns2","w");

	if((cellgridFID!=NULL)){
	
		fprintf(cellgridFID," %11i %11i %11i\n",nCellKeep-nBound,nEdgeKeep,nVertKeep);
		for (ii=0;ii<nEdgeKeep;ii++){
			
				fprintf(cellgridFID," %11i",vertMatchList[edgeArray[keepEdgePos[ii]*4]]);
				fprintf(cellgridFID," %11i",vertMatchList[edgeArray[keepEdgePos[ii]*4+1]]);
				fprintf(cellgridFID," %11i",cellMatchList[-cellPosMin+edgeArray[keepEdgePos[ii]*4+2]]);
				fprintf(cellgridFID," %11i\n",cellMatchList[-cellPosMin+edgeArray[keepEdgePos[ii]*4+3]]);
			
		}
		for (ii=0;ii<nVertKeep;ii++){
			fprintf(cellgridFID," %11i",vertMatchList[vertInd[keepVertPos[ii]]]);
			fprintf(cellgridFID," %30.20E",vertCoord[keepVertPos[ii]*2]);
			fprintf(cellgridFID," %30.20E\n",vertCoord[keepVertPos[ii]*2+1]);
		}
		for (ii=nBound;ii<nCellKeep;ii++){
			fprintf(cellgridFID," %30.20E\n",cellVol[keepCellPos[ii]-1]);
		}
	} else {
		perror("Output file failed to open!");
		exit(-1);
	}
	fclose(cellgridFID);
}


void WriteOutSymPlane(){
	int ii,jj,kk;
	FILE *cellgridFID;
	int cellPosMin,cellPosMax,nVertWritten,currVert;
	int *vertWrittenLog;


	vertWrittenLog=(int*)calloc(nVertKeep,sizeof(int));
	nVertWritten=0;
	
	cellPosMin=min_array(keepCellPos, nCellKeep);

	cellgridFID=fopen("symplane.xyz","w");

	if((cellgridFID!=NULL)){
	
		
		for (ii=0;ii<nEdgeKeep;ii++){
			
			if ((cellMatchList[-cellPosMin+edgeArray[keepEdgePos[ii]*4+2]]==-3)
				| (cellMatchList[-cellPosMin+edgeArray[keepEdgePos[ii]*4+3]]==-3) ){

				for (jj=0;jj<2;jj++){
					currVert=vertMatchList[edgeArray[keepEdgePos[ii]*4+jj]]-1;
					if (vertWrittenLog[currVert]==0){
						fprintf(cellgridFID," %30.20E",vertCoord[keepVertPos[currVert]*2]);
						fprintf(cellgridFID," %30.20E",vertCoord[keepVertPos[currVert]*2+1]);
						fprintf(cellgridFID," %30.20E\n",0);
						vertWrittenLog[currVert]=1;
						nVertWritten++;
					}
				}
			}
			
			
		}
		printf("Symmetry File written: %i vertices\n",nVertWritten);
		
	} else {
		perror("Output file failed to open!");
		exit(-1);
	}
	fclose(cellgridFID);
	free(vertWrittenLog);
}

void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
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

void RemoveIdenticalEntries_int(int **array, int nArray, int *nNewArray){
	
	int jj,ii;
	
	qsort((*array),nArray,sizeof(int),SortOrder);
	jj=1;
	/*printf("\n"); */
	for (ii=1;ii<nArray;ii++){
		if ((*array)[ii]>(*array)[jj-1]){
			
			/*printf("\n%i %i %i",ii,jj,nArray); */
			(*array)[jj]=(*array)[ii];
			jj++;
			/*printf("%i ", jj); */
		}
	}
	if (jj<nArray){
		(*array)=(int*)realloc((*array),jj*sizeof(int));
	}
	*nNewArray=jj;
}

int SortOrder(const void * elem1, const void * elem2) {
    int f = *((int*)elem1);
    int s = *((int*)elem2);
	
    return (f > s) - (f < s);
}

	/* Min Max functions  */
	/* (see Macro) */

double fsum(double *array, int nArray){
	int ii;
	double sum;
	sum=0;
	for(ii=0;ii<nArray;ii++){
		sum=sum+array[ii];
	}
	
	return(sum);
}
double fprod(double *array, int nArray){
	int ii;
	double sum;
	sum=0;
	for(ii=0;ii<nArray;ii++){
		sum=sum*array[ii];
	}
	
	return(sum);
}
int isum(int *array, int nArray){
	int ii;
	int sum;
	sum=0;
	for(ii=0;ii<nArray;ii++){
		sum=sum+array[ii];
	}
	
	return(sum);
}
int iprod(int *array, int nArray){
	int ii;
	int sum;
	sum=0;
	for(ii=0;ii<nArray;ii++){
		sum=sum*array[ii];
	}
	
	return(sum);
}

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


