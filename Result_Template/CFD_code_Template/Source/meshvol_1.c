#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "meshvol.h"

#ifndef MESHSYM_VAR_INCLUDED
#define MESHSYM_VAR_INCLUDED
int nCell, nEdge, nVert;

int *edgeArray=NULL,*vertInd=NULL;
double *cellVol=NULL,*cellVol2=NULL,*vertCoord=NULL,*volEdge=NULL;

#endif



int main(){
	int ii,kk;
	//int nVertNeg,nEdgeDel;
	//int *posVertNeg=NULL,*posEdgeDel=NULL;
	
	DataIn();
	CalculateVolEdge();
	CalculateCellVol();
	WriteOutGridNewVolume();

	// for (ii=0;ii<nCell;ii++){
	// 	kk=(fabs((cellVol[ii]-cellVol2[ii])/cellVol[ii])>pow(10,-8));
	// }
	// printf("The difference in volumes are %i",kk);
	//free(posVertNeg);
	//free(posEdgeDel);
	//FreeAllMeshSym();
	return(0);
}

void FreeAllMeshSym(){
free(edgeArray);
free(vertInd);
free(cellVol);
free(volEdge);
free(cellVol2);
free(vertCoord);

}

void CalculateVolEdge(){

	int ii,currVert[2];
	double x1,y1,x2,y2;

	volEdge=(double*)calloc(nEdge,sizeof(double));
	
	for (ii=0;ii<nEdge;ii++){
		// currVert[1]=edgeArray[ii*4];
		// currVert[2]=edgeArray[ii*4+1];
		FindObjNum(&(edgeArray[ii*4]), vertInd, currVert, 2, nVert);

		x1=vertCoord[currVert[0]*2];
		y1=vertCoord[currVert[0]*2+1];
		x2=vertCoord[currVert[1]*2];
		y2=vertCoord[currVert[1]*2+1];

		// THis line implements green's theorem for each edge
		volEdge[ii]=0.5*((y1-y2)*(x2+x1)/2.0+(x2-x1)*(y2+y1)/2.0);

	}
	printf("Calculation of Volume Edges done!\n");

}

void CalculateCellVol(){

	int ii;

	cellVol2=(double*)calloc(nCell,sizeof(double));
	for (ii=0;ii<nEdge;ii++){
		
		
		// assumes the first cell is the one on the left and the second the one on the right of the edge.
		if (edgeArray[ii*4+2]>0){
			cellVol2[edgeArray[ii*4+2]-1]=cellVol2[edgeArray[ii*4+2]-1]-volEdge[ii];
		}
		if (edgeArray[ii*4+3]>0){
			cellVol2[edgeArray[ii*4+3]-1]=cellVol2[edgeArray[ii*4+3]-1]+volEdge[ii];
		}


	}
	printf("Calculation of Cell Volumes done!\n");
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
	
}


void WriteOutGridNewVolume(){
	int ii,jj,kk;
	FILE *cellgridFID;
	

	cellgridFID=fopen("griduns2","w");

	if((cellgridFID!=NULL)){
	
		fprintf(cellgridFID," %11i %11i %11i\n",nCell,nEdge,nVert);
		for (ii=0;ii<nEdge;ii++){
			
				fprintf(cellgridFID," %11i",edgeArray[ii*4]);
				fprintf(cellgridFID," %11i",edgeArray[ii*4+1]);
				fprintf(cellgridFID," %11i",edgeArray[ii*4+2]);
				fprintf(cellgridFID," %11i\n",edgeArray[ii*4+3]);
			
		}
		for (ii=0;ii<nVert;ii++){
			fprintf(cellgridFID," %11i",vertInd[ii]);
			fprintf(cellgridFID," %30.20E",vertCoord[ii*2]);
			fprintf(cellgridFID," %30.20E\n",vertCoord[ii*2+1]);
		}
		for (ii=0;ii<nCell;ii++){
			fprintf(cellgridFID," %30.20E\n",cellVol2[ii]);
		}
	} else {
		perror("Output file failed to open!");
		exit(-1);
	}
	fclose(cellgridFID);
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


