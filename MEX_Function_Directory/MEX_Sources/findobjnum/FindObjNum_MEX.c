#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "mex.h" /* Always include this */
#include "matrix.h"
/* Constant declaration */
/*int dim() = 2; */
#define MEX_COMPILE

void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList, int *nFound);
void FindObjNum_Deprecated(int *lookup, int *listInd, int *dest, int nLook, int nList);

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){
	
	int *lookup, *listInd,*dest;
	int nList, nLook,ii,nFound;
	double *outPtr;
	double *inputArrayPtr;
	/* Read Integers */
	nLook=(int)mxGetScalar(prhs[2]);
	nList=(int)mxGetScalar(prhs[3]);
	
	/*if(nList>0){*/
		/* Allocation */
		lookup=(int*)calloc(nLook,sizeof(int));
		dest=(int*)calloc(nLook*(nList+1),sizeof(int));
		listInd=(int*)calloc(nList+1,sizeof(int));
	
		/* Assignement */
		inputArrayPtr=mxGetPr(prhs[0]);
		for (ii=0;ii<nLook;ii++){
			lookup[ii]=(int)(*(inputArrayPtr+(ii)));
			/*printf("%i ", cellrefineInd[ii]); */
		}
		inputArrayPtr=mxGetPr(prhs[1]);
		for (ii=0;ii<nList;ii++){
			listInd[ii]=(int)(*(inputArrayPtr+(ii)));
			/*printf("%i ", cellrefineInd[ii]); */
		}
		/* ACTION */
		FindObjNum(lookup, listInd, dest, nLook, nList,&nFound);
		/*dest=(int*)realloc(dest,nFound*sizeof(int));*/
	
		/* Output */
		plhs[0]=(mxCreateDoubleMatrix(1,nFound,mxREAL));
		outPtr=mxGetPr(plhs[0]);
		for (ii=0;ii<nFound;ii++){outPtr[ii]=(double)(dest[ii]+1);}
	
		free(lookup);
		free(listInd);
		free(dest);
	/*} else {plhs[0]=(mxCreateDoubleMatrix(1,0,mxREAL));}*/
	
	return;
}


void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList, int *nFound){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
	int ii,jj,kk,kkStart;
	kk=0;
	for(ii=0;ii<nLook;ii++){
		dest[kk]=-1;
		jj=0;
		kkStart=kk;
		do{
			
			dest[kk]=(-(listInd[jj]!=lookup[ii])*(jj+1))+jj;
			/*dest[ii]=(-(memcmp((listInd+jj),(lookup+ii),sizeof(int))!=0)*(jj+1))+jj; */
			kk=kk+(listInd[jj]==lookup[ii]);
			jj++;
		} while((jj<nList));
		kk=kk+(kkStart==kk);
	}
	*nFound=kk;
}

void FindObjNum_Deprecated(int *lookup, int *listInd, int *dest, int nLook, int nList){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
	int ii,jj;
	
	for(ii=0;ii<nLook;ii++){
		dest[ii]=-1;
		jj=0;
		do{
			
			dest[ii]=(-(listInd[jj]!=lookup[ii])*(jj+1))+jj;
			/*dest[ii]=(-(memcmp((listInd+jj),(lookup+ii),sizeof(int))!=0)*(jj+1))+jj; */
			jj++;
		} while((dest[ii]!=(jj-1)) & (jj<nList));
	}

}
