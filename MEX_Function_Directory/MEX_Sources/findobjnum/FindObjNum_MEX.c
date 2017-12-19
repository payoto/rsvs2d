#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "mex.h" /* Always include this */
#include "matrix.h"
/*#include "llist.h"*/
/* Constant declaration */
/*int dim() = 2; */
#define MEX_COMPILE

int FindObjNum(double *lookup, double *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound);
void FindObjNum_Deprecated(double *lookup, double *listInd, int *dest, int nLook, int nList);
int FindObjNumHash(double *lookup, double *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound);
int BuildFindObjHashTable(double *lookup, double *listInd, int nLook, int nList,int *hashKey,int *hashValue, int nHash);

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){
	
	int *dest;
	int nList, nLook,ii,nFound,nMatch,isErr;
	double *outPtr,*lookup, *listInd;
	double *inputArrayPtr;
	/* Read Integers */
	nLook=(int)mxGetScalar(prhs[2]);
	nList=(int)mxGetScalar(prhs[3]);
	
	/*if(nList>0){*/
		/* Allocation */
	nMatch=(nLook+(nList+1))*2;
	dest=(int*)calloc(nMatch,sizeof(int));
	
		/* Assignement */
	lookup=mxGetPr(prhs[0]);
	
	listInd=mxGetPr(prhs[1]);
	
		/* ACTION */

	isErr=(nList+3*nLook)>(nList*nLook/2);
	if (!isErr){
		FindObjNumHash(lookup, listInd, &dest, nLook, nList,&nMatch,&nFound);
	}
	if (isErr){
		FindObjNum(lookup, listInd, &dest, nLook, nList,&nMatch,&nFound);
	}
		/*dest=(int*)realloc(dest,nFound*sizeof(int));*/
	
		/* Output */
	plhs[0]=(mxCreateDoubleMatrix(1,nFound,mxREAL));
	outPtr=mxGetPr(plhs[0]);
	for (ii=0;ii<nFound;ii++){outPtr[ii]=(double)(dest[ii]+1);}
		
	/*free(lookup);*/
	/*free(listInd);*/
	free(dest);
	/*} else {plhs[0]=(mxCreateDoubleMatrix(1,0,mxREAL));}*/
	
	return;
}


int FindObjNum(double *lookup, double *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
	int ii,jj,kk,kkStart;

	if (nList>0){
		kk=0;
		for(ii=0;ii<nLook;ii++){
			(*dest)[kk]=-1;
			jj=0;
			kkStart=kk;
			if ((*nMatch - kk)<nList){
				*nMatch=*nMatch+nList+nLook;
				*dest=realloc(*dest,(*nMatch)*sizeof(int));
				if (*dest==NULL){
					return(-1);
				}
			}
			do{
			
				(*dest)[kk]=(-(*(listInd+jj)!=*(lookup+ii))*(jj+1))+jj;
				/*dest[ii]=(-(memcmp((listInd+jj),(lookup+ii),sizeof(int))!=0)*(jj+1))+jj; */
				kk=kk+(*(listInd+jj)==*(lookup+ii));
				jj++;
			} while((jj<nList));
			kk=kk+(kkStart==kk);
		}
		*nFound=kk;
	} else {
		for (ii=0;ii<nLook;ii++){
			(*dest)[ii]=-1;
		}
		*nFound=ii;
	}
	return(0);
}

int BuildFindObjHashTable(double *lookup, double *listInd, int nLook, int nList,int *hashKey,int *hashValue,int nHash){

	int ii,jjErr,ind,kk;

	jjErr=0;
	for(ii=0;ii<nList;ii++){
		ind=(int)listInd[ii]%nHash;
		kk=0;
		while(hashValue[ind]!=-1  && kk<=nHash){
			if (hashKey[ind]==(int)listInd[ii] && hashValue[ind]<(nList)){
				hashValue[ind]=hashValue[ind]+nList;
				/*printf("%i %i %i \n",listInd[ii],hashKey[ind],hashValue[ind]);*/	
			}
			ind=(ind+1)%nHash;

			kk++;
		}
		/*printf("%i\n",kk);*/
		jjErr=(kk>jjErr)*kk+(kk<=jjErr)*jjErr;
		hashKey[ind]=(int)listInd[ii];		
		hashValue[ind]=ii;
		/*printf("%i %i %i \n",listInd[ii],hashKey[ind],hashValue[ind]);*/
		
	}
	return(jjErr);

}

int ReadFindObjHashTable(int **dest,double *lookup, double *listInd, int nLook, int nList,
	int *hashKey,int *hashValue,int nHash, int *nMatch, int *nFound){

	int ii,jjErr,ind,kk,ll,jj,kkStart,isFound;

	ll=0;
	jjErr=0;
	for(ii=0;ii<nLook;ii++){
		(*dest)[ll]=-1;
		jj=0;
		kkStart=ll;
		if ((*nMatch - ll)<nList){
			*nMatch=*nMatch+nList+nLook;
			*dest=realloc(*dest,(*nMatch)*sizeof(int));
			if (*dest==NULL){
				return(-1);
			}
		}

		ind=(int)lookup[ii]%nHash;
		kk=0;
		while((hashKey[ind]!=(int)lookup[ii] || hashValue[ind]>=nList) && hashValue[ind]!=-1 && kk<=nHash){
			if(hashKey[ind]==(int)lookup[ii]){
				(*dest)[ll]=hashValue[ind]%nList;
				ll++;
			}
			ind=(ind+1)%nHash;
			kk++;
		}
		/*printf("%i\n",kk);*/
		isFound=kkStart!=ll;
		jjErr=(kk>jjErr)*kk+(kk<=jjErr)*jjErr;
		(*dest)[ll]=hashValue[ind];

		ll=ll+1;
	}
	*nFound=ll;
	return(jjErr);


}

int FindObjNumHash(double *lookup, double *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
	int ii,jj,kk,ll,kkStart,nHash,ind,jjErr,nSaved,isFound;
	int *hashKey;
	int *hashValue;
	double hashSize,hashMult;

	if (nList>0){
		hashSize=3;
		nHash=(int)ceil(hashSize*nList);
		hashKey=(int*)calloc(nHash,sizeof(int));
		hashValue=(int*)malloc(nHash*sizeof(int));
		for(ii=0;ii<nHash;ii++){hashValue[ii]=-1;}
		if (!nList || !nHash){
			printf("Exiting nList=%i nHash=%i \n",nList,nHash);
			return(1);
		}

		jjErr=BuildFindObjHashTable(lookup, listInd, nLook, nList,hashKey,hashValue,nHash);
		
		if(jjErr>nHash){
			printf("Hashed FindObjNum Failed; using other version\n");
			return(1);
		}

		jjErr=ReadFindObjHashTable(dest,lookup, listInd, nLook, nList,hashKey,
			hashValue,nHash,nMatch, nFound);
		

		if(jjErr>nHash){
			printf("Hashed FindObjNum Failed; using other version");
			return(1);
		}
		
		free(hashKey);
		free(hashValue);
	} else {
		for (ii=0;ii<nLook;ii++){
			(*dest)[ii]=-1;
		}
		*nFound=ii;
	}

	return(0);
}

void FindObjNum_Deprecated(double *lookup, double *listInd, int *dest, int nLook, int nList){
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
