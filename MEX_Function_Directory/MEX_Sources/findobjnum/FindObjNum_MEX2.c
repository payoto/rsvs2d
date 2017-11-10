#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "mex.h" /* Always include this */
#include "matrix.h"
//#include "llist.h"
/* Constant declaration */
/*int dim() = 2; */
#define MEX_COMPILE

int FindObjNum(int *lookup, int *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound);
void FindObjNum_Deprecated(int *lookup, int *listInd, int *dest, int nLook, int nList);
int FindObjNumHash(int *lookup, int *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound);
int BuildFindObjHashTable(int *lookup, int *listInd, int nLook, int nList,int *hashKey,int *hashValue, int nHash);

void mexFunction(int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray *prhs[]){
	
	int *lookup, *listInd,*dest;
	int nList, nLook,ii,nFound,nMatch,isErr;
	double *outPtr;
	double *inputArrayPtr;
	/* Read Integers */
	nLook=(int)mxGetScalar(prhs[2]);
	nList=(int)mxGetScalar(prhs[3]);
	
	/*if(nList>0){*/
		/* Allocation */
	nMatch=(nLook+(nList+1))*2;
	lookup=(int*)calloc(nLook,sizeof(int));
	dest=(int*)calloc(nMatch,sizeof(int));
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
		
	free(lookup);
	free(listInd);
	free(dest);
	/*} else {plhs[0]=(mxCreateDoubleMatrix(1,0,mxREAL));}*/
	
	return;
}


int FindObjNum(int *lookup, int *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
	int ii,jj,kk,kkStart;
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
			
			(*dest)[kk]=(-(listInd[jj]!=lookup[ii])*(jj+1))+jj;
			/*dest[ii]=(-(memcmp((listInd+jj),(lookup+ii),sizeof(int))!=0)*(jj+1))+jj; */
			kk=kk+(listInd[jj]==lookup[ii]);
			jj++;
		} while((jj<nList));
		kk=kk+(kkStart==kk);
	}
	*nFound=kk;
}

int BuildFindObjHashTable(int *lookup, int *listInd, int nLook, int nList,int *hashKey,int *hashValue,int nHash){

	int ii,jjErr,ind,kk;

	jjErr=0;
	for(ii=0;ii<nList;ii++){
		ind=listInd[ii]%nHash;
		kk=0;
		while(hashValue[ind]!=-1  && kk<=nHash){
			if (hashKey[ind]==listInd[ii] && hashValue[ind]<(nList)){
				hashValue[ind]=hashValue[ind]+nList;
				//printf("%i %i %i \n",listInd[ii],hashKey[ind],hashValue[ind]);
			}
			ind=(ind+1)%nHash;

			kk++;
		}
		//printf("%i\n",kk);
		jjErr=(kk>jjErr)*kk+(kk<=jjErr)*jjErr;
		hashKey[ind]=listInd[ii];		
		hashValue[ind]=ii;
		//printf("%i %i %i \n",listInd[ii],hashKey[ind],hashValue[ind]);
		
	}
	return(jjErr);

}


int ReadFindObjHashTable(int **dest,int *lookup, int *listInd, int nLook, int nList,
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

		ind=lookup[ii]%nHash;
		kk=0;
		while((hashKey[ind]!=lookup[ii] || hashValue[ind]>=nList) && hashValue[ind]!=-1 && kk<=nHash){
			if(hashKey[ind]==lookup[ii]){
				(*dest)[ll]=hashValue[ind]%nList;
				ll++;
			}
			ind=(ind+1)%nHash;
			kk++;
		}
		//printf("%i\n",kk);
		isFound=kkStart!=ll;
		jjErr=(kk>jjErr)*kk+(kk<=jjErr)*jjErr;
		(*dest)[ll]=hashValue[ind];

		ll=ll+1;
	}
	*nFound=ll;
	return(jjErr);


}

int FindObjNumHash(int *lookup, int *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound){
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

int FindObjNumHash2(int *lookup, int *listInd, int **dest, int nLook, int nList, int *nMatch, int *nFound){
	/* Look up is the Indices we need the position of */
	/* listInd is the list of all indices in the right order */
	/* dest is where the position will be stored (-1 means the index wasn't found) */
	/* While not the most efficient yet this should be pretty good */
	int ii,jj,kk,ll,kkStart,nHash,ind,jjErr,nSaved,isFound;
	int *hashKey;
	double *hashValue;
	double hashSize,hashMult;

	hashSize=3;
	nHash=(int)ceil(hashSize*nList);
	hashKey=(int*)calloc(nHash,sizeof(int));
	hashValue=(double*)calloc(nHash,sizeof(double));
	//for(ii=0;ii<nHash;ii++){hashValue[ii]=-1;}

	jjErr=0;
	for(ii=0;ii<nList;ii++){
		ind=listInd[ii]%nHash;
		kk=0;
		while(hashValue[ind]!=0 && hashKey[ind]!=listInd[ii] && kk<=nHash){
			ind=(ind+1)%nHash;
			kk++;
		}
		jjErr=(kk>jjErr)*kk+(kk<=jjErr)*jjErr;
		hashKey[ind]=listInd[ii];
		// Hash value uses a trick to store multiple values in a single double slot
		hashMult=1;
		if (hashValue[ind]!=0){
			hashMult=pow(nList+1,ceil(log(hashValue[ind])/log((double)(nList+1))));
		}
		hashValue[ind]=hashValue[ind]+hashMult*(ii+1);
		jjErr=jjErr+(nHash+2)*(hashMult>(int)pow(10,30));
	}

	if(jjErr>nHash){printf("Hashed FindObjNum Failed; using other version\n");return(1);}

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

		ind=lookup[ii]%nHash;
		kk=0;
		while(hashKey[ind]!=lookup[ii] && hashValue[ind]!=0 && kk<=nHash){
			ind=(ind+1)%nHash;
			kk++;
		}
		isFound=hashValue[ind]!=0;
		jjErr=(kk>jjErr)*kk+(kk<=jjErr)*jjErr;

		if (isFound){
			nSaved=(int)floor(log(hashValue[ind])/log((double)(nList+1)))+1;
			//(*dest)[ll]=((int)(hashValue[ind])%(nList+1))-1;ll++;
			for(jj=0;jj<nSaved;jj++){
				(*dest)[ll]=(int)fmod(floor(hashValue[ind]/pow((double)nList+1.0,(double)jj)),(double)(nList+1))-1;
				jjErr=jjErr*((*dest)[ll]>-2)+(nHash+1)*((*dest)[ll]<-2);
				ll++;
			}

		}

		ll=ll+(!isFound);
	}

	if(jjErr>nHash){printf("Hashed FindObjNum Failed; using other version");return(1);}
	*nFound=ll;
	return(0);
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
