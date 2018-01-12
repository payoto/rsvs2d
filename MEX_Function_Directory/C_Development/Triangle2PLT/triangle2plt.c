#include <stdio.h>
#include <string.h>
#include <math.h>
#include "triangle2plt.h"



int main ( int argc, char *argv[] )
{
    int errorMsg;
    errorMsg=PerformChecks(argc,argv);

    if (errorMsg){
        return(-1);
    }
    if (!strcmp(argv[2],"griduns")){
        errorMsg=OutputGridUns(argv);
    } else if (!strcmp(argv[2],"plt")){
        errorMsg=OutputFelineSeg(argv);
    } else if (!strcmp(argv[2],"su2")){
        errorMsg=BuildSu2(argv);
    } else if (!strcmp(argv[2],"su2e")){
        errorMsg=BuildSu2Edge(argv);
    }
    if (errorMsg){
        return(-1);
    }
    return(0);
}


int PerformChecks(int argc, char *argv[] ){

    char fileext[5][7];
    char fullfilename[400];
    int ii;
    FILE *file;

    strcpy(fileext[0], ".node");
    strcpy(fileext[1], ".ele");
    strcpy(fileext[2], ".poly");
    strcpy(fileext[3], ".edge");
    strcpy(fileext[4], ".neigh");

    if ( argc != 3 ) /* argc should be 2 for correct execution */{
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s triangle_mesh_name outtype", argv[0] );
    }
    else {
        // Check that Each required file can be opened succesfully to build the plt file
        for (ii=0;ii<5;ii++){
            strcpy(fullfilename,argv[1]);
            strcat(fullfilename,fileext[ii]);
            file = fopen( fullfilename, "r" );

            /* fopen returns 0, the NULL pointer, on failure */
            if ( file == 0 ){
                printf( "Could not open file: %s \n",fullfilename);
                return(1);
            }
            else {
                fclose( file );
            }
        }
    }


    return(0);
}


int OutputFelineSeg(char *argv[])
{

    int ii,jj;
    FILE *filein[2];
    FILE *fileout;
    char fileext[3][7];
    char fullfilename[400];
    int nPts,nSeg,nQuant,discardDat,nDisc;
    double dataPt;

    strcpy(fileext[0], ".node");
    strcpy(fileext[1], ".edge");
    strcpy(fileext[2], ".plt");
    for (ii=0;ii<2;ii++){
        strcpy(fullfilename,argv[1]);
        strcat(fullfilename,fileext[ii]);
        filein[ii]= fopen( fullfilename, "r" );
    }
    strcpy(fullfilename,argv[1]);
    strcat(fullfilename,fileext[2]);
    fileout= fopen( fullfilename, "w" );


    if ( fileout != NULL ){
        fscanf(filein[0]," %i %i %i ",&nPts,&discardDat,&nQuant);
        nDisc=1;for (jj = 0; jj < nDisc; ++jj){fscanf(filein[0],"%i",&discardDat);}

        fscanf(filein[1],"%i %i ",&nSeg,&discardDat);

        // Print Header to FELINESEG file
        fprintf(fileout, "VARIABLES = \"X\" ,\"Y\"");
        for(ii=0;ii<nQuant;ii++){
            fprintf(fileout, " ,\"var%i\" ",ii);
        }
        fprintf(fileout, "\nZONE \nVARLOCATION=([1-%i]=NODAL)",nQuant+2);
        fprintf(fileout, "\nNODES=%i \nELEMENTS=%i",nPts,nSeg);
        fprintf(fileout, "\nDATAPACKING=POINT \nZONETYPE=FELINESEG \n");


        // Print Vertex Data
        for(ii=0;ii<nPts;ii++){
            nDisc=1;for (jj = 0; jj < nDisc; ++jj){fscanf(filein[0],"%i",&discardDat);}
            for (jj = 0; jj < (2+nQuant); ++jj){
                fscanf(filein[0],"%lf",&dataPt);
                fprintf(fileout,"%.16e ",dataPt);
            }
            nDisc=1;for (jj = 0; jj < nDisc; ++jj){fscanf(filein[0],"%i",&discardDat);}
            fprintf(fileout,"\n");
        }
        // Print Segment Data
        for(ii=0;ii<nSeg;ii++){
            nDisc=1;for (jj = 0; jj < nDisc; ++jj){fscanf(filein[1],"%i",&discardDat);}
            for (jj = 0; jj < (2); ++jj){
                fscanf(filein[1],"%i",&discardDat);
                fprintf(fileout,"%i ",discardDat);
            	if(ii==38933 || ii==39263 || ii==39278){printf("%i ",discardDat);}
            }
            
            nDisc=1;for (jj = 0; jj < nDisc; ++jj){fscanf(filein[1],"%i",&discardDat);}
            fprintf(fileout,"\n");
            if(ii==38933 || ii==39263 || ii==39278){printf("      %i\n ",ii);}
        }
        fclose(filein[0]);
        fclose(filein[1]);
        fclose(fileout);
    }
    else{
        printf( "Could not open file: %s \n",fullfilename);
        return(1);
    }
    return(0);
}


int OutputGridUns(char *argv[]){

    int ii,jj,kk;
    FILE *filein[5];
    FILE *fileout;
    char fileext[5][9];

    int nPts,nCell,nEdge,nQuant,discardInt,nDisc;
    double dataPt,discardDouble;
    int *edgeDat=NULL,*vertBound=NULL;
    double *edgeVal=NULL,*vertCoord=NULL,*cellVolume=NULL;

    strcpy(fileext[4], ".griduns");
    OpenInFiles(filein, &fileout, fileext,argv);

    // Find the sizes of arraysfileext
    fscanf(filein[0],"%i %i %i ",&nPts,&discardInt,&nQuant);
    fscanf(filein[0],"%i  ",&discardInt);
    fscanf(filein[1],"%i %i ",&nEdge,&discardInt);
    fscanf(filein[2],"%i %i %i ",&nCell,&discardInt,&nQuant);
    fscanf(filein[3],"%i  ",&discardInt);
    fscanf(filein[3],"%i  ",&discardInt);
    printf("nVert: %i ; nEdge: %i ; nCell: %i\n",nPts,nEdge,nCell);

    edgeDat=(int*)malloc(4*(nEdge+1)*sizeof(int));
    vertBound=(int*)malloc(nPts*sizeof(int));

    vertCoord=(double*)malloc(nPts*2*sizeof(double));
    cellVolume=(double*)calloc(nCell,sizeof(double));

    for (ii=0;ii<nPts;ii++){
        fscanf(filein[0],"%i %lf %lf %i ",&discardInt, &vertCoord[2*ii],&vertCoord[2*ii+1],&vertBound[ii]);
    }

    BuildEdgeDataFromTriangle(filein,edgeDat,vertCoord,cellVolume,nQuant,nPts,nEdge,nCell);

    for (ii=0;ii<nEdge;ii++){
        //if ((edgeDat[4*ii]-1)>nPts){printf("%6i %6i %6i \n",ii,edgeDat[4*ii]-1,nPts);}

        edgeDat[4*ii+3]=(edgeDat[4*ii+3]<0)*vertBound[edgeDat[4*ii]-1]
            +(edgeDat[4*ii+3]>=0)*edgeDat[4*ii+3];
    }
    WriteGridUns(fileout, nPts,nEdge, nCell, edgeDat, vertCoord ,cellVolume);
    // Exit cleanly
    free(edgeDat);
    free(vertBound);
    free(vertCoord);
    free(cellVolume);
    for (int i = 0; i < 4; ++i)
    {
        fclose(filein[i]);
    }
    fclose(fileout);
    return(0);
}

void BuildEdgeDataFromTriangle(FILE *filein[4],int *edgeDat,double *vertCoord,
    double *cellVolume,int nQuant,int nPts, int nEdge, int nCell ){

    int ii,kk,discardInt,jj,jjj,nDisc,cEdPosTest;
    int vertInd[3],cellInd[3],currCell,cEdPos,cEdSaved,flagInc,edgeHashInd,nConflicts,flagMatch,nEdgeHash,tDistrib;
    double discardDouble,x1,y1,x2,y2,edgeHash,nHashRatio,saveHash;
    double *edgeVal=NULL,*hashDistrib=NULL;
    int *edgePos=NULL;
    nConflicts=0;
    cEdSaved=0;
    nHashRatio=2;
    tDistrib=20+2;
    nEdgeHash=ceil((nEdge+nCell)*nHashRatio);
    edgeVal=(double*)calloc(nEdgeHash+1,sizeof(double));
    edgePos=(int*)calloc(nEdgeHash+1,sizeof(int));
    hashDistrib=(double*)calloc(tDistrib,sizeof(double));

    for (ii = 0; ii< nCell; ++ii)
    {   
        fscanf(filein[2]," %i ",&currCell);
        fscanf(filein[3]," %i ",&discardInt);
        for(jj=0;jj<3;jj++){
            fscanf(filein[2]," %i ",&vertInd[jj]);
            fscanf(filein[3]," %i ",&cellInd[jj]);
        }

        nDisc=nQuant;for (jj = 0; jj < nDisc; ++jj){fscanf(filein[2]," %lf ",&discardDouble);}

        for(jj=0;jj<3;jj++){
            cEdPos=cEdSaved*4-1; // -1 is because i'm silly and forgot about 0 indexing
            edgeDat[cEdPos+1]=vertInd[jj];
            edgeDat[cEdPos+2]=vertInd[(jj+1)%3];
            edgeDat[cEdPos+3]=currCell;
            edgeDat[cEdPos+4]=cellInd[(3+jj-1)%3];

            
            // Calculate cell volume with greens theorem
            x1=vertCoord[(vertInd[jj]-1)*2];
            y1=vertCoord[(vertInd[jj]-1)*2+1];
            x2=vertCoord[(vertInd[(jj+1)%3]-1)*2];
            y2=vertCoord[(vertInd[(jj+1)%3]-1)*2+1];
            cellVolume[ii]=cellVolume[ii]+((x1+x2)*(y2-y1)-(y1+y2)*(x2-x1))*0.5*0.5;
            //if (ii==0){printf("%i %i %i \n",jj,(jj+1)%3,(3+jj-1)%3);}
            // Calculation of the edges hash
            //edgeVal[cEdSaved]
            // Hash table implementation
            edgeHash=HashEdge(edgeDat,cEdPos,nCell, nPts);

            hashDistrib[(int)floor(1+edgeHash*(tDistrib-2))]++;
            edgeHashInd=ceil((edgeHash)*nEdgeHash);

            // This tests if an edge with the same hash has been saved
            flagInc=1;kk=0;flagMatch=0;
            while (kk<nEdgeHash && flagInc && !flagMatch){
                if ((edgeVal[(kk+edgeHashInd)%nEdgeHash]!=0) && (edgeVal[(kk+edgeHashInd)%nEdgeHash]==edgeHash)){
                	cEdPosTest=(edgePos[(kk+edgeHashInd)%nEdgeHash])*4-1;
                	flagMatch=((edgeDat[cEdPos+1]==edgeDat[cEdPosTest+1]) || (edgeDat[cEdPos+1]==edgeDat[cEdPosTest+2])) &&
                			  ((edgeDat[cEdPos+2]==edgeDat[cEdPosTest+1]) || (edgeDat[cEdPos+2]==edgeDat[cEdPosTest+2])) &&
                			  ((edgeDat[cEdPos+3]==edgeDat[cEdPosTest+3]) || (edgeDat[cEdPos+3]==edgeDat[cEdPosTest+4])) &&
                			  ((edgeDat[cEdPos+4]==edgeDat[cEdPosTest+3]) || (edgeDat[cEdPos+4]==edgeDat[cEdPosTest+4]));
                }

                flagInc=(!flagMatch) 
                    && (edgeVal[(kk+edgeHashInd)%nEdgeHash]!=0);
                kk++;
            }
            nConflicts=nConflicts+kk;
            /*
            if ((cEdSaved%50)==0)
                printf("%20.10lf %10i %10i %10i %10i \n",edgeHash,edgeHashInd,kk-1,nConflicts,cEdSaved);*/
            /*
            if (cEdSaved<50){
                printf("%7i %7i %7i %7i %20.10lf ",edgeDat[cEdPos+1],edgeDat[cEdPos+2],
                    edgeDat[cEdPos+3],edgeDat[cEdPos+4],edgeVal[cEdSaved]);
                printf("%6i %6i \n",cEdSaved,flagInc );
            }*/
   //          if(ii==24063 || ii==24050 || ii==23829 || ii==24066){
	  //           for (jjj=1;jjj<5;++jjj){
			// 		printf("%i ",edgeDat[cEdPos+jjj]);
			// 	}
			// 	printf("%i   %i    %i    %i   %i    %i    %li\n",cEdPos,ii,flagMatch,cEdSaved,
			// 		edgeHashInd,(kk+edgeHashInd)%nEdgeHash,edgeHash);
			// }
			// if(edgeHashInd==56049){
	  //           for (jjj=1;jjj<5;++jjj){
			// 		printf("%8i ",edgeDat[cEdPos+jjj]);
			// 	}
			// 	printf("%8i %8i %8i %8i %8i %8i %.20lf\n",cEdPos,ii,flagMatch,cEdSaved,
			// 		edgeHashInd,(kk+edgeHashInd)%nEdgeHash,edgeHash);
			// }
            if (!flagMatch){
                edgeVal[(kk-1+edgeHashInd)%nEdgeHash]=edgeHash;
                edgePos[(kk-1+edgeHashInd)%nEdgeHash]=cEdSaved;
                cEdSaved=cEdSaved+1;
            }
        }

    }

    printf("%20.10lf %10i %10i %10i %10i \n",edgeHash,edgeHashInd,kk-1,nConflicts,cEdSaved);
    ShowHashDistrib( tDistrib, hashDistrib, 50);
    printf("Hash Efficiency : O(%.2lf)\n",((double)nConflicts/(double)cEdSaved));
    for (ii=0;ii<40000;++ii){
		if(ii==38933 || ii==39263 || ii==39278){
			for (jj=0;jj<4;++jj){
				printf("%i ",edgeDat[4*ii+jj]);
			}
			printf("%i \n",ii);
		}
	}
    free(edgeVal);
    free(hashDistrib);
}

double HashEdge(int *edgeDat,int cEdPos,int nCell, int nPts){
    //  Returns a Hash between 0 and 1  depending on the properties of the edgeDat

    double edgeHash;

    // edgeHash=2.0*(
    //     (double)(abs((edgeDat[cEdPos+1]<edgeDat[cEdPos+2])*edgeDat[cEdPos+1]+
    //      (edgeDat[cEdPos+1]>edgeDat[cEdPos+2])*edgeDat[cEdPos+2]))
    //     *(1.0/(double)(nPts+1)/(double)(nPts+1)/(double)(nCell+1)/(double)(nCell+1))+
    //     (double)(abs((edgeDat[cEdPos+1]>edgeDat[cEdPos+2])*edgeDat[cEdPos+1]+
    //      (edgeDat[cEdPos+1]<edgeDat[cEdPos+2])*edgeDat[cEdPos+2]))
    //     *((double)1.0/((double)(nPts+1)*(double)(nCell+1)*(double)(nCell+1)))+
    //     (double)(abs((edgeDat[cEdPos+3]<edgeDat[cEdPos+4])*edgeDat[cEdPos+3]+
    //      (edgeDat[cEdPos+3]>edgeDat[cEdPos+4])*edgeDat[cEdPos+4]))
    //     *((double)1.0/((double)(nCell+1)*(double)(nCell+1)))+
    //     (double)(abs((edgeDat[cEdPos+3]>edgeDat[cEdPos+4])*edgeDat[cEdPos+3]+
    //      (edgeDat[cEdPos+3]<edgeDat[cEdPos+4])*edgeDat[cEdPos+4]))
    //     *((double)(1.0)/((nCell+1))));
    
    // edgeHash=2.0*(((((((((
    //     ((double)(abs((edgeDat[cEdPos+1]<edgeDat[cEdPos+2])*edgeDat[cEdPos+1]+ // min vert value
    //      (edgeDat[cEdPos+1]>edgeDat[cEdPos+2])*edgeDat[cEdPos+2]))))/(double)(nPts+1))+

    //     ((double)(abs((edgeDat[cEdPos+1]>edgeDat[cEdPos+2])*edgeDat[cEdPos+1]+ // max vert value
    //      (edgeDat[cEdPos+1]<edgeDat[cEdPos+2])*edgeDat[cEdPos+2]))))/(double)(nPts+1))+

    //     ((double)(abs((edgeDat[cEdPos+3]<edgeDat[cEdPos+4])*edgeDat[cEdPos+3]+
    //      (edgeDat[cEdPos+3]>edgeDat[cEdPos+4])*edgeDat[cEdPos+4]))))/(double)(nCell+1))+

    //     ((double)(abs((edgeDat[cEdPos+3]>edgeDat[cEdPos+4])*edgeDat[cEdPos+3]+
    //      (edgeDat[cEdPos+3]<edgeDat[cEdPos+4])*edgeDat[cEdPos+4])))))/(double)(nCell+1));

	edgeHash=
        ((double)(((edgeDat[cEdPos+1]<edgeDat[cEdPos+2])*edgeDat[cEdPos+1]+ // min vert value
         (edgeDat[cEdPos+1]>edgeDat[cEdPos+2])*edgeDat[cEdPos+2])))/(double)(nPts+1);

    edgeHash=(edgeHash+(double)(((edgeDat[cEdPos+1]>edgeDat[cEdPos+2])*edgeDat[cEdPos+1]+ // max vert value
         (edgeDat[cEdPos+1]<edgeDat[cEdPos+2])*edgeDat[cEdPos+2])))/(double)(nPts+1);

    edgeHash=(edgeHash+(double)(((edgeDat[cEdPos+3]<edgeDat[cEdPos+4])*edgeDat[cEdPos+3]+
         (edgeDat[cEdPos+3]>edgeDat[cEdPos+4])*edgeDat[cEdPos+4])))/(double)(nCell+1);

    edgeHash=(edgeHash+(double)(((edgeDat[cEdPos+3]>edgeDat[cEdPos+4])*edgeDat[cEdPos+3]+
         (edgeDat[cEdPos+3]<edgeDat[cEdPos+4])*edgeDat[cEdPos+4])))/(double)(nCell+1);
	edgeHash=edgeHash;
    edgeHash=2.0*(edgeHash)*(edgeHash<=0.5)+2.0*(edgeHash>0.5)*(edgeHash-0.5);
            //edgeHash=(edgeHash*0.5);
    return(edgeHash);

}

void ShowHashDistrib(int tDistrib, double *hashDistrib, int nbar){
    int ii,jj;
    double factorReduce;
    printf("\nPRINT HASH DISTRIBUTION %i\n\n",tDistrib);
    for(ii=0;ii<tDistrib;ii++){
        factorReduce=factorReduce*(factorReduce>hashDistrib[ii])+hashDistrib[ii]*(factorReduce<=hashDistrib[ii]);
    }

    for(ii=0;ii<tDistrib;ii++){
        printf("%3i %10.0lf : ",ii,hashDistrib[ii]);
        for(jj=0;jj< (hashDistrib[ii]*nbar/factorReduce);jj++){
            printf("|");
        }
        printf("\n");
    }
    printf("\n\n");

}

void WriteGridUns(FILE* fileout, int nPts,int nEdge, int nCell, int *edgeDat, double *vertCoord , 
    double *cellVolume){
    int ii,jj,nRep;

    fprintf(fileout, "%i %i %i\n",nCell,nEdge,nPts );
    nRep=4;
    for (ii = 0; ii < nEdge; ++ii){
        for (jj= 0; jj < nRep; ++jj){
            fprintf(fileout,"%i ",edgeDat[nRep*ii+jj] );
        }
        fprintf(fileout,"\n");
    }

    nRep=2;
    for (ii = 0; ii < nPts; ++ii){
        fprintf(fileout,"%i ",ii+1);
        for (jj= 0; jj < nRep; ++jj){
            fprintf(fileout,"%.25e ",vertCoord[nRep*ii+jj] );
        }
        fprintf(fileout,"\n");
    }
    for (ii = 0; ii < nCell; ++ii){
        fprintf(fileout,"%.25e ",cellVolume[ii]);
        fprintf(fileout,"\n");
    }
}

void OpenInFiles(FILE *filein[5], FILE **fileout, char fileext[5][9],char *argv[]){

    int ii;
    char fullfilename[400];

    // Opening all the required input and output files
    strcpy(fileext[0], ".node");
    strcpy(fileext[1], ".edge");
    strcpy(fileext[2], ".ele");
    strcpy(fileext[3], ".neigh");

    for (ii=0;ii<4;ii++){
        strcpy(fullfilename,argv[1]);
        strcat(fullfilename,fileext[ii]);
        filein[ii]= fopen( fullfilename, "r" );
    }
    strcpy(fullfilename,argv[1]);
    strcat(fullfilename,fileext[ii]);
    *fileout= fopen( fullfilename, "w" );

}

int BuildSu2(char *argv[]){

    int ii,jj,kk;
    FILE *filein[5];
    FILE *fileout;
    char fileext[5][9];

    int nPts,nCell,nEdge,nQuant,discardInt,nDisc,nBound,currPos,flagMatch;
    double dataPt,discardDouble;
    int *edgeDat=NULL,*vertBound=NULL,boundType[10],nElmBound[10];
    double *edgeVal=NULL,*vertCoord=NULL,*cellVolume=NULL;

    strcpy(fileext[4], ".su2");
    OpenInFiles(filein, &fileout, fileext,argv);

    // Find the sizes of arraysfileext
    fscanf(filein[0],"%i %i %i ",&nPts,&discardInt,&nQuant);
    fscanf(filein[0],"%i  ",&discardInt);
    fscanf(filein[1],"%i %i ",&nEdge,&discardInt);
    fscanf(filein[2],"%i %i %i ",&nCell,&discardInt,&nQuant);
    fscanf(filein[3],"%i  ",&discardInt);
    fscanf(filein[3],"%i  ",&discardInt);
    printf("nVert: %i ; nEdge: %i ; nCell: %i\n",nPts,nEdge,nCell);

    fprintf(fileout,"NDIME= 2\n");
    fprintf(fileout,"NELEM= %i\n",nCell);

    // Write triangular cell Data
    for(ii=0;ii<nCell;ii++){
        fprintf(fileout," 5 ");
        fscanf(filein[2],"%i ",&discardInt);
        for (jj=0;jj<3;jj++){
            fscanf(filein[2],"%i ",&discardInt);
            fprintf(fileout," %i ",discardInt-1);
        }
        fprintf(fileout," %i\n",ii);
    }
    
    // Write Point Data
    fprintf(fileout,"NPOIN= %i\n",nPts);
    for(ii=0;ii<nPts;ii++){
        fscanf(filein[0],"%i ",&discardInt);
        for (jj=0;jj<2;jj++){
            fscanf(filein[0],"%lf ",&discardDouble);
            fprintf(fileout," %lf ",discardDouble);
        }
        fscanf(filein[0],"%i ",&discardInt);
        fprintf(fileout," %i\n",ii);
    }

    // Read Boundary regions
    edgeDat=(int*) malloc(3*nEdge*sizeof(int));
    kk=0;
    nBound=0;
    printf("nVert: %i ; nEdge: %i ; nCell: %i\n",nPts,nEdge,nCell);

    for (ii=0;ii<nEdge;ii++){
        // Read Data from file
        fscanf(filein[1]," %i ",&discardInt);
        currPos=kk*3;
        for(jj=0;jj<3;jj++){
            fscanf(filein[1]," %i ",&edgeDat[currPos+jj]);
            
        }
        kk=kk+(edgeDat[currPos+2]<0);

        // MAtch to previous recognised boundaries
        jj=0;
        flagMatch=0;
        
        while(edgeDat[currPos+2]<0 && jj<nBound && jj<10 && !flagMatch){
            flagMatch=boundType[jj]==edgeDat[currPos+2];
            jj++;
        }
        /*
        if (edgeDat[currPos+2]<0){
            for (int i = 0; i < 3; ++i){
                printf("%i ",edgeDat[currPos+i]);
            }
            printf("%i %i ",kk,jj,flagMatch);
            printf("\n");
        }*/
        if (!flagMatch && edgeDat[currPos+2]<0){
            boundType[nBound]=edgeDat[currPos+2];
            nElmBound[nBound]=1;
            nBound++;
        } else if (edgeDat[currPos+2]<0) {
            nElmBound[jj-1]++;
        }
    }
    // Write Boundary COnditions
    
    fprintf(fileout,"NMARK= %i\n",nBound);
    for(ii=0;ii<nBound;ii++){

        fprintf(fileout, "MARKER_TAG= %i\n", boundType[ii]);
        fprintf(fileout, "MARKER_ELEMS= %i\n", nElmBound[ii]);

        for (jj=0;jj<kk;jj++){
            if(edgeDat[jj*3+2]==boundType[ii]){
                fprintf(fileout, "3 %i %i\n",edgeDat[jj*3+0]-1,edgeDat[jj*3+1]-1);
            }
        }
    }
    
    // Clean Exit
    free(edgeDat);
    for (int i = 0; i < 4; ++i)
    {
        fclose(filein[i]);
    }
    fclose(fileout);
    return(0);
}

int BuildSu2Edge(char *argv[]){

    int ii,jj,kk;
    FILE *filein[5];
    FILE *fileout;
    char fileext[5][9];

    int nPts,nCell,nEdge,nQuant,discardInt,nDisc,nBound,currPos,flagMatch;
    double dataPt,discardDouble;
    int *edgeDat=NULL,*vertBound=NULL,boundType[10],nElmBound[10];
    double *edgeVal=NULL,*vertCoord=NULL,*cellVolume=NULL;

    strcpy(fileext[4], ".e.su2");
    OpenInFiles(filein, &fileout, fileext,argv);

    // Find the sizes of arraysfileext
    fscanf(filein[0],"%i %i %i ",&nPts,&discardInt,&nQuant);
    fscanf(filein[0],"%i  ",&discardInt);
    fscanf(filein[1],"%i %i ",&nEdge,&discardInt);
    fscanf(filein[2],"%i %i %i ",&nCell,&discardInt,&nQuant);
    fscanf(filein[3],"%i  ",&discardInt);
    fscanf(filein[3],"%i  ",&discardInt);
    printf("nVert: %i ; nEdge: %i ; nCell: %i\n",nPts,nEdge,nCell);

    fprintf(fileout,"NDIME= 2\n");
    fprintf(fileout,"NELEM= %i\n",nEdge);

    // Write Edge date for cell Data
    for(ii=0;ii<nEdge;ii++){
        fprintf(fileout," 3 ");
        fscanf(filein[1],"%i ",&discardInt);
        for (jj=0;jj<2;jj++){
            fscanf(filein[1],"%i ",&discardInt);
            fprintf(fileout," %i ",discardInt-1);
        }
        fscanf(filein[1],"%i ",&discardInt);
        fprintf(fileout," %i\n",ii);
    }
    fseek(filein[1],0,SEEK_SET);
    fscanf(filein[1],"%i %i ",&nEdge,&discardInt);
    
    // Write Point Data
    fprintf(fileout,"NPOIN= %i\n",nPts);
    for(ii=0;ii<nPts;ii++){
        fscanf(filein[0],"%i ",&discardInt);
        for (jj=0;jj<2;jj++){
            fscanf(filein[0],"%lf ",&discardDouble);
            fprintf(fileout," %lf ",discardDouble);
        }
        fscanf(filein[0],"%i ",&discardInt);
        fprintf(fileout," %i\n",ii);
    }

    // Read Boundary regions
    edgeDat=(int*) malloc(3*nEdge*sizeof(int));
    kk=0;
    nBound=0;
    printf("nVert: %i ; nEdge: %i ; nCell: %i\n",nPts,nEdge,nCell);

    for (ii=0;ii<nEdge;ii++){
        // Read Data from file
        fscanf(filein[1]," %i ",&discardInt);
        currPos=kk*3;
        for(jj=0;jj<3;jj++){
            fscanf(filein[1]," %i ",&edgeDat[currPos+jj]);
            
        }
        kk=kk+(edgeDat[currPos+2]<0);

        // MAtch to previous recognised boundaries
        jj=0;
        flagMatch=0;
        
        while(edgeDat[currPos+2]<0 && jj<nBound && jj<10 && !flagMatch){
            flagMatch=boundType[jj]==edgeDat[currPos+2];
            jj++;
        }
        /*
        if (edgeDat[currPos+2]<0){
            for (int i = 0; i < 3; ++i){
                printf("%i ",edgeDat[currPos+i]);
            }
            printf("%i %i ",kk,jj,flagMatch);
            printf("\n");
        }*/
        if (!flagMatch && edgeDat[currPos+2]<0){
            boundType[nBound]=edgeDat[currPos+2];
            nElmBound[nBound]=1;
            nBound++;
        } else if (edgeDat[currPos+2]<0) {
            nElmBound[jj-1]++;
        }
    }
    // Write Boundary COnditions
    
    fprintf(fileout,"NMARK= %i\n",nBound);
    for(ii=0;ii<nBound;ii++){

        fprintf(fileout, "MARKER_TAG= %i\n", boundType[ii]);
        fprintf(fileout, "MARKER_ELEMS= %i\n", nElmBound[ii]);

        for (jj=0;jj<kk;jj++){
            if(edgeDat[jj*3+2]==boundType[ii]){
                fprintf(fileout, "3 %i %i\n",edgeDat[jj*3+0]-1,edgeDat[jj*3+1]-1);
            }
        }
    }
    
    // Clean Exit
    free(edgeDat);
    for (int i = 0; i < 4; ++i)
    {
        fclose(filein[i]);
    }
    fclose(fileout);
    return(0);
}