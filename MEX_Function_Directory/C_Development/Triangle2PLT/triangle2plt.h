#include <stdio.h>
#include <stdlib.h>
#ifndef TRIANGLETOPLT_H_INCLUDED
#define TRIANGLETOPLT_H_INCLUDED

	int PerformChecks(int argc, char *argv[] );
	int OutputFelineSeg(char *argv[]);
	int OutputGridUns(char *argv[]);
	void OpenInFiles(FILE *filein[], FILE **fileout, char fileext[5][9],char *argv[]);
	void BuildEdgeDataFromTriangle(FILE *filein[4],int *edgeDat, double *vertCoord,
	    double *cellVolume,int nQuant,int nPts, int nEdge, int nCell );
	void WriteGridUns(FILE* fileout, int nPts,int nEdge, int nCell, int *edgeDat, double *vertCoord , 
	    double *cellVolume);
	int BuildSu2(char *argv[]);
	void ShowHashDistrib(int tDistrib, double *hashDistrib, int nbar);
	double HashEdge(int *edgeDat,int cEdPos,int nCell, int nPts);

#endif