#include <stdio.h>
#ifndef TEMPLATEGRID_H_INCLUDED
#define TEMPLATEGRID_H_INCLUDED
// Type Definitions
typedef struct {
	int index;
	double fill;
	int refineLvl;
	int *refineVec; // This has to be allocated 
	
} cellTemplate;

typedef struct {
	int index;
	int cellind[2];
	int vertex[2];
	
} edgeTemplate;

typedef struct {
	int index;
	double coord[2]; // needs to be changed wiht dim
	
} vertexTemplate;

// Global Extern Variables
extern int dim;

// Function prototypes
void Allocatecelldatstruct();
void AllocateGridStruct(int domSize[dim],
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
void DataIn();
void DataArrayToStruct();
int edgsub(int I, int J, int l, int domSize[dim]);
int vertsub(int I, int J, int domSize[dim]);
int cellsub(int I, int J, int domSize[dim]);
void EdgeIJtoGrid(int domSize[dim],int IJK[dim], int l);
void EdgeIJtoGridLowBound(int domSize[dim],int IJK[dim], int l);
void EdgeIJtoGridHighBound(int domSize[dim],int IJK[dim], int l);
void VertexIJGrid(int domSize[dim],int IJK[dim]);
void CellIJGrid(int domSize[dim],int IJK[dim], int baseRefineLvl);
void AssignEdgestructContent(int domSize[dim]);
void AssignVertextructContent(int domSize[dim]);
void AssignCelltructContent(int domSize[dim], int baseRefineLvl);
void BuildLvlTemplate(int domSize[dim], int baseRefineLvl, int nLevelsInput,
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
//void OutputTemplateGrid(int domSize[dim], int lvlGrid);
//void GenerateTemplateGrids();

#endif

