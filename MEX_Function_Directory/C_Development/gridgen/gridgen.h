#include <stdio.h>
#include <stdlib.h>
#ifndef TEMPLATEGRID_H_INCLUDED
#define TEMPLATEGRID_H_INCLUDED
/* Type Definitions */
typedef struct {
	int index;
	double fill;
	int refineLvl;
	int *refineVec; /* This has to be allocated  */
	
} cellTemplate;

typedef struct {
	int index;
	int cellind[2];
	int vertex[2];
	int orientation; /* 1 is vertical 0 is Horizontal */
} edgeTemplate;

typedef struct {
	int index;
	double coord[2]; /* needs to be changed wiht dim */
	
} vertexTemplate;

/* Global Extern Variables */
/*extern int 2; */
/* Macros */
#if defined(__GNUC__) || defined(__GNUG__)
	/* GNU GCC/G++. --------------------------------------------- */
#define max(a,b) ({ __typeof__(a) _a = (a);  __typeof__(b) _b = (b);  _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a);  __typeof__(b) _b = (b);  _a < _b ? _a : _b; })

#elif defined(_MSC_VER)
	/* Microsoft Visual Studio. --------------------------------- */
	
#endif

#define dim() (2)
/* Function prototypes */
void ClearWorkSpace();
void Allocatecelldatstruct();
void AllocateGridStruct(int domSize[dim()],
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
void DataIn();
void DataArrayToStruct();
int edgsub(int I, int J, int l, int domSize[dim()]);
int vertsub(int I, int J, int domSize[dim()]);
int cellsub(int I, int J, int domSize[dim()]);
void EdgeIJtoGrid(int domSize[dim()],int IJK[dim()], int l);
void EdgeIJtoGridLowBound(int domSize[dim()],int IJK[dim()], int l);
void EdgeIJtoGridHighBound(int domSize[dim()],int IJK[dim()], int l);
void VertexIJGrid(int domSize[dim()],int IJK[dim()]);
void CellIJGrid(int domSize[dim()],int IJK[dim()], int baseRefineLvl);
void AssignEdgestructContent(int domSize[dim()]);
void AssignVertextructContent(int domSize[dim()]);
void AssignCelltructContent(int domSize[dim()], int baseRefineLvl);
void BuildLvlTemplate(int domSize[dim()], int baseRefineLvl, int nLevelsInput,
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
void CalculateNumElements(int domSize[dim()],int *nCellCurr,int *nEdgeCurr,int *nVertCurr);
void MemCopyCellStruct(cellTemplate *original, cellTemplate *destination, int nElm);
void Allocatecelldatstruct();
void InitialiseGridFromFile();
void DataIn();
void DataArrayToStruct();
void OutputTemplateGrid(int domSize[dim()], int lvlGrid);
void GenerateTemplateGrid(int lvlGenerate);
void RefineGrid();
void DeAllocate();
void DeAllocateTemplate(int domSize[dim()], cellTemplate *cellAct, edgeTemplate *edgeAct, vertexTemplate *vertAct);
void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList);
void GridInitialisation();
void OutputGrid(int jj);
int SortOrder(const void * elem1, const void * elem2);
void IdentifyRefineCell(int refinLvl,int **posGridRef,int **indGridRef, int *nRefineCell);
void IdentifyRefineEdge(int *posCellRefine, int *indCellRefine,int nCellRefine,
		int **posEdgeRefine, int **indEdgeRefine,int *nEdgeRefine, edgeTemplate *edgestructAct, int nEdge, int nEdgepCell);
void RemoveIdenticalEntries_int(int **array, int nArray, int *nNewArray);
int imin(int a, int b);
int imax(int a, int b);
double f_min(double a, double b);
double f_max(double a, double b);
double fsum(double *array, int nArray);
double fprod(double *array, int nArray);
int isum(int *array, int nArray);
int iprod(int *array, int nArray);
int min_array(int *array, int nArray);
int pos_min_array(int *array, int nArray);
int max_array(int *array, int nArray);
int pos_max_array(int *array, int nArray);
double fmin_array(double *array, int nArray);
int pos_fmin_array(double *array, int nArray);
double fmax_array(double *array, int nArray);
int pos_fmax_array(double *array, int nArray);
void RefineSelectedEdges(int domSize[dim()],int *posEdgeRefine,int *indEdgeRefine, int nEdgeRefine);
void* ArrayFromStruct(void* ptr, int nStruct, int sizMember, int sizType);
void ExtractEdgeChain(int *posEdge, int *indEdge,int nEdge ,int *ordPosEdge,
	int *ordIndEdge, int *ordIndVert,edgeTemplate *edgestructAct);
int LowerLeftVertex(double *coordList, int nEdge);
double* ExtractCoord(int nEdge, int *ordIndVert,int * ordPosVert, vertexTemplate *vertstructAct, int nVertAct);
double* ExtractAngles(double *vecTest, int nVec, int nDim);
int IsClockWiseChain(int *minVertPosRet,int nEdge, int *ordIndVert, int *ordPosVert, vertexTemplate *vertstructAct, int nVertAct);
void ReorderArray(void *oldArray, int *newOrder ,int nArray, int sizArray);
void OrderEdgeChain(int nEdge, int *ordPosEdge, int *ordIndEdge, int *ordIndVert,int *ordPosVert, vertexTemplate *vertstructAct, int nVertAct);
void GenerateIndMatch(int domSize[dim()], int posCellRefine ,int *posEdgeSideTemp,int *posVertSideTemp,int *posCellAddTemp,
	int *posEdgeAddTemp, int *posVertAddTemp, int *ordPosEdge, int *ordPosVert, int nCurrEdge,
	int **convCellList, int **convEdgeList, int **convVertList);
void ModifBorderEdges(int nCurrEdge, int posCellRefine, int *ordPosEdge, int *posEdgeSideTemp, int *convCellList);
void ExtendGridStructures(int nCurrEdge,int nCellTemplate,int nEdgeTemplate,int nVertTemplate);
void CopyEdgeStruct(int nCurrEdge,int domSize[dim()],int *posEdgeAddTemp,
	int *convEdgeList,int *convCellList,int *convVertList);
void CopyVertStruct(int nCurrEdge,int domSize[dim()],int *posVertAddTemp,
	int *convVertList, int *ordPosVert);
void CopyCellStruct(int domSize[dim()],int posCellRefine,int *convCellList);
void MergeNewCell(int domSize[dim()], int posCellRefine ,int *posEdgeSideTemp,int *posVertSideTemp,int *posCellAddTemp,
	int *posEdgeAddTemp, int *posVertAddTemp, int *ordPosEdge, int *ordPosVert, int nCurrEdge);
void RefineCell(int domSize[dim()],int posCellRefine,int indCellRefine,int *posEdgeSideTemp,
	int *posVertSideTemp,int *posCellAddTemp,int *posEdgeAddTemp, int *posVertAddTemp);
void PrepareTemplateInfo(int domSize[dim()],int **posEdgeSide,int **posVertSide,int **posCellAdd,
	int **posEdgeAdd,int **posVertAdd);
void RefineSelectedCells(int domSize[dim()],int *posCellRefine,int *indCellRefine,int nCellRefine);
int CompareDouble(double elem1, double elem2);
void AllocateSolutionStruct(int nCellCurr, int nEdgeCurr, int nVertCurr, int nLevelsAct,
		cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
void IdentifyCellEdges(int *posCellRefine, int *indCellRefine,int nCellRefine,
		int **posEdgeRefine, int **indEdgeCell,int *nEdgeCell, edgeTemplate *edgestructAct, int nEdge);
void RefineSelectedEdgesRobust(int domSize[dim()],int *posEdgeRefine,int *indEdgeRefine
	,int *posVertRefine,int *indVertRefine, int nEdgeRefine);
void AddVertexToStruct(double *newCoord, int nCoord, int *newVertInd);
void AddEdgeToStruct(int *newVertInd, int *edgeVertInd,int *edgeVertPos,int *nSplitinEdge,
	int *posEdgeRefine,int *indEdgeRefine,int *posVertRefine,int *indVertRefine,
	int nEdgeRefine,int nNewVert);

		/*void OutputTemplateGrid(int domSize[2], int lvlGrid); */
/*void GenerateTemplateGrids(); */

#endif

