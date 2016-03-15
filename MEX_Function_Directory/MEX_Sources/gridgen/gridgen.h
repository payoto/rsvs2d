#include <stdio.h>
#include <stdlib.h>
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
	int orientation; // 1 is vertical 0 is Horizontal
} edgeTemplate;

typedef struct {
	int index;
	double coord[2]; // needs to be changed wiht dim
	
} vertexTemplate;

// Global Extern Variables
//extern int 2;
// Macros
#if defined(__GNUC__) || defined(__GNUG__)
	/* GNU GCC/G++. --------------------------------------------- */
#define max(a,b) ({ typeof(a) _a = (a);  typeof(b) _b = (b);  _a > _b ? _a : _b; })
#define min(a,b) ({ typeof(a) _a = (a);  typeof(b) _b = (b);  _a < _b ? _a : _b; })

#elif defined(_MSC_VER)
	/* Microsoft Visual Studio. --------------------------------- */
	
#endif

#define dim() (2)
// Function prototypes
void Allocatecelldatstruct();
void AllocateGridStruct(int domSize[2],
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
void DataIn();
void DataArrayToStruct();
int edgsub(int I, int J, int l, int domSize[2]);
int vertsub(int I, int J, int domSize[2]);
int cellsub(int I, int J, int domSize[2]);
void EdgeIJtoGrid(int domSize[2],int IJK[2], int l);
void EdgeIJtoGridLowBound(int domSize[2],int IJK[2], int l);
void EdgeIJtoGridHighBound(int domSize[2],int IJK[2], int l);
void VertexIJGrid(int domSize[2],int IJK[2]);
void CellIJGrid(int domSize[2],int IJK[2], int baseRefineLvl);
void AssignEdgestructContent(int domSize[2]);
void AssignVertextructContent(int domSize[2]);
void AssignCelltructContent(int domSize[2], int baseRefineLvl);
void BuildLvlTemplate(int domSize[2], int baseRefineLvl, int nLevelsInput,
	cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);
void CalculateNumElements(int domSize[2],int *nCellCurr,int *nEdgeCurr,int *nVertCurr);
void MemCopyCellStruct(cellTemplate *original, cellTemplate *destination, int nElm);
void Allocatecelldatstruct();
void InitialiseGridFromFile();
void DataIn();
void DataArrayToStruct();
void OutputTemplateGrid(int domSize[2], int lvlGrid);
void GenerateTemplateGrid(int lvlGenerate);
void RefineGrid();
void DeAllocate();
void DeAllocateTemplate(int domSize[2], cellTemplate *cellAct, edgeTemplate *edgeAct, vertexTemplate *vertAct);
void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList);
void GridInitialisation();
void OutputGrid(int jj);
int SortOrder(const void * elem1, const void * elem2);
void IdentifyRefineCell(int refinLvl,int **posGridRef,int **indGridRef, int *nRefineCell);
void IdentifyRefineEdge(int *posCellRefine, int *indCellRefine,int nCellRefine,
		int **posEdgeRefine, int **indEdgeRefine,int *nEdgeRefine, edgeTemplate *edgestructAct, int nEdge);
void RemoveIdenticalEntries_int(int **array, int nArray, int *nNewArray);
int imin(int a, int b);
int imax(int a, int b);
double f_min(double a, double b);
double f_max(double a, double b);
int min_array(int *array, int nArray);
int pos_min_array(int *array, int nArray);
int max_array(int *array, int nArray);
int pos_max_array(int *array, int nArray);
double fmin_array(double *array, int nArray);
int pos_fmin_array(double *array, int nArray);
double fmax_array(double *array, int nArray);
int pos_fmax_array(double *array, int nArray);
void RefineSelectedEdges(int domSize[2],int *posEdgeRefine,int *indEdgeRefine, int nEdgeRefine);
void* ArrayFromStruct(void* ptr, int nStruct, int sizMember, int sizType);
void ExtractEdgeChain(int *posEdge, int *indEdge,int nEdge ,int *ordPosEdge,
	int *ordIndEdge, int *ordIndVert,edgeTemplate *edgestructAct);
int LowerLeftVertex(double *coordList, int nEdge);
double* ExtractCoord(int nEdge, int *ordIndVert,int * ordPosVert, vertexTemplate *vertstructAct, int nVertAct);
double* ExtractAngles(double *vecTest, int nVec, int nDim);
int IsClockWiseChain(int *minVertPosRet,int nEdge, int *ordIndVert, int *ordPosVert, vertexTemplate *vertstructAct, int nVertAct);
void ReorderArray(void *oldArray, int *newOrder ,int nArray, int sizArray);
void OrderEdgeChain(int nEdge, int *ordPosEdge, int *ordIndEdge, int *ordIndVert,int *ordPosVert, vertexTemplate *vertstructAct, int nVertAct);
void GenerateIndMatch(int domSize[2], int posCellRefine ,int *posEdgeSideTemp,int *posVertSideTemp,int *posCellAddTemp,
	int *posEdgeAddTemp, int *posVertAddTemp, int *ordPosEdge, int *ordPosVert, int nCurrEdge,
	int **convCellList, int **convEdgeList, int **convVertList);
void ModifBorderEdges(int nCurrEdge, int posCellRefine, int *ordPosEdge, int *posEdgeSideTemp, int *convCellList);
void ExtendGridStructures(int nCurrEdge,int nCellTemplate,int nEdgeTemplate,int nVertTemplate);
void CopyEdgeStruct(int nCurrEdge,int domSize[2],int *posEdgeAddTemp,
	int *convEdgeList,int *convCellList,int *convVertList);
void CopyVertStruct(int nCurrEdge,int domSize[2],int *posVertAddTemp,
	int *convVertList, int *ordPosVert);
void CopyCellStruct(int domSize[2],int posCellRefine,int *convCellList);
void MergeNewCell(int domSize[2], int posCellRefine ,int *posEdgeSideTemp,int *posVertSideTemp,int *posCellAddTemp,
	int *posEdgeAddTemp, int *posVertAddTemp, int *ordPosEdge, int *ordPosVert, int nCurrEdge);
void RefineCell(int domSize[2],int posCellRefine,int indCellRefine,int *posEdgeSideTemp,
	int *posVertSideTemp,int *posCellAddTemp,int *posEdgeAddTemp, int *posVertAddTemp);
void PrepareTemplateInfo(int domSize[2],int **posEdgeSide,int **posVertSide,int **posCellAdd,
	int **posEdgeAdd,int **posVertAdd);
void RefineSelectedCells(int domSize[2],int *posCellRefine,int *indCellRefine,int nCellRefine);
int CompareDouble(double elem1, double elem2);
void AllocateSolutionStruct(int nCellCurr, int nEdgeCurr, int nVertCurr, int nLevelsAct,
		cellTemplate **cellstructTempOut,vertexTemplate **vertstructTempOut,edgeTemplate **edgestructTempOut);

//void OutputTemplateGrid(int domSize[2], int lvlGrid);
//void GenerateTemplateGrids();

#endif

