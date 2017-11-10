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


void ReadInGriduns();
void DataIn();
int* FindNegVertex(int *nVertNeg);
int *FindDelEdge(int *nEdgeDel,int *posVertNeg,int nVertNeg);
void FindObjNum(int *lookup, int *listInd, int *dest, int nLook, int nList);
void RemoveIdenticalEntries_int(int **array, int nArray, int *nNewArray);
int SortOrder(const void * elem1, const void * elem2) ;
void BuildCellMatchList();
void BuildVertMatchList();
void WriteOutGrid();
void WriteOutSymPlane();
void FreeAllMeshSym();

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
#endif