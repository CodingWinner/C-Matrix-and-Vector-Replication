#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef int Rows;
typedef int Columns;

typedef struct
{
    double **values;
    Rows r;
    Columns c;
} Matrix;

typedef struct
{
    double *values;
    Columns c;
} Vector;