#include "Math.h"

/*
    r is the number of rows
    c is the number of columns
    vals is NULL if the matrix is to be filled with zeroes, otherwise it is a 2d array with r rows and c columns
    these values then replace those of the matrix
 */
Matrix create_matrix(Rows r, Columns c, double **vals)
{
    Matrix m;
    m.r = r;
    m.c = c;
    m.values = malloc(sizeof(double *) * r);
    for (int i = 0; i < r; i++)
    {
        m.values[i] = calloc(c, sizeof(double));
    }
    if (vals != NULL)
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                m.values[i][j] = vals[i][j];
            }
        }
    }
    return m;
}

/*
    c is the size of the vector
    vals is NULL if the vector is to be filled with zeroes otherwise it is a 1d array with c columns
    these values then make up the values of the vector
*/
Vector create_vector(Columns c, double *vals)
{
    Vector v;
    v.c = c;
    v.values = calloc(c, sizeof(double));
    if (vals != NULL)
    {
        for (int i = 0; i < c; i++)
        {
            v.values[i] = vals[i];
        }
    }
}

/*
    Call after done using a matrix
*/
void free_matrix(Matrix m)
{
    for (int i = 0; i < m.r; i++)
    {
        free(m.values[i]);
    }
    free(m.values);
}

/*
    Call after done using a vector
*/
void free_vector(Vector v)
{
    free(v.values);
}

/*
    m1 is the matrix to be copied and m2 is the matrix to get the values of m1
    *Assumes that both matrices have already been created
*/
void copy_matrix(Matrix m1, Matrix *m2)
{
    m2->r = m1.r;
    m2->c = m1.c;
    m2->values = realloc(m2->values, sizeof(double *) * m2->r);

    for (int i = 0; i < m2->r; i++)
    {
        m2->values[i] = calloc(m2->c, sizeof(double));
        for (int j = 0; j < m2->c; j++)
        {
            m2->values[i][j] = m1.values[i][j];
        }
    }
}

/*
    v1 is the vector to be copied and v2 is the target vector of the copying
    *Assumes that both vectors have already been created
*/
void copy_vector(Vector v1, Vector *v2)
{
    v2->c = v1.c;
    v2->values = realloc(v2->values, sizeof(double) * v2->c);

    for (int i = 0; i < v2->c; i++)
    {
        v2->values[i] = v1.values[i];
    }
}

/*
    Flips the shape of the matrix
    ex: {{2, 3, 2}, {4, 3, 4}} turns into {{2, 4}, {3, 3}, {2, 4}}
*/
Matrix transpose_matrix(Matrix m)
{
    Matrix new_m = create_matrix(m.c, m.r, NULL);
    for (int i = 0; i < m.r; i++)
    {
        for (int j = 0; j < m.c; j++)
        {
            new_m.values[j][i] = m.values[i][j];
        }
    }

    copy_matrix(new_m, &m);
    free_matrix(new_m);
}

/*
    Returns a matrix with values equivalent to the sum of the two matrices
*/
Matrix add_matrices(Matrix m1, Matrix m2)
{
    if (m1.r != m2.r || m1.c != m2.c)
    {
        printf("Dimensions do not match, function not executed");
    }
    else
    {
        Matrix m3 = create_matrix(m1.r, m2.c, m1.values);
        for (int i = 0; i < m3.r; i++)
        {
            for (int j = 0; j < m3.c; j++)
            {
                m3.values[i][j] += m2.values[i][j];
            }
        }
        return m3;
    }
}

/*
    Returns a matrix with values equivalent to the difference of the two matrices
*/
Matrix subtract_matrices(Matrix m1, Matrix m2)
{
    if (m1.r != m2.r || m1.c != m2.c)
    {
        printf("Dimensions do not match, function not executed");
    }
    else
    {
        Matrix m3 = create_matrix(m1.r, m2.c, m1.values);
        for (int i = 0; i < m3.r; i++)
        {
            for (int j = 0; j < m3.c; j++)
            {
                m3.values[i][j] -= m2.values[i][j];
            }
        }
        return m3;
    }
}

/*
    Returns sum of each value at index i of vector 1 multiplied by the value at index i of the second vector
*/
double dot_vectors(Vector v1, Vector v2)
{
    if (v1.c != v2.c)
    {
        printf("Error, dimesions do not match");
    }
    else
    {
        double sum = 0;
        for (int i = 0; i < v1.c; i++)
        {
            sum += v1.c * v2.c;
        }
        return sum;
    }
}

/*
    If axis = 0 then converts the row at specified index of matrix into a vector
    else if axis = 1 then converts the column at specified index of matrix into a vector
*/
Vector matrix_part_to_vector(Matrix m, int axis, int index)
{
    Vector v;
    if (axis = 0)
    {
        v = create_vector(m.c, m.values[index]);
    }
    else
    {
        v = create_vector(m.r, m.values[index]);
    }
    return v;
}

/*
    Returns dot product of the two matrices
*/
Matrix dot_matrices(Matrix m1, Matrix m2)
{
    if (m1.r != m2.c || m1.c != m2.r)
    {
        printf("Error, dimensions do not match");
    }
    else
    {
        Vector temp_v;
        Vector temp_v2;
        Matrix m3 = create_matrix(m1.r, m2.c, NULL);
        for (int i = 0; i < m1.r; i++)
        {
            for (int j = 0; j < m2.c; j++)
            {
                temp_v = matrix_part_to_vector(m1, 0, i);
                temp_v2 = matrix_part_to_vector(m2, 1, j);
                m3.values[i][j] = dot_vectors(temp_v, temp_v2);
                free_vector(temp_v);
                free_vector(temp_v2);
            }
        }
        return m3;
    }
}

/*
    Returns each value in m1 multiplied by the corresponding value at the same index in m2
*/
Matrix cross_matrices(Matrix m1, Matrix m2)
{
    if (m1.r != m2.r || m1.c != m2.c)
    {
        printf("Error, dimensions do not match");
    }
    else
    {
        Matrix m3 = create_matrix(m1.r, m1.c, m1.values);
        for (int i = 0; i < m1.r; i++)
        {
            for (int j = 0; j < m1.c; j++)
            {
                m3.values[i][j] *= m2.values[i][j];
            }
        }
        return m3;
    }
}

/*
    Add number to matrix
*/
Matrix matrix_add_num(Matrix m, double v)
{
    Matrix m3 = create_matrix(m.r, m.c, m.values);
    for (int i = 0; i < m3.r; i++)
    {
        for (int j = 0; j < m3.c; j++)
        {
            m3.values[i][j] += v;
        }
    }
    return m3;
}

/*
    Subtract number from matrix
*/
Matrix matrix_subtract_num(Matrix m, double v)
{
    Matrix m3 = create_matrix(m.r, m.c, m.values);
    for (int i = 0; i < m3.r; i++)
    {
        for (int j = 0; j < m3.c; j++)
        {
            m3.values[i][j] -= v;
        }
    }
    return m3;
}

/*
    Divide matrix by number
*/
Matrix matrix_divide_by_num(Matrix m, double v)
{
    Matrix m3 = create_matrix(m.r, m.c, m.values);
    for (int i = 0; i < m3.r; i++)
    {
        for (int j = 0; j < m3.c; j++)
        {
            m3.values[i][j] /= v;
        }
    }
    return m3;
}

/*
    Multiply matrix by num
*/
Matrix matrix_multiply_by_num(Matrix m, double v)
{
    Matrix m3 = create_matrix(m.r, m.c, m.values);
    for (int i = 0; i < m3.r; i++)
    {
        for (int j = 0; j < m3.c; j++)
        {
            m3.values[i][j] *= v;
        }
    }
    return m3;
}

/*
    Raise matrix to exponent
    Positive num only
*/
Matrix matrix_exponentiate(Matrix m, int v)
{
    Matrix m3 = create_matrix(m.r, m.c, m.values);
    for (int i = 0; i < m3.r; i++)
    {
        for (int j = 0; j < m3.c; j++)
        {
            double temp = m3.values[i][j];
            for (int k = 1; k < v; k++)
            {
                temp *= m3.values[i][j];
            }
            m3.values[i][j] = temp;
        }
    }
    return m3;
}

/*
    Apply a function to the matrix
    Function must have return type double and take in matrix value as sole parameter
*/
Matrix function_matrix(Matrix m, double (*func)(double))
{
    Matrix m2 = create_matrix(m.r, m.c, m.values);
    for (int i = 0; i < m2.r; i++)
    {
        for (int j = 0; j < m2.c; j++)
        {
            m2.values[i][j] = func(m2.values[i][j]);
        }
    }
    return m2;
}

/*
    Fills the matrix with the value
*/
void fill_matrix(Matrix *m, double val)
{
    for (int i = 0; i < m->r; i++)
    {
        for (int j = 0; j < m->c; j++)
        {
            m->values[i][j] = val;
        }
    }
}

/*
    Turn matrix into a vector
*/
Vector flatten_matrix(Matrix m)
{
    Vector v = create_vector(m.r * m.c, NULL);
    for (int i = 0; i < m.r; i++)
    {
        for (int j = 0; j < m.c; j++)
        {
            v.values[(i * m.c) + j] = m.values[i][j];
        }
    }
    return v;
}

/*
    Multiplies each value of v1 in column i by the ith value of column 2
*/
Vector cross_vectors(Vector v1, Vector v2)
{
    if (v1.c != v2.c)
    {
        printf("Dimensions do not match\n");
    }
    else
    {
        Vector v3 = create_vector(v1.c, v1.values);
        for (int i = 0; i < v1.c; i++)
        {
            v3.values[i] *= v2.values[i];
        }
        return v3;
    }
}

/*
    Adds each value of v1 in column i by the ith value of column 2
*/
Vector add_vectors(Vector v1, Vector v2)
{
    if (v1.c != v2.c)
    {
        printf("Dimensions do not match\n");
    }
    else
    {
        Vector v3 = create_vector(v1.c, v1.values);
        for (int i = 0; i < v1.c; i++)
        {
            v3.values[i] += v2.values[i];
        }
        return v3;
    }
}

/*
    Subtracts each value of v1 in column i by the ith value of column 2
*/
Vector subtract_vectors(Vector v1, Vector v2)
{
    if (v1.c != v2.c)
    {
        printf("Dimensions do not match\n");
    }
    else
    {
        Vector v3 = create_vector(v1.c, v1.values);
        for (int i = 0; i < v1.c; i++)
        {
            v3.values[i] -= v2.values[i];
        }
        return v3;
    }
}

/*
    Divides each value of v1 in column i by the ith value of column 2
*/
Vector divide_vectors(Vector v1, Vector v2)
{
    if (v1.c != v2.c)
    {
        printf("Dimensions do not match\n");
    }
    else
    {
        Vector v3 = create_vector(v1.c, v1.values);
        for (int i = 0; i < v1.c; i++)
        {
            v3.values[i] /= v2.values[i];
        }
        return v3;
    }
}

/*
    Adds the value specified to each element in the vector
*/
Vector vector_add_num(Vector v, double val)
{
    Vector v2 = create_vector(v.c, v.values);
    for (int i = 0; i < v.c; i++)
    {
        v.values[i] += val;
    }
    return v2;
}

/*
    Subtracts the value specified from each element in the vector
*/
Vector vector_subtract_num(Vector v, double val)
{
    Vector v2 = create_vector(v.c, v.values);
    for (int i = 0; i < v.c; i++)
    {
        v.values[i] -= val;
    }
    return v2;
}

/*
    Multiplies by the value specified to each element in the vector
*/
Vector vector_multiply_num(Vector v, double val)
{
    Vector v2 = create_vector(v.c, v.values);
    for (int i = 0; i < v.c; i++)
    {
        v.values[i] *= val;
    }
    return v2;
}

/*
    Raises each value in the vector v to the valth power
    Positive num only
*/
Vector vector_exponentiate_num(Vector v, int val)
{
    Vector v2 = create_vector(v.c, v.values);
    for (int i = 0; i < v.c; i++)
    {
        double sum = v.values[i];
        for (int k = 1; k < val; k++)
        {
            sum *= v.values[i];
        }
        v.values[i] = sum;
    }
    return v2;
}

/*
    Applies the double function that takes in value double to each element in the vector
*/
Vector function_vector(Vector v, double (*func)(double))
{
    Vector v2 = create_vector(v.c, NULL);
    for (int i = 0; i < v.c; i++)
    {
        v.values[i] = func(v.values[i]);
    }
    return v;
    2;
}

/*
    Returns a matrix where a multiplication table is made by lining up the vectors across the x and y axis
    A 0 for axis specifies that v1 goes on top otherwsie v1 goes on the y axis
*/
Matrix vectors_to_matrix(Vector v1, Vector v2, int axis)
{
    if (axis == 0)
    {
        Matrix m3 = create_matrix(v1.c, v2.c, NULL);
        for (int i = 0; i < v1.c; i++)
        {
            for (int j = 0; j < v2.c; j++)
            {
                m3.values[i][j] = v1.values[i] * v2.values[j];
            }
        }
        return m3;
    }
    else
    {
        Matrix m3 = create_matrix(v2.c, v1.c, NULL);
        for (int i = 0; i < v2.c; i++)
        {
            for (int j = 0; j < v1.c; j++)
            {
                m3.values[i][j] = v2.values[i] * v1.values[j];
            }
        }
        return m3;
    }
}

Vector vector_divide_matrix(Matrix m, Vector v);
Vector vector_multiply_matrix(Matrix m, Vector v);
Matrix matrix_divide_vector(Matrix m, Vector v);
Matrix matrix_multiply_vector(Matrix m, Vector v);
Matrix matrix_add_vector(Matrix m, Vector v);
Matrix matrix_subtract_vector(Matrix m, Vector v);

/*
    If axis is zero then it squashes the matrix by summing the rows otherwise squashes the columns up
*/
Vector vectorize_matrix(Matrix m, int axis)
{
    if (axis == 0)
    {
        Vector v = create_vector(m.c, NULL);
        for (int i = 0; i < m.r; i++)
        {
            for (int j = 0; j < m.c; j++)
            {
                v.values[j] += m.values[i][j];
            }
        }
        return v;
    }
    else
    {
        Vector v = create_vector(m.r, NULL);
        for (int i = 0; i < m.c; i++)
        {
            for (int j = 0; j < m.r; j++)
            {
                v.values[j] += m.values[i][j];
            }
        }
        return v;
    }
}
