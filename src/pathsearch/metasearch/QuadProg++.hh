#ifndef _QUADPROGPP
#define _QUADPROGPP
#ifndef MATRIX_DIM
#define MATRIX_DIM 500
#endif

double solve_quadprog(double G[][MATRIX_DIM], double g0[], int n, 
                      double CE[][MATRIX_DIM], double ce0[], int p, 
                      double CI[][MATRIX_DIM], double ci0[], int m,
                      double x[]);
void print_matrix(const char* name, double A[][MATRIX_DIM], int n);//YU const
void print_vector(char* name, double v[], int n);//YU
#endif // #define _QUADPROGPP
