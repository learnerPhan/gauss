#include <stdio.h>
#include <stdlib.h>

void getdata (FILE *file, int *dim, double **matrix, double *vector)
{
	file = fopen ("gauss2.dat", "rw");	
	if (file == NULL)
	{
		perror ("fopen");
		exit (EXIT_FAILURE);
	}

	fscanf (file, "%d", dim);		
	printf ("%d\n", *dim);
	
	int i = 0;
	for (i = 0; i < *dim; ++i)
	{
		matrix[i] = (double*) malloc ((*dim) * sizeof (double));
	}

	int j = 0;
	for (i = 0; i < *dim; ++i)
	{
		for (j = 0; j < *dim; ++j)
		{
			fscanf (file, "%lf", &matrix[i][j]);
			printf ("%lf  ", matrix[i][j]);
		}
		printf ("\n");
	}

	for (i = 0; i < *dim; ++i)
	{
		fscanf (file, "%lf", &vector[i]);
		printf ("%lf\n", vector[i]);
	}
	
}

int max_row (double **matrix, int diag, int dim)
{
	int maxrow = diag;
	double max = matrix[maxrow][maxrow];
	int i;
	for (i = 1; i < dim - diag; ++i)
	{
		if (matrix[diag+i][diag] > max)
		{
			maxrow = diag + i;
			max = matrix[maxrow][diag];
			
		}
	}
	printf ("maxrow %d : %lf\n", maxrow, max);
	return maxrow;
}
void interchange_row (double **matrix, int diag, int maxrow, int dim)
{
	int j;
	double tmp;
	for (j = 0; j < dim; ++j)
	{
		tmp = matrix[diag][j];
		matrix[diag][j] = matrix[maxrow][j];
		matrix[maxrow][j] = tmp;
	}
}

void print_matrix (double **matrix, int dim)
{
	int i,j;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; j++)
		{
			printf ("%lf ", matrix[i][j]);
		}
		printf ("\n");
	}
}

void zerooperate_row (double **matrix, int diag, int row, int dim)
{
	int i;
	double tmp;
	tmp = matrix[row][0]/matrix[diag][diag];
	for (i = 0; i < dim; ++i)
	{
		matrix[row][i] -= tmp*matrix[row][i];
	}
}

void privateprivot_gauss (double **matrix, double *vector, double *solution, int dim)
{
	//find max row 
	//row interchange
	//make 0 for all row behin diagonal
	int maxrow;
	for (int i = 0; i < dim; ++i)
	{
		//find max row
		maxrow = max_row (matrix, i, dim);

		// row interchange
		if (maxrow != i)
		{
			interchange_row (matrix, i, maxrow, dim);
		}

		print_matrix (matrix, dim);
		printf ("\n");

		// operate rows below diagonal
		for (int j = i + 1; j < dim; ++j)
		{
			zerooperate_row (matrix, i, j, dim);
		}

		print_matrix (matrix, dim);
		printf ("\n");

	}
	
}

int main (void)
{
	FILE *stream;
	int dim;
	// allocate matrix A
	double **A = (double **) malloc (dim * sizeof (double*));
	//allocate vector b
	double *b = (double *) malloc (dim * sizeof (double));
	//allocate vector x
	double *x = (double *) malloc (dim * sizeof (double));
	getdata (stream, &dim, A, b);
	privateprivot_gauss (A, b, x, dim);
	return 0;
}
