#include <stdio.h>
#include <stdlib.h>

#define ABS(number) (number > 0 ? number : -number)

void getdata (FILE *file, int *rows, int *cols, double ***matrix, double **vector)
{
	file = fopen ("gauss4.dat", "rw");	
	if (file == NULL)
	{
		perror ("fopen");
		exit (EXIT_FAILURE);
	}

	fscanf (file, "%d", rows);
	fscanf (file, "%d", cols);		

	*matrix = (double **)malloc ((*rows)*sizeof(double*));
	*vector = (double *)malloc ((*rows)*sizeof(double));	
	int i = 0;
	for (i = 0; i < *rows; ++i)
	{
		(*matrix)[i] = (double*) malloc ((*cols) * sizeof (double));
	}

	int j = 0;
	for (i = 0; i < *rows; ++i)
	{
		for (j = 0; j < *cols; ++j)
		{
			fscanf (file, "%lf", &(*matrix)[i][j]);
			//printf ("%lf  ", matrix[i][j]);
		}
		//printf ("\n");
	}

	for (i = 0; i < *rows; ++i)
	{
		fscanf (file, "%lf", &(*vector)[i]);
		//printf ("%lf\n", vector[i]);
	}
	
}

int max_row (double ***matrix, int diag, int rows)
{
	int tmp;
	tmp = diag;
	double max, abs_val;
	max = ABS((*matrix)[tmp][tmp]);
	int i;
	for (i = 0; i < rows - diag; ++i)
	{
		if (ABS ((*matrix)[diag+i][diag]) > max)
		{
			tmp = diag + i;
			max = ABS ((*matrix)[diag+i][diag]);
		}
	}
	return tmp;
}

void print_matrix (double ***matrix, int rows, int cols)
{
	int i,j;
	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < cols; j++)
		{
			printf ("%lf ", (*matrix)[i][j]);
		}
		printf ("\n");
	}
}

void print_vector (double **vec, int dim)
{
	for (int i = 0; i < dim; ++i)
	{
		printf ("%lf ", (*vec)[i]);
	}
	printf ("\n");
}

void interchange_row (double ***matrix, int diag, int maxrow, int cols)
{
	int j;
	double tmp;
	for (j = 0; j < cols; ++j)
	{
		tmp = (*matrix)[diag][j];
		(*matrix)[diag][j] = (*matrix)[maxrow][j];
		(*matrix)[maxrow][j] = tmp;
	}
}
void zerooperate_row (double ***matrix, int diag, int row, int cols)
{
	int i;
	double tmp;
	tmp = (*matrix)[row][diag]/(*matrix)[diag][diag];
	printf ("factor : %lf/%lf\n", (*matrix)[row][0], (*matrix)[diag][diag]);
	for (i = 0; i < cols; ++i)
	{
		(*matrix)[row][i] -= tmp*(*matrix)[diag][i];
	}
}

void solve (double ***matrix, int rows, int cols, double **solution)
{
	*solution = (double *) malloc ((cols-1) * sizeof(double));
	(*solution)[cols - 2] = 0;
	double tmp;
	for (int i = rows - 1; i > -1; --i)
	{
		tmp = (*matrix)[i][cols - 1];
		for (int j = cols - 2; j > i; --j)
		{
			 tmp -= (*matrix)[i][j]*(*solution)[j];
		}
		(*solution)[i] = tmp / (*matrix)[i][i];
	}
	printf ("\n");
}

void privateprivot_gauss (double ***matrix, double **vector, double **solution, int rows, int cols)
{
	//find max row 
	//row interchange
	//make 0 for all row behin diagonal
	int maxrow;
	int num_iter = rows < cols ? rows : cols;
	for (int i = 0; i < num_iter; ++i)
	{
		//find max row
		maxrow = max_row (matrix, i, rows);
		printf ("\nmax_row : #%d\n", maxrow + 1);

		// row interchange
		if (maxrow != i)
		{
			printf ("interchange line #%d and #%d\n", maxrow + 1, i+ 1);
			interchange_row (matrix, i, maxrow, cols);
			print_matrix (matrix, rows, cols);
			printf ("\n");
		}

		// operate rows below diagonal
		for (int j = i + 1; j < rows; ++j)
		{
			zerooperate_row (matrix, i, j, cols);
			print_matrix (matrix, rows, cols);
		}

		// backward substitution
		solve (matrix, rows, cols, solution);

	}
	
}

int main (void)
{
	FILE *stream;
	int rows, cols;
	double **A;
	double *b;
	//double *x = (double *) malloc (dim * sizeof (double));
	double *x;
	getdata (stream, &rows, &cols, &A, &b);
	print_matrix (&A, rows, cols);
	privateprivot_gauss (&A, &b, &x, rows, cols);
	printf ("\nsolution :\n");
	print_vector (&x, cols-1);
	free (A);
	free (b);
	free (x);
	return 0;
}
