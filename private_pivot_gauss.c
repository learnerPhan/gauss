/* This is a code for the gaussian elimination
* that is based on lecture https://ocw.mit.edu/courses/chemical-engineering/10-34-numerical-methods-applied-to-chemical-engineering-fall-2005/lecture-notes/lecturenotes123.pdf
*/

#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

/*
* mesure execution time
*/
double gettime()
{
	struct timeval time;
	gettimeofday (&time, NULL);
	return time.tv_sec + (time.tv_sec * 1.0e-6L);
}

#define ABS(number) (number > 0 ? number : -number)

/*
* Read matrix from a file. The last column on the matrix is the vector b
*/
void getdata (FILE **file, int *rows, int *cols, double ***matrix)
{
	*file = fopen ("gauss5.dat", "rw");	
	if (*file == NULL)
	{
		perror ("fopen");
		exit (EXIT_FAILURE);
	}

	fscanf (*file, "%d", rows);
	fscanf (*file, "%d", cols);		

	*matrix = (double **)malloc ((*rows)*sizeof(double*));
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
			fscanf (*file, "%lf", &(*matrix)[i][j]);
		}
	}
	fclose (*file);
}

/*
* find the row that contains the largest magnitude's elements in a column
*/
int max_row (double ***matrix, int diag, int rows)
{
	int tmp;
	tmp = diag;
	double max;
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
	if (max == 0.)
	{
		perror ("Matrix is singular");
		exit (EXIT_FAILURE);
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

/*
* swap the current row with the row found after the max_rox function
*/
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

/*
* row operation in order to have 0 below the current diagonal
*/
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

/*
* having the matrix as the upper triangular, we solve by backward substitution
*/
void back_sub (double ***matrix, int rows, int cols, double **solution)
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

/*
* the private pivot gaussian elimination 
*/
void privateprivot_gauss (double ***matrix, double **solution, int rows, int cols)
{
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
		back_sub (matrix, rows, cols, solution);
	}
}

int main (void)
{
	FILE *stream;
	int rows, cols;
	double **A;
	double *x;
	double begin=0.0, end=0.0;
	getdata (&stream, &rows, &cols, &A);
	print_matrix (&A, rows, cols);
	
	/*start timing*/
	begin = gettime();
	
	/*do gauss*/
	privateprivot_gauss (&A, &x, rows, cols);

	/*end timing*/
	end = gettime();

	printf ("compution time : %g s\n", end - begin);
	printf ("\nsolution :\n");
	print_vector (&x, cols-1);
	for (int i = 0; i < rows; ++i)
	{
		free (A[i]);
	}
	free (A);
	free (x);
	return 0;
}
