#include <stdio.h>
#include <stdlib.h>

void getdata (FILE *file, int *rows, int *cols, double ***matrix, double **vector)
{
	file = fopen ("gauss3.dat", "rw");	
	if (file == NULL)
	{
		perror ("fopen");
		exit (EXIT_FAILURE);
	}

	fscanf (file, "%d", rows);
	fscanf (file, "%d", cols);		
	printf ("%d\n", *rows);
	printf ("%d\n", *cols);

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

int max_row (double **matrix, int diag, int rows)
{
	int tmp;
	tmp = diag;
	double max, abs_val;
	max = matrix[tmp][tmp] > 0 ? matrix[tmp][tmp] : -matrix[tmp][tmp];
	int i;
	for (i = 1; i < rows - diag; ++i)
	{
		abs_val = matrix[diag+i][diag] > 0 ? matrix[diag+i][diag] : -matrix[diag+i][diag];
		if (abs_val > max)
		{
			tmp = diag + i;
			max = abs_val;
			
		}
	}
	printf ("max_row %d : %f\n", tmp, max);
	return tmp;
}
void interchange_row (double ***matrix, int diag, int maxrow, int rows)
{
	int j;
	double tmp;
	for (j = 0; j < rows; ++j)
	{
		tmp = (*matrix)[diag][j];
		(*matrix)[diag][j] = (*matrix)[maxrow][j];
		(*matrix)[maxrow][j] = tmp;
	}
}

void print_matrix (double ***matrix, int rows, int cols)
{
	int i,j;
	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < cols; j++)
		{
			printf ("%f ", (*matrix)[i][j]);
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
		maxrow = max_row (*matrix, i, rows);

		// row interchange
		if (maxrow != i)
		{
			interchange_row (matrix, i, maxrow, rows);
			print_matrix (matrix, rows, cols);
			printf ("\n");
		}

		/*
		// operate rows below diagonal
		for (int j = i + 1; j < dim; ++j)
		{
			zerooperate_row (matrix, i, j, dim);
		}

		print_matrix (matrix, dim);
		printf ("\n");
		*/

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
	printf ("rows = %d, colum = %d \n", rows, cols);
	privateprivot_gauss (&A, &b, &x, rows, cols);
	return 0;
}
