#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<conio.h>
//#include "winbgi2.h"

void UpperTriangleMatrix(double **M1,double **UM, int n);   //Conversion to Upper Traingular Matrix
void determinant(double **umatrix, int n,double *det); // Determinant of the matrix
void main()
{
	srand((unsigned)time(0));
	double **A,**AT,*K,**C,**UMA,**UMAT,**COM,*x,sum;   //where A- A Matix, AT - transponse of the Matrix A, UMA -upper matrix A
	double M=1.;                                        // UMAT- upper matrix AT, COM- combination of the matrix to solve the system
	double min=-10.0;                                   // x- solution
	double max=10.0;
	int i,j,n;
	printf("\n enter the value of n:");
	scanf("%d",&n);
	K=(double*)malloc(n*sizeof(double));
	x=(double*)malloc(n*sizeof(double));
	UMA=(double**)malloc(n*sizeof(double)); 
	for(i=0;i<n;i++)
	{
		UMA[i]=(double*)malloc(n*sizeof(double));
	}
	UMAT=(double**)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
	{
		UMAT[i]=(double*)malloc(n*sizeof(double));
	}
	C=(double**)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
	{
		C[i]=(double*)malloc(n*sizeof(double));
	}
	A=(double**)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
	{
		A[i]=(double*)malloc(n*sizeof(double));
	}
	AT=(double**)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
	{
		AT[i]=(double*)malloc(n*sizeof(double));
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			A[i][j]= min + rand() / (RAND_MAX / (max - min + 1) + 1);
		}
	}
	// Displaying the Random of matrix A
	 printf("\nRandom Matrix A: \n");
	for(i=0; i<n; ++i)
	{
        for(j=0; j<n; ++j)
        {
            printf("%lf  ", A[i][j]);
            if (j == n-1)
                printf("\n\n");
        }
	}
	for(i=0; i<n; ++i)
	{
		for(j=0; j<n; ++j)
        {
            AT[j][i] = A[i][j];
        }
	}
    // Displaying the transpose of matrix A
    printf("\nTranspose of Matrix A:\n");
    for(i=0; i<n; ++i)
	{
        for(j=0; j<n; ++j)
        {
            printf("%lf  ",AT[i][j]);
            if(j==n-1)
                printf("\n\n");
        }
	}
	UpperTriangleMatrix(AT,UMAT,n);
	double det;
	determinant(UMAT,n,&det);
	printf("\nDeterminant of Matrix AT: \t%lf\n",det);
	 for(i=0; i<n; ++i)
	{
        for(j=0; j<n; ++j)
        {
            M=M*AT[i][j];
            if(j==n-1)
			{
				K[i]=M;
                M=1.;
			}
        }
	 }
	 printf("\n");
	 printf("Evaluating Matrix C {Row wise Multiplication[AT]/det AT}:\n");   // Display Matrix C
	  for(i=0; i<n; ++i)
	  {
        for(j=0; j<1; ++j)
        {
			C[i][j]=K[i]/det;
			printf("%lf",C[i][j]);
            printf("\n\n");
        }
	  }
	  COM=(double**)malloc(n*sizeof(double));  //Combining Matrix A and C to Solve the system
	  for(i=0;i<n;i++)
	  {
		  COM[i]=(double*)malloc(n*sizeof(double));
	  }
	  for(i=0; i<n; i++)
	  {
		  for(j=0; j<=n; j++)
		  {
			  COM[i][j]=A[i][j];
			  if(j==n)
			  {
				  COM[i][j]=C[i][0];
			  }
		  }
	  }
	  printf("\n");
	  printf("Combined Matrix for solving system A*X=C:\n");   
	  for(i=0; i<n; ++i)
	  {
        for(j=0; j<=n; ++j)
        {
			printf("%lf ",COM[i][j]);
			if(j==n)
			{
            printf("\n\n");
			}
         }
	  }
	  double ratio;
	  for(i=0; i<n; i++)
	  {
		  for(j=0; j<n; j++)
		  {
			  if(j>i)
			  {
				  ratio=COM[j][i]/COM[i][i];
				  for(int k=0; k<=n; k++)
				  {
					  COM[j][k]=COM[j][k]-ratio*COM[i][k];
				  }
			  }
		  }
	  }
	  // initializing matrix x to zeros //
  for(int i=0; i<n; i++)
  {
    x[i]=0;
  }
  // Backward substitution loop //
  printf("\n\nBackward Substitution: \n");
  for(i=n-1; i>=0; i--)
  {
    sum=0.;
    for(j=0; j<n; j++)
    {
      if(i!=j)
        sum=sum+COM[i][j]*x[j];
    }
    x[i]=(COM[i][n]-sum)/COM[i][i];
  }
  printf("\nThe solution is: \n");
  for(i=0; i<n; i++)
  {
    printf("\nx%d=%lf\t",i,x[i]);
  }
	  getch();
}
void UpperTriangleMatrix(double **M1,double **UM, int n)   // Conversion to upper traingular method
{
	double ratio;
		for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
        {
            UM[i][j] = M1[i][j];
        }
	}
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
            if(j>i)
			{
                ratio = UM[j][i]/UM[i][i];
                for(int k = 0; k < n; k++)
				{
                    UM[j][k] -= ratio * UM[i][k];
                }
			}
		}
	}
}
void determinant(double **umatrix, int n,double *det)
{
	double d = 1.; 
    for(int i = 0; i < n; i++)
	{
        d *= umatrix[i][i];
	}
	*det=d;
}
