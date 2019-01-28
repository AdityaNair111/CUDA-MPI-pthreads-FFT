#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "complex.h"
#include "input_image.h"


#include <mpi.h>

const float PI = 3.14159265358979f;
using namespace std;

void separate (Complex *array,int n)
{
	Complex temp[n/2];
	for (int i = 0; i < n/2; i++)
	{
		temp[i]=array[i*2+1];
	}
	for (int i = 0; i < n/2; i++)
	{
		array[i]= array[i*2];
	}
	for (int i = 0; i < n/2; i++)
	{
		array[i+n/2] = temp[i];
	}
}
void fft (Complex *array,int N)
{
	if (N<2)
	{
		//Very happy/
	}
	else
	{
		separate(array,N);
		fft(array,N/2);
		fft(array+N/2,N/2);
		for (int k = 0; k < N/2; k++)
		{
			Complex even = array[k];
			Complex odd = array[k+N/2];
			Complex w = Complex(cos((float)2*PI*k/N),sin((float)-2*PI*k/N));
			array[k]=even + w*odd;
			array[k+N/2]=even - w*odd;
		}
	}
}

int main(int argc, char *argv[])
{
	//cout<<"in here"<<endl;
	InputImage im = InputImage(argv[2]);
	int Dim = im.get_height();


	//for (int i = 0; i < Dim; i++)
	//{
		//cout<<indexx[i]<<" ";
	//}
	//cout<<endl;
	string direction_string;
	direction_string = argv[1];
	//cout<<direction_string<<endl;
	string f ("forward");
	int direction = 0;
	for (int i = 0; i < direction_string.length(); i++)
	{
		if (int(direction_string.at(i)) != int(f.at(i)))
		{
			direction = 1;
		}
	}
	Complex *Input = im.get_image_data();
	Complex twiddle=Complex(0,0);
	float theta;
	Complex sum;

	MPI_Init(&argc,&argv);
	int numProcs,proc_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	Complex current_grid[Dim*Dim/numProcs];
	Complex backUp_grid[Dim*Dim/numProcs];

	for (int i = 0; i < Dim/numProcs ; i++)
	{
		for (int j = 0; j < Dim; j++)
		{
			current_grid[i*Dim+j]=Input[(proc_rank*Dim/numProcs + i)*Dim + j];
		}
	}
	if (direction==0)
	{
		for (int m = 0; m < Dim/numProcs; m++)
		{
			fft(m*Dim + current_grid, Dim);
		}
	}
	else
	{
		for (int m = 0; m < Dim/numProcs; m++)
		{
			for (int j = 0; j < Dim; j++)
			{
				sum=Complex(0,0);
				theta=2*PI/Dim*j;
				
				for (int i = 0; i < Dim; i++)
				{
					twiddle=Complex(cos(theta*i),sin(theta*i));
					sum=sum+twiddle*Input[(proc_rank*Dim/numProcs + m)*Dim + i];
				}
				current_grid[m*Dim+j]=sum*Complex((float)1/Dim,0);
			}
			
		}
	}
	

		MPI_Barrier(MPI_COMM_WORLD);

		if (proc_rank==0)
		{
			Complex Trans[Dim*Dim];
			
				for (int i = 0; i < Dim*Dim/numProcs; i++)
				{
					Input[i]=current_grid[i];
				}
				
			for (int i = 1; i < numProcs; i++)
			{
				MPI_Status status;
				MPI_Recv(&Input[i*Dim*Dim/numProcs],Dim*Dim/numProcs,MPI_COMPLEX,i,1,MPI_COMM_WORLD,&status);
			}

			int i,j;
			for (int n = 0; n < Dim*Dim; n++) //transpose
			{
				i=n/Dim;j=n%Dim;
				Trans[n]=Input[j*Dim+i];
			}
				
			
				for (int i = 0; i < Dim*Dim/numProcs; i++)
				{
					backUp_grid[i]=Trans[i];
				}
				for (int i = 1; i < numProcs; i++)
				{
					MPI_Send(&Trans[i*Dim*Dim/numProcs],Dim*Dim/numProcs,MPI_COMPLEX,i,1,MPI_COMM_WORLD);
				}
			
	
				/*if (direction==0)
			    {
			    	im.save_image_data(argv[3], &TempTrans[0], Dim, Dim);
			    }
			    else
			    {
			    	im.save_image_data_real(argv[3], &TempTrans[0], Dim, Dim);
			    }
				for (int j = 0; j < Dim; j++)
				{
					for (int i = 0; i < Dim; i++)
					{
						cout<<TempTrans[j*Dim + i]<< " ,";
					}
					cout<<endl;
				}*/

		}
		else
		{
		
			MPI_Send(&current_grid[0],Dim*Dim/numProcs,MPI_COMPLEX,0,1,MPI_COMM_WORLD);
		
			MPI_Status status;
			MPI_Recv(&backUp_grid[0],Dim*Dim/numProcs,MPI_COMPLEX,0,1,MPI_COMM_WORLD,&status);
		}
	MPI_Barrier(MPI_COMM_WORLD);



	if(direction==0)
	{
		for (int m = 0; m < Dim/numProcs; m++)
		{
			fft(m*Dim + backUp_grid, Dim);
		}
	}
	else
	{
		for (int m = 0; m < Dim/numProcs; m++)
			{
				for (int j = 0; j < Dim; j++)
				{
					sum = Complex(0,0);
					theta=2*PI/Dim*j;
					
					for (int i = 0; i < Dim; i++)
					{
						twiddle=Complex(cos(theta*i),sin(theta*i));
						sum=sum+twiddle*backUp_grid[m*Dim + i];
					}
					current_grid[m*Dim+j]=sum*Complex((float)1/Dim,0);
				}
				
			}
	}

		MPI_Barrier(MPI_COMM_WORLD);

		if (proc_rank==0)
		{
			Complex Trans[Dim*Dim];
			if (direction==0)
			{
				for (int i = 0; i < Dim*Dim/numProcs; i++)
				{
					Input[i]=backUp_grid[i];
				}
			}
			else
			{
				for (int i = 0; i < Dim*Dim/numProcs; i++)
				{
					Input[i]=current_grid[i];
				}
			}
				
				
			for (int i = 1; i < numProcs; i++)
			{
				MPI_Status status;
				MPI_Recv(&Input[i*Dim*Dim/numProcs],Dim*Dim/numProcs,MPI_COMPLEX,i,1,MPI_COMM_WORLD,&status);
			}

			int i,j;
			for (int n = 0; n < Dim*Dim; n++) //transpose
			{
				i=n/Dim;j=n%Dim;
				Trans[n]=Input[j*Dim+i];
			}
			
	
				if (direction==0)
			    {
			    	im.save_image_data(argv[3], &Trans[0], Dim, Dim);
			    }
			    else
			    {
			    	im.save_image_data_real(argv[3], &Trans[0], Dim, Dim);
			    }
				/*for (int j = 0; j < Dim; j++)
				{
					for (int i = 0; i < Dim; i++)
					{
						cout<<Trans[j*Dim + i]<< " ,";
					}
					cout<<endl;
				}*/

		}
		else
		{
			if (direction==0)
			{
				MPI_Send(&backUp_grid[0],Dim*Dim/numProcs,MPI_COMPLEX,0,1,MPI_COMM_WORLD);
			}
			else
			{
				MPI_Send(&current_grid[0],Dim*Dim/numProcs,MPI_COMPLEX,0,1,MPI_COMM_WORLD);
			}
			
		}
	


	MPI_Finalize();
	return EXIT_SUCCESS;
}