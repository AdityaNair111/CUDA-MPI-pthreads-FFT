#include<iostream>
#include<fstream>
#include<string>
#include<algorithm>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits>
#include<iomanip>
#include <chrono>
#include "input_image.h"
#include "complex.h"

#define BLOCK_SIZE 16
using namespace std;
const float PI = 3.14159265358979f;

class MyComplex{
	public:
		__device__ __host__ MyComplex() : real(0.0f), imag(0.0f) {} 
		__device__ __host__ MyComplex(float r, float i) : real(r), imag(i) {}
		__device__ __host__ MyComplex operator+(const MyComplex &b) const 
		{
		MyComplex a;
		a.real = real + b.real;
		a.imag = imag + b.imag;
		return a;
		}
		__device__ __host__ MyComplex operator*(const MyComplex &b) const 
		{
		MyComplex a;
		a.real=real*b.real-imag*b.imag;
		a.imag=imag*b.real+real*b.imag;
		return a;
		}
		float real;
		float imag;
	};
std::ostream& operator<< (std::ostream& os, const MyComplex& rhs) 
{
	MyComplex c(rhs);
	if(fabsf(rhs.imag) < 1e-10) c.imag = 0.0f;
	if(fabsf(rhs.real) < 1e-10) c.real = 0.0f;

	if(c.imag == 0) os << c.real;
	else os << "(" << c.real << "," << c.imag << ")";
	return os;
}

__global__ void matrix1DFFT(MyComplex *MatOld,MyComplex *MatNew,MyComplex *Temp,int Rows,int Cols,int Operation,int direction)
{
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	if (i<Cols && j <Rows)
	{
		if (direction==0)  //Forward
		{
			if (Operation==0) // operating on rows
			{
				MyComplex sum=MyComplex(0, 0);
				float theta=2*PI/Cols*j;
				MyComplex twiddle=MyComplex(0,0);			
				for (int jj = 0; jj < Cols; jj++)
				{
					twiddle=MyComplex(cos(theta*jj),sin(-theta*jj));
					sum=sum+MatOld[i*Cols+jj]*twiddle;
			 	}
			 	Temp[i*Cols+j]=sum;
			}
			else // operating on columns
			{
				MyComplex sum=MyComplex(0, 0);
				float theta=2*PI/Rows*i;
				MyComplex twiddle=MyComplex(0,0);			
				for (int ii = 0; ii < Rows; ii++)
				{
					twiddle=MyComplex(cos(theta*ii),sin(-theta*ii));
					sum=sum+Temp[ii*Cols+j]*twiddle;
			 	}
			 	MatNew[i*Cols+j]=sum;
			}
		}
		else // Reverse
		{
			if (Operation==0) // operating on rows
			{
				MyComplex sum=MyComplex(0, 0);
				float theta=2*PI/Cols*j;
				MyComplex twiddle=MyComplex(0,0);			
				for (int jj = 0; jj < Cols; jj++)
				{
					twiddle=MyComplex(cos(theta*jj),sin(theta*jj));
					sum=sum+MatOld[i*Cols+jj]*twiddle;
			 	}
			 	Temp[i*Cols+j]=sum*MyComplex(((float)1/Cols),0);
			}
			else // operating on columns
			{
				MyComplex sum=MyComplex(0, 0);
				float theta=2*PI/Rows*i;
				MyComplex twiddle=MyComplex(0,0);			
				for (int ii = 0; ii < Rows; ii++)
				{
					twiddle=MyComplex(cos(theta*ii),sin(theta*ii));
					sum=sum+Temp[ii*Cols+j]*twiddle;
			 	}
			 	MatNew[i*Cols+j]=sum*MyComplex(((float)1/Cols),0);
			}
		}	
	}
}
int main(int argc, char const *argv[])
{

	InputImage im=InputImage(argv[2]);
	int Rows=im.get_height();
	int Cols=im.get_width();

	string direction_string;
	direction_string=argv[1];
	string f ("forward");
	int direction =0;
	for (int i = 0; i < direction_string.length(); i++)
	{
		if (int(direction_string.at(i))!=int(f.at(i)))
		{
			direction=1;
		}
	}
	cout<<"Direction : "<<direction<<endl;
	
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	MyComplex *d_MatOld,*d_MatNew,*d_Temp;
	int matSize= Rows*Cols*2*sizeof(float);
	cudaMalloc((void **)&d_MatOld,matSize);
	cudaMalloc((void **)&d_MatNew,matSize);
	cudaMalloc((void **)&d_Temp,matSize);
	cudaMemcpy(d_MatOld,im.get_image_data(),matSize,cudaMemcpyHostToDevice);

	dim3 Block(BLOCK_SIZE,BLOCK_SIZE);
	dim3 Grid(ceil((Rows+BLOCK_SIZE-1)/BLOCK_SIZE),ceil((Cols+BLOCK_SIZE-1)/BLOCK_SIZE));
	if (direction==0)
	{
		matrix1DFFT<<<Grid,Block>>>(d_MatOld,d_MatNew,d_Temp,Rows,Cols,0,direction);
		matrix1DFFT<<<Grid,Block>>>(d_MatOld,d_MatNew,d_Temp,Rows,Cols,1,direction);
	}
	else
	{
		matrix1DFFT<<<Grid,Block>>>(d_MatOld,d_MatNew,d_Temp,Rows,Cols,0,direction);
		matrix1DFFT<<<Grid,Block>>>(d_MatOld,d_MatNew,d_Temp,Rows,Cols,1,direction);
	}
	
	cudaMemcpy(im.get_image_data(),d_MatNew,matSize,cudaMemcpyDeviceToHost);
	cudaFree(d_Temp);cudaFree(d_MatNew);cudaFree(d_MatOld);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start; 
    std::time_t end_time = std::chrono::system_clock::to_time_t(end); 
    std::cout << "finished computation at " << std::ctime(&end_time) 
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    if (direction==0)
    {
    	im.save_image_data(argv[3], im.get_image_data(), Cols, Rows);
    }
    else
    {
    	im.save_image_data_real(argv[3], im.get_image_data(), Cols, Rows);
    }
	return 0;
}
