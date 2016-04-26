#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cuda_runtime.h>

__device__ int position;			//index of the largest value
__device__ int largest;				//value of the largest value
int lenString=594;
int maxNumStrings = 1000000;                           
int threshold = 2;

__global__ void kernel(int set) { position = set;}

__global__ void search(int *d_b, int *d_c, int max_count)
{
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	if(d_c[my_id]==0 && (d_b[my_id] == largest )&&(my_id < max_count))
	{
		position = my_id;
	}
}


__global__ void populate (int *d_b, int *copy_db, int *d_c, int size) {
	int n = 0;
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if (my_id < size) {
		n = abs((bool)d_c[my_id] - 1);
		copy_db[my_id] = d_b[my_id] * n;
	}		
}
__device__ void cuda_select(int *db, int size)
{
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	if(my_id < size)
	{
		if(db[2*my_id] > db[2*my_id + 1])
			db[my_id] = db[2*my_id];
		else
			db[my_id] = db[2*my_id + 1];
	}
	
}
__global__ void select(int *db, int size)
{
	int height = (int)ceil(log2((double)size));
	int i = 0;
	
	for(i = 0; i < height; i++)
	{
		size = (int)ceil((double) size/2);
		//int threads_num = 512, blocks_num;
		//blocks_num = (int)ceil((float)size/threads_num);
		cuda_select(db, size);
	}
	largest = db[0];
}

__global__ void Compare(char *d_a, int *d_b, int *d_c, int max_count, int lenString, int threshold){

	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	
		
	if ((my_id < max_count) && (d_c[my_id] == 0) && (my_id != position)){
		
		int x, i, diffs = 0, stop =0;

		for (x=0;x<lenString;x++){
			diffs += (bool)(d_a[(lenString*position)+x]^d_a[(my_id*lenString)+x]);
			
			if (diffs > threshold){
				break;}
		}
		
		if (diffs <= threshold){
			d_b[position] += d_b[my_id];
			d_c[position] = 2;
			d_c[my_id] = 1;
			
		}
		else {
			d_c[position] = 2;
		}
	}

}


int main(int argc, char** argv) {//allocation of variables
	char *strings, *d_a;
	int *counts, *merged, *d_b, *d_c;	//host copy of a
	int *largest, *copy_db, *copy_copy_db;
	char copy[lenString+1]; //string to copy in info
	int numbers=0;
	int i=0, actual_count=0;
	int size_string = maxNumStrings*sizeof(char)*(lenString+1);
	int size_int = maxNumStrings*sizeof(int);
	int size_int_2 = sizeof(int);
	struct timeval start, end; 				//using time
	double wallTime;
	cudaError_t status = (cudaError_t)0;

	int set_value;


	//opening the file
	FILE *fp;
	fp=fopen("/cluster/home/charliep/courses/cs360/single-linkage-clustering/Iceland2014.trim.contigs.good.unique.good.filter.unique.count.fasta", "r");


	if (!(strings= (char *)malloc(size_string))) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);}
	if (!(counts= (int*)malloc(size_int))) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);}
	merged = (int *)malloc(size_int);
	copy_copy_db = (int *)malloc(size_int);
	cudaMemset(&position,0,sizeof(int));
	cudaMemset(&largest,0,sizeof(int));

	while( fscanf(fp,"%s %d", copy, &numbers) != EOF && actual_count <100){
		strcpy(&strings[i],copy);
		counts[actual_count]=numbers;
		//printf("%s\n", copy);
		//printf("%s\n", &a[i]);
		i=i+lenString;
		actual_count++;
		
		}
	fclose(fp);
	
	cudaMalloc(&d_a, size_string);
	cudaMalloc(&d_b, size_int);
	cudaMalloc(&d_c, size_int);
	cudaMalloc(&copy_db, size_int);
	
	for(int i =0; i<actual_count;i++)
	{
		printf("%d , %d\n", counts[i], merged[i]);
	}
	
	
	cudaMemcpy(d_a, strings, size_string, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, counts, size_int, cudaMemcpyHostToDevice);
	cudaMemcpy(d_c, merged, size_int, cudaMemcpyHostToDevice);
	
	int threads_num = 512, blocks_num;
	blocks_num = (int)ceil((float)actual_count/threads_num);


	
	populate<<<blocks_num, threads_num>>>(d_b, copy_db, d_c, actual_count); //this does what it is suppose to do so far checking next part
	select<<<blocks_num, threads_num>>>(copy_db, actual_count);							//not selecting the largest value
	
	kernel<<<1,1>>>(set_value);
	search<<<blocks_num, threads_num>>>(d_b, d_c, actual_count);
	Compare<<<blocks_num, threads_num>>>(d_a, d_b, d_c, actual_count, lenString, threshold);

		
	cudaMemcpy(counts, d_b, size_int, cudaMemcpyDeviceToHost);
	cudaMemcpy(merged, d_c, size_int, cudaMemcpyDeviceToHost);
	
	printf("\n");
	
	
	for(int i =0; i<actual_count;i++)
	{
		printf("%d , %d\n", counts[i], merged[i]);
	}
	
		
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	free(strings);
	free(counts);
	free(merged);
}
