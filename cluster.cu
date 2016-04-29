#include <stdio.h>
#include <stdlib.h>
#include <cmath>

__device__ int position;			//index of the largest value
__device__ int largest;				//value of the largest value
int lenString = 593;
int maxNumStrings = 1000000;                           
int threshold = 2;

__global__ void anyLeft(int *d_c, int *remaining, int max_count) {
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	if((d_c[my_id] == 0) && (my_id < max_count)) {
		*remaining = 0;
	}
}

__global__ void search(int *d_b, int *d_c, int max_count) {
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	if((d_c[my_id] == 0) && (d_b[my_id] == largest) && (my_id < max_count)) {
		position = my_id;
	}
}

__global__ void populate (int *d_b, int *copy_db, int *d_c, int size, int *left) {
	int n = 0;
	*left = 1;

	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if (my_id < size) {
		n = abs((bool)d_c[my_id] - 1);
		copy_db[my_id] = d_b[my_id] * n;
	}		
}

__device__ void cuda_select(int *db, int size) {
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if(my_id < size) {
		if(db[2 * my_id] > db[2 * my_id + 1])
			db[my_id] = db[2 * my_id];
		else
			db[my_id] = db[2 * my_id + 1];
	}	
}

__global__ void select(int *db, int size) {
	int height = (int)ceil(log2((double)size));
	int i = 0;
	
	for(i = 0; i < height; i++) {
		size = (int)ceil((double) size/2);
		cuda_select(db, size);
	}
	largest = db[0];
}

__global__ void compare(char *d_a, int *d_b, int *d_c, int max_count, int lenString, int threshold) {

	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if (my_id == position) 
		d_c[my_id] = 2;
	
		
	if ((my_id < max_count) && (d_c[my_id] == 0) && (my_id != position)) {	
		int x, diffs = 0;

		for (x = 0; x < lenString; x++) {
			diffs += (bool)(d_a[(lenString*position)+x]^d_a[(my_id*lenString)+x]);
			
			if (diffs > threshold)
				break;
		}
		
		if (diffs <= threshold) {
			d_b[position] += d_b[my_id];
			d_c[my_id] = 1;
		}
	} 
}

int main(int argc, char** argv) {
	char *strings, *d_a;
	int *counts, *merged, *d_b, *d_c;
	int *largest, *copy_db;
	char copy[lenString+1];
	int numbers = 0;
	int i = 0, actual_count = 0;
	int size_string = maxNumStrings*sizeof(char)*(lenString+1);
	int size_int = maxNumStrings*sizeof(int);
	
	int *any_left, *left;

	//opening the file
	FILE *fp;
	fp = fopen("/cluster/home/charliep/courses/cs360/single-linkage-clustering/Iceland2014.trim.contigs.good.unique.good.filter.unique.count.fasta", "r");

	if (!(strings= (char *)malloc(size_string))) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);
	}

	if (!(counts= (int*)malloc(size_int))) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);
	}

	merged = (int *)malloc(size_int);
	any_left = (int *)malloc(sizeof(int));
	cudaMemset(&position, 0, sizeof(int));
	cudaMemset(&largest, 0, sizeof(int));

	while(fscanf(fp, "%s %d", copy, &numbers) != EOF && actual_count <1000){
		strcpy(&strings[i], copy);
		counts[actual_count] = numbers;
		
		i = i + lenString;
		actual_count++;	
	}
	
	fclose(fp);
	
	cudaMalloc(&d_a, size_string);
	cudaMalloc(&d_b, size_int);
	cudaMalloc(&d_c, size_int);
	cudaMalloc(&copy_db, size_int);
	cudaMalloc(&left, size_int);
	
	cudaMemcpy(d_a, strings, size_string, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, counts, size_int, cudaMemcpyHostToDevice);
	cudaMemcpy(d_c, merged, size_int, cudaMemcpyHostToDevice);
	//cudaMemcpy(left, any_left, sizeof(int), cudaMemcpyHostToDevice);
	
	int threads_num = 512, blocks_num;
	blocks_num = (int)ceil((float)actual_count/threads_num);
	
	do {
	populate<<<blocks_num, threads_num>>>(d_b, copy_db, d_c, actual_count,left); 	
	select<<<blocks_num, threads_num>>>(copy_db, actual_count);	
	search<<<blocks_num, threads_num>>>(d_b, d_c, actual_count);
	compare<<<blocks_num, threads_num>>>(d_a, d_b, d_c, actual_count, lenString, threshold);
        anyLeft<<<blocks_num, threads_num>>>(d_c, left, actual_count);
	cudaMemcpy(any_left, left, sizeof(int), cudaMemcpyDeviceToHost);
	} while (*any_left == 0);
		
	cudaMemcpy(counts, d_b, size_int, cudaMemcpyDeviceToHost);
	cudaMemcpy(merged, d_c, size_int, cudaMemcpyDeviceToHost);
	
	int counter=0;
	for(int i =0; i<actual_count;i++) {
		if (merged[i] == 2)
			counter++;
	}

	printf("%d\n", counter);
	
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	free(strings);
	free(counts);
	free(merged);
}
