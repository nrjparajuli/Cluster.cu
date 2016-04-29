#include <stdio.h>
#include <stdlib.h>
#include <cmath>

__device__ int position;			//index of the largest value
__device__ int largest;				//value of the largest value
int lenString = 593;				
int maxNumStrings = 1000000;                           
int threshold = 2;

// Checks if there are any unmerged tuples (Parallelized)
__global__ void anyLeft(int *d_c, int *remaining, int size) {
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	if((d_c[my_id] == 0) && (my_id < size)) {
		*remaining = 0;
	}
}

// Searches for the index of the largest count (Parallelized)
__global__ void search(int *d_b, int *d_c, int size) {
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;
	if((d_c[my_id] == 0) && (d_b[my_id] == largest) && (my_id < size)) {
		position = my_id;
	}
}

// Populates copy_db such that the counts for merged tuple is 0
// but count for unmerged tuples is unchanged (Parallelized)
__global__ void populate (int *d_b, int *copy_db, int *d_c, int size, int *left) {
	int n = 0;
	*left = 1;	// reinitalized to false to check if all strings are merged

	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if (my_id < size) {
		n = abs((bool)d_c[my_id] - 1);
		copy_db[my_id] = d_b[my_id] * n;
	}		
}

// Reduction-type tree implementation to find largest count (Parallelized)
__device__ void cuda_select(int *db, int size) {
	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if(my_id < size) {
		if(db[2 * my_id] > db[2 * my_id + 1])
			db[my_id] = db[2 * my_id];
		else
			db[my_id] = db[2 * my_id + 1];
	}	
}

// Loops cuda_select function until largest value is at index 0
__global__ void select(int *db, int size) {
	int height = (int)ceil(log2((double)size));
	int i = 0;
	
	for(i = 0; i < height; i++) {
		size = (int)ceil((double) size/2);
		cuda_select(db, size);
	}
	largest = db[0];
}

// Compares target string to all other unmerged strings with lesser count
__global__ void compare(char *d_a, int *d_b, int *d_c, int size, int lenString, int threshold) {

	int my_id = blockDim.x * blockIdx.x + threadIdx.x;

	if (my_id == position) 
		d_c[my_id] = 2;
	
		
	if ((my_id < size) && (d_c[my_id] == 0) && (my_id != position)) {	
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
	char *strings, *d_a;	// host and device copy of strings
	int *counts, *d_b;	// host and device copy of counts
	int *merged, *d_c;	// host and device copy of bools
	int *copy_db;		// device copy of counts (counts of merged is 0)
	char copy[lenString+1]; // intermediate variable to load strings into array from file
	int numbers;		// intermediate variable to load counts into array from file
	int *any_left, *left;   // host and device copies to check if all tuples are merged
	int size = 0;		// keeps track of number of tuples in the file
	int i = 0;		// loop variable
	int size_string = maxNumStrings*sizeof(char)*(lenString+1);
	int size_int = maxNumStrings*sizeof(int);

	// Open the file
	FILE *fp;
	fp = fopen("/cluster/home/charliep/courses/cs360/single-linkage-clustering/Iceland2014.trim.contigs.good.unique.good.filter.unique.count.fasta", "r");

	// Allocate space for arrays on the host
	if (!(strings = (char *)malloc(size_string))) {
		fprintf(stderr, "malloc() FAILED (Block)\n");
		exit(0);
	}

	if (!(counts = (int*)malloc(size_int))) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);
	}

	if (!(merged = (int*)malloc(size_int))) {
                fprintf(stderr, "malloc() FAILED (Block)\n");
                exit(0);
        }
	
	any_left = (int *)malloc(sizeof(int));
	
	// Set the values of global variables on the device
	cudaMemset(&position, 0, sizeof(int));
	cudaMemset(&largest, 0, sizeof(int));

	// Load strings and counts into array
	while(fscanf(fp, "%s %d", copy, &numbers) != EOF && size < 1000){
		strcpy(&strings[i], copy);
		counts[size] = numbers;
		
		i = i + lenString;
		size++;	
	}
	
	// Close file
	fclose(fp);
	
	// Allocate space for arrays on the device
	cudaMalloc(&d_a, size_string);
	cudaMalloc(&d_b, size_int);
	cudaMalloc(&d_c, size_int);
	cudaMalloc(&copy_db, size_int);
	cudaMalloc(&left, sizeof(int));
	
	// Copy arrays from host to device
	cudaMemcpy(d_a, strings, size_string, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, counts, size_int, cudaMemcpyHostToDevice);
	cudaMemcpy(d_c, merged, size_int, cudaMemcpyHostToDevice);

	// Determine number of threads and blocks needed	
	int threads_num = 512, blocks_num;
	blocks_num = (int)ceil((float)size/threads_num);
	
	// Cluster the strings for the given threshold
	do {
	populate<<<blocks_num, threads_num>>>(d_b, copy_db, d_c, size,left); 	
	select<<<blocks_num, threads_num>>>(copy_db, size);	
	search<<<blocks_num, threads_num>>>(d_b, d_c, size);
	compare<<<blocks_num, threads_num>>>(d_a, d_b, d_c, size, lenString, threshold);
        anyLeft<<<blocks_num, threads_num>>>(d_c, left, size);
	cudaMemcpy(any_left, left, sizeof(int), cudaMemcpyDeviceToHost);
	} while (*any_left == 0);
	
	// Copy results back from device to host
	cudaMemcpy(strings, d_a, size_string, cudaMemcpyDeviceToHost);	
	cudaMemcpy(counts, d_b, size_int, cudaMemcpyDeviceToHost);
	cudaMemcpy(merged, d_c, size_int, cudaMemcpyDeviceToHost);
	
	int counter = 0;
	
	FILE *output = fopen("output2.txt", "w+");
	
	for(i = 0; i < size; i++) {
		strncpy(copy, &strings[i*lenString], lenString);
		fprintf(output, "%s %d\n", copy, counts[i]);
		if (merged[i] == 2)
			counter++;
	}
	
	fclose(output);

	printf("%d\n", counter);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);
	cudaFree(copy_db);
	cudaFree(left);

	free(strings);
	free(counts);
	free(merged);
	free(any_left);
}
