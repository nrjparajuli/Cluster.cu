#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <cstring>


int lenString=594;
int maxNumStrings = 1000000;                       
int main(int argc, char** argv) {
	
//allocation of variables
	char *strings,*counts;	//host copy of a
	char *d_a; //device copy of a
	char copy[lenString+1]; //string to copy in info
	int i=0, k=0;
	int size_string = maxNumStrings*sizeof(char)*(lenString+1);
	int size_int = maxNumStrings*sizeof(int)
	struct timeval start, end; 				//using time
	double wallTime;
	cudaError_t status = (cudaError_t)0;



	//opening the file
	FILE *fp;
	//fp=fopen("/cluster/home/charliep/courses/cs360/single-linkage-clustering/Iceland2014.trim.contigs.good.unique.good.filter.unique.count.fasta", "r");


	if (!(strings= (char *)malloc(size_string) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);}
	if (!(counts= (int *)malloc(size_int))) {
		fprintf(stderr, "malloc() FAILED (Block)\n"); 
		exit(0);}


	while( fscanf(fp,"%s %d", copy, count) != EOF){
		strncpy(&strings[i],copy,lenString);
		b[k]=count
		//printf("%s\n", copy);
		//printf("%s\n", &a[i]);
		i=i+lenString;
		k++;
		
		}

