#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "svm.h"

#define ARRAYLENA 12
#define ARRAYLENB 14

#define VLENGTH1 500
#define VLENGTH2 1000
#define VLENGTH3 2500

//Check if a value is contained within an array.
int contains(int val, int* columns, int numColumns){
	for(int i = 0; i < numColumns; i++){
		if(columns[i] == val){
			return 1;
		}
	}
	return 0;
}

//Used to zero out a column in the data pool, effectively eliminating it logically.
int zero_columns(svm_node *px, int *columns, int numColumns){
	
	int skip_idx = 0;
	int use_skips = numColumns > 0;

	while(px->index != -1){
		if(contains(px->index,columns,numColumns)){
			printf("Zeroing out index %d.\n",px->index);
			px->value = 0;
			++px;
		}else{
			++px;
		}
	}
	return 1;
}
	
	
//Used to remove a column from the data pool, eliminating it both logically and virtually.
svm_node* remove_columns(svm_node *px, int *columns, int numColumns){

	int idx = 0;
	int removal_idx = 0;
	int skip_idx = 0;

	//Stores the indices for removal during the second pass.
	int* removeList = (int*)malloc(numColumns * sizeof(int*));

	while(px->index != -1){
		if(contains(px->index,columns,numColumns)){
			printf("Adding index %d to the removal list.\n", px->index);
			
			//Store the current index in the removeList
			removeList[removal_idx] = px->index;
			
			++px;;
			idx++;
			removal_idx++;			
		}else{
			++px;
		}	
	}
	free(removeList);
	return px;
}

//Returns an array that is a sample of the full data set.
svm_node* sample_set(svm_node *px, int len, float percentage){

	//Calculate the desired sample size.
	float sampleFloat = len * percentage;
	int sampleSize = (sampleFloat - floor(sampleFloat) > 0.5) ? ceil(sampleFloat) : floor(sampleFloat);
	
	//Malloc a new array the length of the sample size.
	svm_node* sample = (svm_node*)malloc(sampleSize * sizeof *sample);

	for(int i = 0; i < sampleSize; i++){
		//Begin modifying the data within the array.
		//The values returned from this can return duplicates.
		int targetIndex = rand() % ((len-1) + 1 - 0) + 0;
		sample[i] = px[targetIndex];
	}

	return sample;
}

//Generates an array of svm_nodes with increasing index and random value (0-100)
svm_node* generate_svm_node_array(int LEN){

	svm_node* nodes = (svm_node*)malloc(LEN * sizeof *nodes);

	int minVal = 0;
	//Random range is rand() % (max_num + 1 - min_num) + min_num
	int offset = rand() % (10 + 1 - 2) + 2;
	
	for(int i = 0; i < LEN - 1; i++){
		int tempIndex = rand() % ((minVal + offset) + 1 - (minVal + 1)) + (minVal + 1);
		minVal = tempIndex;
		nodes[i].index = tempIndex;
		offset = rand() % (10 + 1 - 2) + 2;
		nodes[i].value = (float)rand()/(float)((RAND_MAX / 100) + 1);
		printf("Adding node with index: %03d, value: %05.2f\n", nodes[i].index, nodes[i].value);
	}

	//Append the final node with the (-1,-1) pair to denote the end.
	nodes[LEN-1].index = -1;
	nodes[LEN-1].value = -1;

	return nodes;
}

int main(void){

	int skips[7] = {1,4,6,8,16,22,31};

	//sample_set(100,0.1);
	//sample_set(100,0.25);
	//sample_set(100,0.50);
	
	svm_node *px = generate_svm_node_array(ARRAYLENA);
	svm_node *py = generate_svm_node_array(ARRAYLENB);

	sample_set(px,100,0.1);
	sample_set(px,100,0.25);
	sample_set(px,100,0.50);

	//Print the first array
	printf("\n\nFirst array:\n[ ");
	//while(px->index != -1){
	for(int i = 0; i < ARRAYLENA; i++){
		printf("%03d:%05.2f, ", px[i].index, px[i].value);
	}
	printf("]\n\n");

	//Print the second array
	printf("\n\nSecond array:\n[ ");
	for(int i = 0; i < ARRAYLENB; i++){
		printf("%03d:%05.2f, ", py[i].index, py[i].value);
	}
	printf("]\n\n");

	//Zero out the columns
	zero_columns(px,skips,7);
	zero_columns(py,skips,7);

	//Print the first array
	printf("\n\nFirst array:\n[ ");
	//while(px->index != -1){
	for(int i = 0; i < ARRAYLENA; i++){
		printf("%03d:%05.2f, ", px[i].index, px[i].value);
	}
	printf("]\n\n");

	//Print the second array
	printf("\n\nSecond array:\n[ ");
	for(int i = 0; i < ARRAYLENB; i++){
		printf("%03d:%05.2f, ", py[i].index, py[i].value);
	}
	printf("]\n\n");

	//Remove the columns
	printf("Removing the columns.\n");
	remove_columns(px,skips,7);
	remove_columns(px,skips,7);


}
