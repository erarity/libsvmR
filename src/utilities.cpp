#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "svm.h"
#include <sys/time.h>

#define VLENGTH1 32		//500
#define VLENGTH2 2500
#define VLENGTH3 5000

//Globals (Parameters)
int VL;
float per_sparsity;
float per_skips;
int numRows;
int numReps;
//function pointer to test type

//Globals (Internals)
svm_node** data_set;
int* skips;
int sparsityLength;
int skipsLength;

//Timing variables
static struct timeval tv;
static struct timezone tz;
static double time_start, time_end;

//****************************\\
//       Timing Methods       \\
//****************************\\

void start_time(){
	gettimeofday(&tv, &tz);
	time_start = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

double elapsed_time(){
	gettimeofday(&tv, &tz);
	time_end   = (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
	return time_end - time_start;
}

//****************************\\
//       Utility Methods      \\
//****************************\\

//Check if a value is contained within an array.
int contains(int val, int* columns, int numColumns){
	for(int i = 0; i < numColumns; i++){
		if(columns[i] == val){
			return 1;
		}
	}
	return 0;
}

// Compare two integers pointed to by x and y; used in invocation of
// the qsort() routine for sorting indices
int compare_ints(const void *x, const void *y){
  return *((int*) x) - *((int*) y);
}

//Remove a value from an array
void remove_val(svm_node* px, int idx, int len){

	// printf("\nList before removal:\n[");
	// for(int i = 0; px[i].index!=-1; i++){
	// 	printf("%03d:%05.3f, ",px[i].index, px[i].value);
	// }
	// printf("]\n");

  for(int i = idx; px[i].index != -1; i++){
    px[i] = px[i+1];
  }

	//svm_node* after_removal = realloc(px, (array_length - 1) * sizeof(svm_node*));

	// printf("List after removal:\n[");
	// while(px->index != -1){
	// 	printf("%03d:%05.3f, ",px->index, px->value);
	// 	++px;
	// }
	// printf("]\n\n");
}

int get_length(svm_node* px){
	if(px[0].index == -1){
		return 0;
	}

	int len = 0;
	while(px->index != -1){
		len++;
		px++;
	}
	return len+1; //Append the dummy node to the length for removal purposes.
}

//****************************\\
//       Testing Methods      \\
//****************************\\

//Used to zero out a column in the data pool, effectively eliminating it logically.
int zero_columns(svm_node *px, int *columns, int numColumns){
	while(px->index != -1){
		if(contains(px->index,columns,numColumns)){
			//printf("Zeroing out index %d.\n",px->index);
			px->value = 0;
			++px;
		}else{
			++px;
		}
	}
	return 1;
}
	
	
//Remove a value from an array
void shift_left_from(svm_node* px, int idx){
  for(int i = idx; px[i].index != -1; i++){
    px[i] = px[i+1];
  }
}


//Used to remove a column from the data pool, eliminating it both logically and virtually.
svm_node* remove_columns(svm_node *px, int *skips, int skip_len){
  int idx = 0;
  //int removal_idx = 0;
  int skip_idx = 0;
  while(px[idx].index != -1 && skip_idx < skip_len){
    if(px[idx].index == skips[skip_idx]){
      shift_left_from(px,idx);
    }
    else if(px[idx].index < skips[skip_idx]){
      idx++;
    }
    else{
      skip_idx++;
    }
  }
  return px;
}


//   int array_length = get_length(px);
//   //Stores the indices for removal during the second pass.
//   //int* removeList = (int*)malloc(numColumns * sizeof(int*));

//   while(px->index != -1){
//     if(contains(px->index,columns,numColumns)){
//       // printf("Removing index %d from the array,\n", px->index);
			
//       //Store the current index in the removeList
//       //removeList[removal_idx] = px->index;
//       remove_val(px, idx, array_length);
//       array_length--;
			
//       ++px;;
//       idx++;
//       //removal_idx++;			
//     }else{
//       ++px;
//     }	
//   }

//   return px;
// }

//****************************\\
//    Fisher-Yates Shuffle    \\
//****************************\\

void shuffle(int* indices, int LEN){
	
	//Fisher-Yates variables.
	// int max = LEN;	//Current logical "maximum" of the array.
	int ran = 0;	//Random index to swap with max.

	//Populate the indices array.
	for(int i = 0; i < LEN; i++){
		indices[i] = i+1;
	}
        // printf("Unshuffled array:\n");
	// for(int i = 0; i < LEN; i++){
        //   printf("%d:%d\n",i,indices[i]);
	// }
        // printf("\n");

        // // Kauffman-Karypis shuffling
        // int shuffles = 3 * LEN;
	// for(int s=0; s<shuffles; s++){
        //   int i = rand() % LEN;
        //   int j = rand() % LEN;
        //   int tempi = indices[i];
        //   indices[i] = indices[j];
        //   indices[j] = tempi;
        // }

	//Fisher-Yates shuffle.	
        int max = LEN-1;
        while(max > 1){
          ran = rand() % (max);
          int temp = indices[ran];
          indices[ran] = indices[max];
          indices[max] = temp;
          max--;
	}

	//printf("Indices: \n[");
	//for(int i = 0; i < LEN; i++){
	//	printf("%03d, ",indices[i]);
	//}
	//printf("]\n");

        // printf("Shuffled array:\n");
	// for(int i = 0; i < LEN; i++){
        //   printf("%d:%d\n",i,indices[i]);
	// }
        // printf("\n");


}

//****************************\\
//    Generating Functions    \\
//****************************\\

//Generates an array of svm_nodes with increasing index and random value (0-99)
svm_node* generate_svm_node_array(int LEN){

        svm_node* nodes = (svm_node*)malloc((LEN+1)   * sizeof *nodes);	//Array to be returned containing a sorted list of svm_nodes.
	int* indices = (int*)malloc(VL    * sizeof(int*));		//Array that is populated with the indices for the fake svm_nodes*
	// int* working = (int*)malloc(LEN-1 * sizeof(int*));		//Working array for applying mergesort to indices.

	//Shuffle the indices
	shuffle(indices, VL);

	//Sort indices by ascending order.
        qsort(indices, LEN, sizeof(int), compare_ints);

	//Transfer the values from indices[0] to indices[LEN-1] over to nodes
	for(int i = 0; i < LEN; i++){
		nodes[i].index = indices[i];
		nodes[i].value = rand() % (99 + 1 - 0) + 0;
	}

	//Append the final node with the (-1,-1) pair to denote the end.
	nodes[LEN].index = -1;
	nodes[LEN].value = -1;

	return nodes;
}

//Generate an array of skip infices with increasing index.
int* generate_skip_list(int LEN){
	
	int* skips = (int*)malloc(VL * sizeof(int*));

	shuffle(skips, VL);

        qsort(skips, LEN, sizeof(int), compare_ints);

	return skips;
}

//****************************\\
//        Dot Functions       \\
//****************************\\

double dot_skip1(const svm_node *px, const svm_node *py, int *skips, int numSkips) {
	double sum = 0;
	int skip_idx = 0;
	int use_skips = numSkips > 0;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			
			while(use_skips && (skip_idx < numSkips) && (skips[skip_idx ] < px->index)){
				skip_idx++;
			}
			if(use_skips && skip_idx >= numSkips){
				use_skips = 0;
			}
			if(!use_skips || skips[skip_idx] != px->index){
				sum += px->value * py->value;
				++px;
				++py;
			}
			else{
				// skipping
				++px;
			       	++py;
			       	skip_idx++;
			}
		}

		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}		
	}
	return sum;
}

double dot_skip2(const svm_node *px, const svm_node *py, int *skips, int numSkips) {
	double sum = 0;
	int skip_idx = 0;
	int use_skips = numSkips > 0;
	numSkips = numSkips - 1;
	while(px->index != -1 && py->index != -1)
	{
		if(px->index == py->index)
		{
			while(use_skips && (skip_idx < numSkips) && (skips[skip_idx ] < px->index)){
				skip_idx++;
			}
			
			if(!use_skips || skips[skip_idx] != px->index){
				sum += px->value * py->value;
				++px;
				++py;
			}
			else{
				++px;
				++py;
			}
		}
		else
		{
			if(px->index > py->index)
				++py;
			else
				++px;
		}		
	}
	return sum;
}

double dot_skip3(const svm_node *px, const svm_node *py, int *skips, int numSkips) {
	double sum = 0;
	int skip_idx = 0;
	int use_skips = numSkips > 0;
	int pxi = px->index;
	int pyi = py->index;
	numSkips = numSkips - 1;
	while(pxi != -1 && pyi != -1)
	{
		if(pxi == pyi)
		{
			while(use_skips && (skip_idx < numSkips) && (skips[skip_idx ] < px->index)){
				skip_idx++;
			}

			if(!use_skips || skips[skip_idx] != px->index){
				sum += px->value * py->value;
				++px;
				++py;
			}else{
				++px;
				++py;			
			}
		}
		else
		{
			if(pxi > pyi)
				++py;
			else
				++px;
		}		
	}
	return sum;
}

// double dot(const svm_node *px, const svm_node *py, int *skips, int numSkips){
	
//   double sum = 0;
//   int skip_idx = 0;
//   int use_skips = numSkips > 0;

//   printf("numSkips %d\n",numSkips);
//   while(px->index != -1 && py->index != -1){
//     if((skip_idx < numSkips) && (skips[skip_idx] < px->index || skips[skip_idx] < py->index)){
//       skip_idx++;
//     }else if(px->index == py->index){
//       printf("px->index %d py->index %d skip_idx %d\n",px->index,py->index,skip_idx);
//       if(px->index != skips[skip_idx]){
//         sum+= px->value * py->value;
//       }
//       ++px;
//       ++py;
//     }else{
//       if(px->index > py->index){
//         ++py;
//       }else{
//         ++px;
//       }
//     }
//   }
// }

double dot_old(const svm_node *px, const svm_node *py){
	
	double sum = 0;
	while(px->index != -1 && py->index != -1){
		if(px->index == py->index){
			sum += px->value * py->value;
			++px;
			++py;
		}else{
			if(px->index > py->index){
				++py;
			}else{
				++px;
			}
		}
	}
	return sum;
}

//****************************\\
//     Printing Functions     \\
//****************************\\

//Printing function for displaying the data set in a formatted way.
void print_data_set(void){
	//Print the data set.
	printf("\n\nDataset:\n");
	for(int i = 0; i < numRows; i++){
		svm_node *row = data_set[i];
		for(int j = 0; j < sparsityLength; j++){

			printf("%03d:%05.2f ", row[j].index, row[j].value);
		}
		printf("\n");
	}
	printf("\n\n");
  return;
}

//Printing function for displaying the skiplist in a formatted way.
void print_skip_list(void){
	//Print the skips list.
	printf("\n\nSkip list:\n");
	for(int i = 0; i < skipsLength; i++){
		printf("%03d ",skips[i]);
	}
	printf("\n\n");	
}

#define RAND_SEED 123456789


int main(int argc, char **argv){

  	srand(RAND_SEED);

	//********************************
	// 	    Parameters
	//********************************
	//  [0] VL 		- virtual length
	//  [1]	per_sparsity 	- how much of the virtual length is sparse
	//  [2]	per_skips	- how much of the virtual length is skipped
	//  [3]	numRows		- how many times an svm_node_array is made (rows)
	//  [4]	numReps		- how many times the algorithm is run on the dot method
	//  [5]	dot_method	- which type of method is used (optional)
	
	//Parse the input.
	if(argc < 6){
		//printf("Too few arguments. Use './<Program> <VirtualLength> <SparsityPercentage> <SkipsPercentage> <NumberOfRows> <NumberOfRepetitions> <DotMethod>'\n");
		printf("Too few arguments. Use './<Program> <VirtualLength> <SparsityPercentage> <SkipsPercentage> <NumberOfRows> <NumberOfRepetitions>'\n");
		exit(EXIT_FAILURE);
	}else if(argc == 6){
		// printf("Proper input received!\n");

	//TODO: Implement a possible function pointer/enum that would allow for a specified testing method.
	//}else if(argc == 7){
	//	printf("Proper input received! (With optional argument for dot_method)\n");
	
	}else{
		//printf("Too many arguments. Use './<Program> <VirtualLength> <SparsityPercentage> <SkipsPercentage> <NumberOfRows> <NumberOfRepetitions> <DotMethod>'\n");
		printf("Too many arguments. Use './<Program> <VirtualLength> <SparsityPercentage> <SkipsPercentage> <NumberOfRows> <NumberOfRepetitions>'\n");
		exit(EXIT_FAILURE);
	}

	//Convert the program arguments.
	VL 		= atoi(argv[1]);
	per_sparsity 	= atof(argv[2]);
	per_skips 	= atof(argv[3]);
	numRows		= atoi(argv[4]);
	numReps		= atoi(argv[5]);
	//dot_method 	= atoi(argv[6]);

        if(per_sparsity > 1.0 || per_skips > 1.0){
          printf("Sparsity and skips must be <= 1.0\n");
          return 0;
        }

	//VL = VLENGTH1;
	//per_sparsity = 0.50;
	//per_skips = 0.10;
	//numRows = 16;
	//numReps = 1;
	
	//Identify which method of dot to use
	//TODO: Implement a function poitner that is later invoked in the dot testing.
	//switch(dot_method){
	//	case 1:
	//		
	//}
	
	//Calculate the proper length of the adjusted skips length.
	float skips_perf = VL * per_skips;
	skipsLength      = (skips_perf - floor(skips_perf) > 0.5) ? ceil(skips_perf) : floor(skips_perf); 

	//Construct the skip list.
	skips = generate_skip_list(skipsLength);

	//Allocate space for the data set.
	data_set = (svm_node**)malloc(numRows * sizeof(svm_node**));

	//Calculate the proper length of the svm_node arrays.
	float sparsityFloat = VL * per_sparsity;
	sparsityLength      = (sparsityFloat - floor(sparsityFloat) > 0.5) ? ceil(sparsityFloat) : floor(sparsityFloat);

	//Construct the data set.
	for(int i = 0; i < numRows; i++){
		//TODO: Implement a enum that allows for the selection of square or
		//	"jagged" data tables.
		svm_node* row = generate_svm_node_array(sparsityLength);
		data_set[i] = row;
	}

	print_skip_list();
	print_data_set();

	//Testing
	
	//Start the timer
	double sum1 = 0.0;
	start_time();
	for(int r = 1; r <= numReps; r++){
		//FIXME: Possible error here? .txt specifies i being set to 1.
		for(int i = 0; i < numRows; i++){
			//FIXME: Possible error here? .txt specifies j being set to 1 and being less than number of columns.
			for(int j = 0; j < numRows; j++){
                          double product = dot_skip1(data_set[i], data_set[j],skips,skipsLength);
			  //printf("Sum 1 (Ver A) so far: %f\n", product);
		  	  sum1 += product;	  
			}
		}
	}
	//End timer. (Report time)
	double skip1_time = elapsed_time();
	printf("Sum 1 (Ver A) total: %f\n", sum1);

	//Start the timer
	float sum2 = 0.0f;
	start_time();
	for(int r = 1; r <= numReps; r++){
		//FIXME: Possible error here? .txt specifies i being set to 1.
		for(int i = 0; i < numRows; i++){
			//FIXME: Possible error here? .txt specifies j being set to 1 and being less than number of columns.
			for(int j = 0; j < numRows; j++){
                          float product = dot_skip2(data_set[i], data_set[j],skips,skipsLength);
			  //printf("Sum 1 (Ver B) so far: %f\n", product);
			  sum2 += product;	
			}
		}
	}
	//End timer. (Report time)
	double skip2_time = elapsed_time();
	printf("Sum 2 (Ver B) total: %f\n", sum2);

	/*	
	//Start the timer
	float sum5 = 0.0f;
	start_time();
	for(int r = 1; r <= numReps; r++){
		//FIXME: Possible error here? .txt specifies i being set to 1.
		for(int i = 0; i < numRows; i++){
			//FIXME: Possible error here? .txt specifies j being set to 1 and being less than number of columns.
			for(int j = 0; j < numRows; j++){
                          float product = dot_skip3(data_set[i], data_set[j],skips,skipsLength);
			  //printf("Sum 1 (Ver B) so far: %f\n", product);
			  sum5 += product;	
			}
		}
	}
	//End timer. (Report time)
	double skip3_time = elapsed_time();
	printf("Sum 3 (Ver C) total: %f\n", sum5);
	*/
	

	// printf("\nTime taken for test 1: %lf\n",skip_list_time);

	//Zero out the columns
	for(int i = 0; i < numRows; i++){
		svm_node* px = data_set[i];	
		zero_columns(px,skips,skipsLength);
	}

	//Start timer again.
	float sum3 = 0.0f;
	start_time();
	for(int r = 1; r <= numReps; r++){
		//FIXME: Possible error here? .txt specifies i being set to 1.
		for(int i = 0; i < numRows; i++){
			//FIXME: Possible error here? .txt specifies j being set to 1 and being less than number of columns.
			for(int j = 0; j < numRows; j++){
				float product = dot_old(data_set[i], data_set[j]);
				//printf("Sum 2 so far: %f\n", product);
				sum3 += product;
			}
		}
	}
	//End timer. (Report time)
	double zeroed_node_time = elapsed_time();
	printf("Sum 3 total: %f\n", sum3);
	
	// printf("\nTime taken for test 2: %lf\n",zeroed_node_time);

	//Remove the columns
	for(int i = 0; i < numRows; i++){
		svm_node* px = data_set[i];
		remove_columns(px,skips,skipsLength);
	}
	
	float sum4 = 0.0f;
	start_time();
	for(int r = 1; r <= numReps; r++){
		//FIXME: Possible error here? .txt specifies i being set to 1.
		for(int i = 0; i < numRows; i++){
			//FIXME: Possible error here? .txt specifies j being set to 1 and being less than number of columns.
			for(int j = 0; j < numRows; j++){
				float product = dot_old(data_set[i], data_set[j]);
				//printf("Sum 3 so far: %f\n", product);
				sum4 += product;
			}
		}
	}
	//End timer. (Report time)
	double removed_node_time = elapsed_time();
	printf("Sum 4 total: %f\n", sum4);

	// printf("\nTime taken for test 3: %lf\n",removed_node_time);

        printf("HEADER : ");
        printf("%6s ","vlen");
        printf("%6s ","sparse");
        printf("%6s ","skips");
        printf("%6s ","rows");
        printf("%6s ","reps");
        printf("%10s ","skip1_time");
        printf("%10s ","skip2_time");
        printf("%10s ","zero_time");
        printf("%10s ","rem_time");
        printf("\n");

        printf("RESULTS: ");
        printf("%6d ",VL);
        printf("%6.2f ",per_sparsity);
        printf("%6.2f ",per_skips);
        printf("%6d ",numRows);
        printf("%6d ",numReps);
        printf("%10.2e ",skip1_time);
        printf("%10.2e ",skip2_time);
        printf("%10.2e ",zeroed_node_time);
        printf("%10.2e ",removed_node_time);
        printf("\n");
}
