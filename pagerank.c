#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"
#define limit EPSILON * EPSILON

size_t g_chunk;
size_t g_npages;
size_t g_nthreads;


typedef struct {
	double result;
	const size_t id;
	const double* pzero;
	const double* pone;
}arg_con;

typedef struct {
	const double* mmatrix;
	const double* pmatrix;
	double* ponematrix;
}arg_it;

size_t nextRow;


void* workerIt(void* args){
	arg_it* wargs = (arg_it*) args;

	int work;

	while (nextRow < g_npages){
		const size_t start = __sync_fetch_and_add(&nextRow, 1);
		work = g_npages*start;
		double value = 0;
		//printf("doing %zu\n", start);
		for (int k = 0; k < g_npages; k++){
			value += wargs->mmatrix[work + k]*(wargs->pmatrix[k]);
		}
		wargs->ponematrix[start] = value;
	}
	return NULL;

}


double* iterate(double* matrix, double* pmatrix, size_t npages, pthread_t* thread){
	nextRow = 0;
	double* ponematrix =  calloc(npages, sizeof(double));
	arg_it* args = malloc(g_nthreads*sizeof(arg_it));
	for(size_t i = 0; i < g_nthreads; i++){
		args[i] = (arg_it){
			.mmatrix = matrix,
			.pmatrix = pmatrix,
			.ponematrix = ponematrix
		};
		pthread_create(thread + i, NULL, workerIt, args + i);
	}
	for (int i = 0; i < g_nthreads; i++) {
		        // wait for threads to finish
		pthread_join(thread[i], NULL);
	}
	return ponematrix;









	for (int j = 0; j < npages; j++){

		ponematrix[j] = 0;
		int work = j*npages;
		for (int k = 0; k < npages; k++){
			//printf("BEFORE p1 += %.3lf, matrix pos is %.3lf, p0 is %.3lf \n", ponematrix[j], matrix[j*npages + k], (pmatrix[j]));
			ponematrix[j] += matrix[work + k]*(pmatrix[k]);
			//printf("AFTER p1 += %.3lf, matrix pos is %.3lf, p0 is %.3lf \n", ponematrix[j], matrix[j*npages + k], (pmatrix[j]));
		}
	}
	//printf("First Iteration \n");
	// printf("P one matrix: \n");
	// int i = 0;
	// for (int j = 0; j <npages; j++){
	// 	printf("%.3lf ", ponematrix[i * npages + j]);
	// 	printf("\n");
	// }
	// printf("\n");


	// printf("P zero matrix: \n");
	// int i = 0;
	// for (int j = 0; j <npages; j++){
	// 	printf("%.3lf ", pmatrix[i * npages + j]);
	// 	printf("\n");
	// }
	// printf("\n");
	// printf("P one matrix: \n");
	// i = 0;
	// for (int j = 0; j <npages; j++){
	// 	printf("%.3lf ", ponematrix[i * npages + j]);
	// 	printf("\n");
	// }
	// printf("\n");


	return ponematrix;

}


void* workerCon(void* args){
	arg_con* wargs = (arg_con*) args;
	const size_t id = wargs->id;
	const size_t start = id * g_chunk;
	const size_t end = (id == g_nthreads - 1) ? g_npages : (id + 1) * g_chunk;

	double diff = 0;
	const double* pz = wargs->pzero;
	const double* po = wargs->pone;
	// double* pn = wargs->id == g_nthreads-1 ? g_npages : (wargs->id+1)*g_chunk;
	for (size_t i = start; i < end; i++){
		diff += pow(po[i] - pz[i], 2);
	}

	wargs->result = diff;
	return NULL;
}


//CORRECT
double con(double* pmatrix, double* p1matrix, size_t npages, size_t nthreads, pthread_t* thread){
	arg_con* args = malloc(g_nthreads*sizeof(arg_con));
	for(size_t i = 0; i < g_nthreads; i++){
		args[i] = (arg_con){
			.id = i,
			.result = 0,
			.pzero = pmatrix,
			.pone = p1matrix
		};
		pthread_create(thread + i, NULL, workerCon, args + i);
	}

	for (int i = 0; i < g_nthreads; i++) {
		        // wait for threads to finish
		pthread_join(thread[i], NULL);
	}
	double diff  = 0;
	for (int i = 0; i < g_nthreads; i++){
		diff += args[i].result;
	}
	return diff;



	//
	// free(args);


	// double tot = 0;
	// for (int i = 0; i < npages; i++){
	// 	//printf("p1 is %.3lf and p0 is %.3lf\n", p1matrix[i], pmatrix[i]);
	// 	tot += pow(p1matrix[i]-pmatrix[i], 2);
	// }
	// return tot;
}

void pagerank(node* list, size_t npages, size_t nedges, size_t nthreads, double dampener) {
	g_nthreads = (nthreads > npages)? npages : nthreads;
	g_npages = npages;
	g_chunk = npages/g_nthreads;
	pthread_t threads[nthreads];
	/*
		TODO

		- implement this function
		- implement any other necessary functions
		- implement any other useful data structures
	*/

	/*
	############################
	### INITIALISED MATRIX M ###
	############################
	*/

		node* current = list;
		node* inlinkcurrent;


		double* matrix= calloc(npages*npages, sizeof(double));

		while (current != NULL){
			inlinkcurrent = current->page->inlinks;
			//printf("Page: %s\n", current->page->name);

			if (current->page->noutlinks == 0){
				for (int i = 0; i < npages; i++){
					matrix[i * npages + current->page->index] = 1.0/npages;
				}

			}

			while (inlinkcurrent != NULL){
				//printf("In page of %s: %s\n", current->page->name, inlinkcurrent->page->name);
				matrix[inlinkcurrent->page->index + (npages*current->page->index)] = 1.0/inlinkcurrent->page->noutlinks;
				inlinkcurrent = inlinkcurrent->next;
			}


			current = current->next;
			// inlinkcurrent = current->page->inlinks;
		}
		//Prints matrix to debug
		// printf("Matrix: \n");
		// for (int i = 0; i < npages; i++){
		// 	for (int j = 0; j <npages; j++){
		// 		printf("%.3lf ", matrix[i * npages + j]);
		// 	}
		// 	printf("\n");
		// }
		// printf("\n");
		/*
		###########################################
		### INITIALISE MATRIX M HAT ONTO MATRIX ###
		###########################################
		*/

		//dampenered

		//double* hatmatrix = calloc(npages*npages, sizeof(double));
		//double* E = calloc(npages*npages, sizeof(double));
		int space = npages*npages;
		const double register damp = ((1.0-dampener)/npages);
		for (int i = 0; i < space; i ++){
			matrix[i] = matrix[i]*dampener + damp;
		}

		//Prints for debugging purposes.
		// printf("matrix: \n");
		// for (int i = 0; i < npages; i++){
		// 	for (int j = 0; j <npages; j++){
		// 		printf("%.3lf ", matrix[i * npages + j]);
		// 	}
		// 	printf("\n");
		// }
		// printf("\n");

		/*
		######################
		### INITIALISE P 0 ###
		######################
		*/

		double* pzeromatrix = calloc(npages, sizeof(double));
		// for (int j = 0; j < npages; j++){
		for (int j = npages-1; j >= 0; j--){
			pzeromatrix[j] = 1.0/npages;
		}

		/*
		######################
		### INITIALISE P 1 ###
		######################
		*/

		// printf("------\n");
		double* ponematrix = iterate(matrix, pzeromatrix, npages, threads);
	    // printf("------\n");

	   // calloc(npages, sizeof(double));
	   // for (int j = 0; j < npages; j++){
	   //
	   // 	ponematrix[j] = 0;
	   // 	for (int k = 0; k < npages; k++){
	   // 		ponematrix[j] += matrix[j*npages + k]*(pzeromatrix[j]);
	   // 	}
	   // }


		//Prints for debugging purposes.
		// printf("P one matrix: \n");
		// int i = 0;
		// for (int j = 0; j <npages; j++){
		// 	printf("%.3lf ", ponematrix[i * npages + j]);
		// 	printf("\n");
		// }
		// printf("\n");
		/*
		##########################
		### INITIALISE P1 - P0 ###
		##########################
		*/
		double tot = 0;
		for (int i = 0; i < npages; i++){
			tot += pow(ponematrix[i]-pzeromatrix[i], 2);
		}
		double converge = con(pzeromatrix,ponematrix, npages, nthreads, threads);


		// printf("%.3lf\n", converge);

		/*
		###############
		### ITERATE ###
		###############
		*/


		// free(pzeromatrix);


		// int counter = 0;

		while (converge > limit){
			pzeromatrix = ponematrix;
			ponematrix = iterate(matrix, pzeromatrix, npages, threads);
			converge = con(pzeromatrix, ponematrix, npages, nthreads, threads);


			// printf("Con is: %.3lf\n", converge);
			// if (counter == 5){
			// 	printf("Break\n");
			// 	break;
			// }
			// counter++;
			// free(pzeromatrix);
			// pzeromatrix = ponematrix;


		}
		// printf("Con is %.3lf < 0.005\n", converge);
		//
		// printf("Final matrix: \n");
		current = list;

		for (int i = 0; i < npages; i++){
			//printf("%s %.8lf\n", current->page->name, ponematrix[0 * npages + i]);
			printf("%s %.8lf\n", current->page->name, ponematrix[i]);
			current = current->next;
		}
		// free(matrix);
		// free(pzeromatrix);




		// printf("A %.8lf\n", ponematrix[0 * npages + 0]);
		// printf("B %.8lf\n", ponematrix[0 * npages + 1]);
		// printf("C %.8lf\n", ponematrix[0 * npages + 2]);
		// printf("D %.8lf\n", ponematrix[0 * npages + 3]);

}

/*
######################################
### DO NOT MODIFY BELOW THIS POINT ###
######################################
*/

int main(int argc, char** argv) {

	/*
	######################################################
	### DO NOT MODIFY THE MAIN FUNCTION OR HEADER FILE ###
	######################################################
	*/

	config conf;

	init(&conf, argc, argv);

	node* list = conf.list;
	size_t npages = conf.npages;
	size_t nedges = conf.nedges;
	size_t nthreads = conf.nthreads;
	double dampener = conf.dampener;

	pagerank(list, npages, nedges, nthreads, dampener);

	release(list);

	return 0;
}
