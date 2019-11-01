#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <immintrin.h>

#include "pagerank.h"

double* iterate(double* matrix, double* pmatrix, size_t npages){
	double* ponematrix =  calloc(npages, sizeof(double));
	for (int j = 0; j < npages; j++){

		//ponematrix[j] = 0;
		int work = j*npages;
		for (int k = 0; k < npages; k++){
			ponematrix[j] += matrix[work+ k]*(pmatrix[k]);
		}
	}

	return ponematrix;

}
//CORRECT
double con(double* pmatrix, double* p1matrix, size_t npages){
	double tot = 0;
	for (int i = 0; i < npages; i++){
		tot += pow(p1matrix[i]-pmatrix[i], 2);
	}
	return sqrt(tot);
}

void pagerank(node* list, size_t npages, size_t nedges, size_t nthreads, double dampener) {

	/*
		TODO

		- implement this function
		- implement any other necessary functions
		- implement any other useful data structures
	*/

		node* current = list;
		node* inlinkcurrent;


		double* matrix= calloc(npages*npages, sizeof(double));

		while (current != NULL){
			inlinkcurrent = current->page->inlinks;

			if (current->page->noutlinks == 0){
				for (int i = 0; i < npages; i++){
					matrix[i * npages + current->page->index] = (double) 1.0/npages;
				}

			}

			while (inlinkcurrent != NULL){

				matrix[inlinkcurrent->page->index + (npages*current->page->index)] = (double) 1.0/inlinkcurrent->page->noutlinks;
				inlinkcurrent = inlinkcurrent->next;
			}


			current = current->next;

		}

		int space = npages*npages;
		double damp = ((1.0-dampener)/npages);
		for (int i = 0; i <space; i ++){
			matrix[i] = matrix[i]*dampener + damp;
		}

		double* pzeromatrix =  calloc(npages, sizeof(double));
		// for (int j = 0; j < npages; j++){
		for (int j = npages-1; j >= 0; j--){
			pzeromatrix[j] = (double) 1.0/npages;
		}


		double* ponematrix = iterate(matrix, pzeromatrix, npages);

		double tot = 0;
		for (int i = 0; i < npages; i++){
			tot += pow(ponematrix[i]-pzeromatrix[i], 2);
		}
		double converge = con(pzeromatrix,ponematrix, npages);

		free(pzeromatrix);
		pzeromatrix = ponematrix;


		while (converge > 0.005){
			ponematrix = iterate(matrix, pzeromatrix, npages);
			converge = con(pzeromatrix, ponematrix, npages);
			free(pzeromatrix);
			pzeromatrix = ponematrix;


		}

		current = list;

		for (int i = 0; i < npages; i++){
			printf("%s %.8lf\n", current->page->name, ponematrix[i]);
			current = current->next;
		}
		free(matrix);
		free(pzeromatrix);


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
