/*!
@file testrw.c -- test reading of ana file

Test program to investigage memory leaks. Seems to work when we simply read a 
file to memory and free the memory directly afterwards.

Compile as:
	gcc -g -Wall -o testrw testrw.c anarw.c anacompress.c anadecompress.c
*/

#include <stdlib.h>  
#include <stdio.h>  
#include <stdint.h>

#include "types.h"
#include "anarw.h"

const int NITER = 5;

int main(int argc, char *argv[]) {
	// Function arguments
	char *filename = argv[1];
	int debug=0;
	// Init ANA IO variables
	char *header = NULL;			// ANA header (comments)
	uint8_t *anaraw = NULL;			// Raw data
	int	nd=-1, type=-1, *ds, size=-1, d; // Various properties
	// Data manipulation
	
	// Read ANA file
	printf("testrw.c: Reading in ANA file a few times\n");
	for (d = 0; d<NITER; d++) {
		printf("iter %d\n", d);
		anaraw = ana_fzread(filename, &ds, &nd, &header, &type, &size);
		free(header);
		free(ds);	
		free(anaraw);
	}
		
	
	return 0;
}
