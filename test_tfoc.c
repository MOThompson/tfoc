/* ===========================================================================
  Simple test of subroutine-oriented TFOC
=========================================================================== */
#include <stdio.h>

/* ------------------------------ */
/* Local include files            */
/* ------------------------------ */
#include "tfoc.h"

int main(int argc, char *argv[]) {

	printf("%f", TFOC_GetRefl("a.sam", 10600, 75, TE, 300.0).R);
	return 0;
}
