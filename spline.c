/* gvmath.c */

/* ------------------------------ */
/* Feature test macros            */
/* ------------------------------ */
#define	_POSIX_SOURCE						/* Always require POSIX standard */

/* ------------------------------ */
/* Standard include files         */
/* ------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <float.h>

/* ------------------------------ */
/* Local include files            */
/* ------------------------------ */
#define TFOC_CODE
#include "tfoc.h"
#include "gcc_help.h"

/* ------------------------------- */
/* My local typedef's and defines  */
/* ------------------------------- */
#define	REAL_MAX	(FLT_MAX)			/* Whatever REAL is set to */

#define	TRUE	(1)
#define	FALSE	(0)

#define	GVTrimToReal(x)		(x)
#define	ERRprintf(str)			fprintf(stderr,str)

typedef struct _spline_element {
	double cf[3];							/* Coefficients of the spline */
	double x,y;								/* x,y values at the starting knot */
} SPLINE;

/* ------------------------------- */
/* My external function prototypes */
/* ------------------------------- */

/* ------------------------------- */
/* My internal function prototypes */
/* ------------------------------- */

/* ------------------------------- */
/* My usage of other external fncs */
/* ------------------------------- */

/* ------------------------------- */
/* Locally defined global vars     */
/* ------------------------------- */


/* ===========================================================================
-- Heap sort for spline coefficients 
=========================================================================== */
static void spl_sort(SPLINE *spl, int npt) {

	int i,j,hire,ir;
	SPLINE tmp;

/* s is the array being sorted (normally X).  This allows treating all arrays
   symmetrically in the copy operations.  Note that the "s" array is always
	left in a sorted order */
	hire = npt/2;												/* Center of data set	*/
	ir   = npt-1;												/* On the way down		*/

	while (TRUE) {												/* Repeat forever			*/
		if (hire > 0) {										/* Hiring, walk down		*/
			hire--;
			tmp = spl[hire];
		} else {
			tmp = spl[ir];
			spl[ir] = spl[0];
			ir--;
			if (ir == 0) {
				spl[0] = tmp;
				return;
			}
		}

		i = hire;
		j = 2*hire+1;

		while (TRUE) {
			if (j > ir) {								/* Last element? */
				spl[i] = tmp;
				break;
			}
			if (j < ir) {								/* Repeat to put this in place */
				if (spl[j].x < spl[j+1].x) j++;
			}
			if (tmp.x < spl[j].x) {
				spl[i] = spl[j];
				i = j;
				j = 2*j+1;
			} else {
				j = ir+1;
			}
		}
	}
}


/* ============================================================================
-- ... Simple spline determination
--
-- Usage:   void *GVFitSpline(void *work, double *x, double *y, int npt, int opts);
--
-- Inputs:  X    - X coordinates - normally strictly sorted, but not always
--          Y    - Y coordinates
--          npt  - number of elements in x and y.  Must be greater than 3
--          work - pointer to memory block to receive SPLINE description
--                 structures.  Basically knots and coefficients.
--          opts - bit-wise options
--                   0x01 -> do piecewise linear instead of cubic, but
--                           keep full structure so can be used elsewhere.
--
-- Output:  work - filled with series of SPLINE structures with X,Y and
--                 cubic coefficients.  This array of structures is to be
--                 passed to evaluation the spline.
--
-- Returns: pointer to work NULL, or an malloc'd space for coefficients.  On
--          error, returns NULL.
--
-- Notes: (1) The work array (or the malloc'd replacement) contains the
--            spline coefficients and data.  Each element of the array
--            points to a structure containing:
--              struct _spline_element {
--					    double x, y, c[3];
--  				  }
--            where the value of the spline at T Is 
--		           s(T) = (((elem.c[2]*d + elem.c[1])*d + elem.c[0])*d) + elem.y;
--					  elem.x <= T < (elem+1).x
--               d = T-elem.x
--        (2) The X value of the last point is reset to be +REAL_MAX so that
--            any X given to FitEvalSpline will correctly stop.
============================================================================ */
void *GVFitSpline(void *work, double *x, double *y, int npt, int opts) { 

	SPLINE *spl;
	double t1,t2,tm0,tm1,tm2,tmp1,tmp2;			/* Bunch of temp variables */
	int bad_sort=FALSE;
	int i;

/* Validity tests */
	if (npt < 3) return NULL;

/* Create valid pointer to an array of spline structures */
	spl = (SPLINE *) work;
	if (spl == NULL) spl = malloc(npt*sizeof(*spl));

/* First elements */
	for (i=0; i<npt; i++) {									/* Copy over the x,y			*/
		if (i>0 && x[i]<=x[i-1]) bad_sort = TRUE;
		spl[i].x = x[i];
		spl[i].y = y[i];
	}
	if (bad_sort) {
		spl_sort(spl, npt);									/* Try to fix */
		for (i=1; i<npt; i++) { if (spl[i].x<=spl[i-1].x) break; }
		if (i != npt) {
			ERRprintf("ERROR: Data for spline fit must have unique X coordinates.\n");
			if (work == NULL) free(spl);
			return NULL;
		}
	}

/* If faking as piecewise linear, do now */
	if (opts & 0x01) {										/* Piecewise linear only */
		for (i=0; i<npt-1; i++) {
			spl[i].cf[0] = (spl[i+1].y-spl[i].y)/(spl[i+1].x-spl[i].x);	/* Start of RHS solution	*/
			spl[i].cf[1] = spl[i].cf[2] = 0;
		}
		spl[npt-1].x = REAL_MAX;							/* Terminating condition */
		spl[npt-1].cf[0] = spl[npt-2].cf[0];
		spl[npt-1].cf[1] = spl[npt-1].cf[2] = 0;
		return spl;
	}

/* ... Compute not-a-knot spline */
	for (i=1; i<npt; i++) {									/* Fill in matrix first		*/
		spl[i].cf[0] = 0.0f;									/* Just so has some value	*/
		spl[i].cf[1] = spl[i].x-spl[i-1].x;				/* h(j) in Stoer notation	*/
		spl[i].cf[2] = (spl[i].y-spl[i-1].y)/(spl[i].x-spl[i-1].x);	/* Start of RHS solution	*/
	}

/* Continue with not-a-knot spline - complete first block */
	spl[0].cf[2] = spl[2].cf[1];							/* Duplicate 2nd point */
	spl[0].cf[1] = spl[1].cf[1] + spl[2].cf[1];		/* And sum of first two */
	spl[0].cf[0] = ((spl[1].cf[1]+2*spl[0].cf[1])*spl[1].cf[2]*spl[2].cf[1] +
		              spl[1].cf[1]*spl[1].cf[1]*spl[2].cf[2])/spl[0].cf[1];

	tm1 = spl[npt-1].cf[1];
	tm2 = spl[npt-1].cf[2];

	for (i=1; i<npt-2; i++) {
		t1 = -spl[i+1].cf[1]/spl[i-1].cf[2];
		spl[i].cf[0] = GVTrimToReal(t1*spl[i-1].cf[0] + 3*(spl[i].cf[1]*spl[i+1].cf[2] + spl[i+1].cf[1]*spl[i].cf[2]));
		spl[i].cf[2] = GVTrimToReal(t1*spl[i-1].cf[1] + 2*(spl[i].cf[1]+spl[i+1].cf[1]));
	}

	t1  = -tm1/spl[npt-3].cf[2];
	spl[npt-2].cf[0] = GVTrimToReal(t1*spl[npt-3].cf[0] + 3*(spl[npt-2].cf[1]*tm2 + tm1*spl[npt-2].cf[2]));
	spl[npt-2].cf[2] = GVTrimToReal(t1*spl[npt-3].cf[1] + 2*(spl[npt-2].cf[1]+tm1));
	t1  = spl[npt-2].cf[1] + tm1;
	tm0 = ((tm1+2*t1)*tm2*spl[npt-2].cf[1]+tm1*tm1*(spl[npt-2].y-spl[npt-3].y)/spl[npt-2].cf[1])/t1;
	t1  = -t1/spl[npt-2].cf[2];
	tm2 = spl[npt-2].cf[1];
	tm2 = t1*spl[npt-2].cf[1]+tm2;
	tm0 = (t1*spl[npt-2].cf[0]+tm0)/tm2;
	spl[npt-2].cf[0] = GVTrimToReal((spl[npt-2].cf[0]-spl[npt-2].cf[1]*tm0)/spl[npt-2].cf[2]);

	for (i=npt-3; i>=0; i--) 
		spl[i].cf[0] = GVTrimToReal((spl[i].cf[0]-spl[i].cf[1]*spl[i+1].cf[0]) / spl[i].cf[2]);

	for (i=1; i<npt-1; i++) {
		t2   = spl[i].cf[1];
		tmp1 = (spl[i].y-spl[i-1].y)/t2;
		tmp2 = spl[i-1].cf[0]+spl[i].cf[0]-2*tmp1;	
		spl[i-1].cf[1] = GVTrimToReal((tmp1-spl[i-1].cf[0]-tmp2)/t2);
		spl[i-1].cf[2] = GVTrimToReal(tmp2/t2/t2);
	}

	t2   = tm1;
	tmp1 = (spl[npt-1].y-spl[npt-2].y)/t2;
	tmp2 = spl[npt-2].cf[0]+tm0-2*tmp1;
	spl[npt-2].cf[1] = GVTrimToReal((tmp1-spl[npt-2].cf[0]-tmp2)/t2);
	spl[npt-2].cf[2] = GVTrimToReal(tmp2/t2/t2);

	spl[npt-1].x = REAL_MAX;							/* Terminating condition */
	return spl;
}


/* ============================================================================
-- Evaluate a spline function
--
-- Usage: void *GVEvalSpline(spl, x);
--
-- Inputs: x   - x value to evaluate spline
--         spl - pointer to workspace with spline data and coefficients
--
-- Output: none
--
-- Return: Returns value of spline at specified point
============================================================================ */
double GVEvalSpline(void *work, double x) {

	SPLINE *spl;

	spl = (SPLINE *) work;						/* Just rename it		 */
	if (x <= spl[0].x) return spl[0].y;		/* Constant extension */

	while (x > spl[1].x) spl++;				/* Should terminate by REAL_MAX */
	if (spl[1].x == REAL_MAX) return spl[0].y;

	x = x - spl->x;								/* Distance from knot */
	return ( ((spl->cf[2]*x + spl->cf[1])*x + spl->cf[0])*x + spl->y );
}
