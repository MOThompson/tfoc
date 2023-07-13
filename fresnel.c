/* ------------------------------ */
/* Feature test macros            */
/* ------------------------------ */
#define _POSIX_SOURCE						/* Always require POSIX standard */
/*	#include "preload.h" */

/* ------------------------------ */
/* Standard include files         */
/* ------------------------------ */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <ctype.h>

/* ------------------------------ */
/* Local include files            */
/* ------------------------------ */
#define TFOC_CODE
#include "tfoc.h"
#include "gcc_help.h"

/* ------------------------------- */
/* My local typedef's and defines  */
/* ------------------------------- */
#define	panic			SysPanic(__FILE__, __LINE__)
/* #define DEBUG */

typedef struct _M_ARRAY {
	COMPLEX A,B,C,D;
} M_ARRAY;

/* ------------------------------- */
/* My external function prototypes */
/* ------------------------------- */
REFL TFOC_ReflN(double theta, POLARIZATION mode, double lambda, TFOC_LAYER layer[]);
REFL TFOC_Refl (double theta, POLARIZATION mode, double lambda, COMPLEX n0, COMPLEX n1, COMPLEX ns, double z);

COMPLEX CADD(COMPLEX a, COMPLEX b);			/* Used occasionally by other routines */
COMPLEX CSUB(COMPLEX a, COMPLEX b);
COMPLEX CMUL(COMPLEX a, COMPLEX b);
COMPLEX CDIV(COMPLEX a, COMPLEX b);
COMPLEX CCSQRT(COMPLEX r);
COMPLEX CPOW(COMPLEX r, double pow);
double  CABS(COMPLEX a);
COMPLEX CSQRT(double r);

/* ------------------------------- */
/* My internal function prototypes */
/* ------------------------------- */
static REFL my_TFOC_ReflN(double theta, POLARIZATION mode, double lambda, TFOC_LAYER layer[]);
static REFL my_TFOC_Refl(double theta, POLARIZATION mode, double lambda, COMPLEX n0, COMPLEX n1, COMPLEX ns, double z);

static BOOL CalcFresnel(double S, POLARIZATION mode, COMPLEX ni, COMPLEX nj, COMPLEX *rij, COMPLEX *tij);
static M_ARRAY CalcInterface(double S, POLARIZATION mode, COMPLEX ni, COMPLEX nj);
static M_ARRAY CalcGap(double z, double S, COMPLEX ni, double lambda);

static M_ARRAY IDENTITY_MATRIX(void);
static M_ARRAY MATMUL(M_ARRAY *a, M_ARRAY *b);
#if 0													/* Not actually used, so comment out for now */
	static M_ARRAY MATINV(M_ARRAY *a);
#endif

/* ------------------------------- */
/* My usage of other external fncs */
/* ------------------------------- */

/* ------------------------------- */
/* Locally defined global vars     */
/* ------------------------------- */

/* ===========================================================================
-- Routine to print out an array for user	interpretation
=========================================================================== */
void Print_M_Array(M_ARRAY m, char *text) {
#ifdef DEBUG
	printf("%s\n", text);
	printf("%g%+gi\t%g%+gi\n", m.A.x, m.A.y, m.B.x, m.B.y);
	printf("%g%+gi\t%g%+gi\n", m.C.x, m.C.y, m.D.x, m.D.y);
	printf("\n");
#endif
	return;
}

/* ===========================================================================
-- Calculate the reflectance off of a complex stack structure.  The
-- initial and substrate media are specified.  The stack is an array
-- of material and thickness structures, terminated by a NULL material.
=========================================================================== */
REFL TFOC_ReflN(double theta, POLARIZATION mode, double lambda, TFOC_LAYER layer[]) {
	REFL te={0.0,0.0},tm;
	
	switch (mode) {
		case TM:
		case TE:
			te = my_TFOC_ReflN(theta, mode, lambda, layer);
			break;
		case UNPOLARIZED:
			te = my_TFOC_ReflN(theta, TE, lambda, layer);
			tm = my_TFOC_ReflN(theta, TM, lambda, layer);
			te.R = 0.5*(te.R + tm.R);
			te.T = 0.5*(te.T + tm.T);
	}
	return te;
}

static REFL my_TFOC_ReflN(double theta, POLARIZATION mode, double lambda, TFOC_LAYER layer[]) {

	double S, sin_out;
	M_ARRAY Ct, Cij, Ciz;
	COMPLEX ni,nj,one={1.0, 0.0};
	REFL rc;
	int i;

	S = layer[0].n.x*sin(theta*pi/180.0f);				/* S factor */

	Ct = IDENTITY_MATRIX();									/* Make an identify matrix */
	Print_M_Array(Ct, "Ct");

	ni = layer[0].n;											/* Current is incident */
	for (i=1; layer[i].type == SUBLAYER; i++) {
		if (layer[i].z <= 0.0) continue;					/* Not really there	*/
		nj  = layer[i].n;
		Cij = CalcInterface(S, mode, ni, nj);
		Print_M_Array(Cij, "Cij");
		Ciz = CalcGap(layer[i].z, S, nj, lambda);
		Print_M_Array(Ciz, "Ciz");
		Ct = MATMUL(&Ct, &Cij);
		Print_M_Array(Ct, "Ct");
		Ct = MATMUL(&Ct, &Ciz);
		Print_M_Array(Ct, "Ct");
		ni = nj;
	}

/* And add the substrate */
	Cij = CalcInterface(S,mode,ni,layer[i].n);
	Ct = MATMUL(&Ct, &Cij);

/* Calculate the reflectivity */
	rc.R = pow(CABS(CDIV(Ct.C,Ct.A)),2);
	sin_out = S/layer[i].n.x;
	if (sin_out > 1.0 || sin_out < 0.0) {
		rc.T = 0;
	} else {
		rc.T = pow(CABS(CDIV(one, Ct.A)),2) *				/* Electric field term				 */
				 layer[i].n.x / layer[0].n.x  *				/* Correct for index of substrate */
				 sqrt(1.0-sin_out*sin_out) / cos(theta*pi/180.0f);	/* Angle correction */
	}

	return rc;
}


/* ===========================================================================
-- Simple routine to return the reflection off a single layer.  Takes
-- incident medium, substrate medium, film properties, thickness, and
-- initial setup values and determines the reflectance coefficient.
=========================================================================== */
REFL TFOC_Refl(double theta, POLARIZATION mode, double lambda, COMPLEX n0, COMPLEX n1, COMPLEX ns, double z) {
	REFL te={0.0,0.0},tm;

	switch (mode) {
		case TM:
		case TE:
			te = my_TFOC_Refl(theta, mode, lambda, n0, n1, ns, z);
			break;
		case UNPOLARIZED:
			te = my_TFOC_Refl(theta, TE, lambda, n0, n1, ns, z);
			tm = my_TFOC_Refl(theta, TM, lambda, n0, n1, ns, z);
			te.R = 0.5*(te.R + tm.R);
			te.T = 0.5*(te.T + tm.T);
	}
	return te;
}

static REFL my_TFOC_Refl(double theta, POLARIZATION mode, double lambda, COMPLEX n0, COMPLEX n1, COMPLEX ns, double z) {

	double S, sin_out;
	M_ARRAY C01, C1z, C12, Ct;
	COMPLEX one={1.0, 0.0};
	REFL rc;

	S = sin(theta*pi/180.0f);								/* S factor */

/* Calculate the front interface element first */
	C01 = CalcInterface(S, mode, n0, n1);

/* Calculate the Z phase */
	C1z = CalcGap(z, S, n1, lambda);

/* Calculate the amorphous/crystalline now */
	C12 = CalcInterface(S, mode, n1, ns);

/* And do the multiplications */
	Ct = MATMUL(&C01, &C1z);
	Ct = MATMUL(&Ct, &C12);

/* Calculate the reflectivity */
	rc.R = pow(CABS(CDIV(Ct.C,Ct.A)),2);
	sin_out = S/ns.x;
	if (sin_out > 1.0 || sin_out < 0.0) {
		rc.T = 0;
	} else {
		rc.T = pow(CABS(CDIV(one, Ct.A)),2) *			/* Electric field term				 */
				 ns.x / n0.x  *								/* Correct for index of substrate */
				 sqrt(1.0-sin_out*sin_out) / cos(theta*pi/180.0f);	/* Angle correction */
	}

	return rc;
}

/* ===========================================================================
-- Routine to return the matrix for a gap through material i
--
-- Inputs: S - Snell constant of motion - n*sin(theta)
--         mode - TE or TM mode (s or p)
=========================================================================== */
static M_ARRAY CalcGap(double z, double S, COMPLEX ni, double lambda) {

	M_ARRAY Cij;
	COMPLEX phase;
	COMPLEX cos_theta_i;												/* cosine of angles	*/

	cos_theta_i = CSQRT(1.0-S*S/(ni.x*ni.x+ni.y*ni.y));	/* If < 0 evanescent */
/*	cos_theta_i = sqrt(1.0-pow(S/ni.x,2)); */

/* --------------------------------------------------------------------
 * 2023.07.13 - Mike Thompson
 * Major screwup that has been in the code for 20+ years.  The phase
 * should be ni DIVIDED by cos_theta_i, not MULTIPLIED.
 * ----------------------------------------------------------------- */
	phase = CDIV(ni,cos_theta_i);
	phase.x *= (2*pi/lambda)*z;
	phase.y *= (2*pi/lambda)*z;
/*	phase.x = ni.x*(2*pi/lambda)*z*cos_theta_i; */
/*	phase.y = ni.y*(2*pi/lambda)*z*cos_theta_i; */

	Cij.A.x = cos(phase.x) * exp(-phase.y);
	Cij.A.y = sin(phase.x) * exp(-phase.y);
	Cij.D.x = cos(-phase.x) * exp(phase.y);
	Cij.D.y = sin(-phase.x) * exp(phase.y);
	Cij.B.x = Cij.B.y = Cij.C.x = Cij.C.y = 0.0;
	return Cij;
}


/* ===========================================================================
-- Routine to return the matrix for an interface from i to j
--
-- Inputs: S - Snell constant of motion - n*sin(theta)
--         mode - TE or TM mode (s or p)
=========================================================================== */
static M_ARRAY CalcInterface(double S, POLARIZATION mode, COMPLEX ni, COMPLEX nj) {

	M_ARRAY Cij;
	COMPLEX rij, tij, one={1.0, 0.0};

	CalcFresnel(S, mode, ni, nj, &rij, &tij);

	Cij.A = Cij.D = CDIV(one, tij);
	Cij.B = Cij.C = CDIV(rij, tij);
	return Cij;
}

/* ===========================================================================
-- Routine to calculate the r,t Fresnel coefficient for an interface.
--
-- Inputs: S - Snell constant of motion - n*sin(theta)
--         mode - TE or TM mode (s or p)
=========================================================================== */
static BOOL CalcFresnel(double S, POLARIZATION mode, COMPLEX ni, COMPLEX nj, COMPLEX *rij, COMPLEX *tij) {
	
	COMPLEX a,b;
	COMPLEX cos_theta_i, cos_theta_j;				/* cosine of angles			*/

#if 0																		/* How should k value be handled??? */
	cos_theta_i = CSQRT(1.0-pow(S/ni.x,2));					/* If < 0 evanescent */
	cos_theta_j = CSQRT(1.0-pow(S/nj.x,2));					/* If < 0 evanescent */
#endif
	cos_theta_i = CSQRT(1.0-S*S/(ni.x*ni.x+ni.y*ni.y));	/* If < 0 evanescent */
	cos_theta_j = CSQRT(1.0-S*S/(nj.x*nj.x+nj.y*nj.y));	/* If < 0 evanescent */
#ifdef DEBUG
	printf("S: %g\tmode: %d\tni: %g%+gi\tnj: %g%+gi\n", S, mode, ni.x, ni.y, nj.x, nj.y);
	printf("cos_theta_i: %g%+gi\tcos_theta_j: %g%+gi\n", cos_theta_i.x, cos_theta_i.y, cos_theta_j.x, cos_theta_j.y);
#endif

	if (mode == TE) {													/* Transverse electric */
		a = CSUB(CMUL(ni,cos_theta_i),CMUL(nj,cos_theta_j));
		b = CADD(CMUL(ni,cos_theta_i),CMUL(nj,cos_theta_j));
		*rij = CDIV(a,b);												/* Reflection coefficient */
		a = CMUL(ni,cos_theta_i); a.x *= 2; a.y *= 2;		/* 2*ni*cos_theta_i */
		*tij = CDIV(a,b);
	} else {
		a = CSUB(CMUL(ni,cos_theta_j),CMUL(nj,cos_theta_i));
		b = CADD(CMUL(ni,cos_theta_j),CMUL(nj,cos_theta_i));
		*rij = CDIV(a,b);												/* Reflection coefficient */
		a = CMUL(ni,cos_theta_i); a.x *= 2; a.y *= 2;		/* 2*ni*cos_theta_i */
		*tij = CDIV(a,b);												/* Transmission coefficient */
	}

#ifdef DEBUG
	printf("rij: %g%+gi\ttij: %g%+gi\n", rij->x, rij->y, tij->x, tij->y);
#endif
	return TRUE;
}
		
#if 0												/* Not actually used ... so comment out */
static M_ARRAY MATINV(M_ARRAY *a) {

	COMPLEX det;
	M_ARRAY ai;

	det = CSUB(CMUL(a->A,a->D),CMUL(a->B,a->D));

	ai.A = CDIV(a->D,det);
	ai.B = CDIV(a->B,det); ai.B.x *= -1; ai.B.y *= -1; 
	ai.C = CDIV(a->C,det); ai.C.x *= -1; ai.C.y *= -1; 
	ai.D = CDIV(a->A,det);

	return ai;
}
#endif

static M_ARRAY IDENTITY_MATRIX(void) {

	M_ARRAY ab;

	memset(&ab, 0, sizeof(ab));
	ab.A.x = ab.D.x = 1.0;

	return ab;
}


static M_ARRAY MATMUL(M_ARRAY *a, M_ARRAY *b) {

	M_ARRAY ab;

	ab.A = CADD(CMUL(a->A,b->A),CMUL(a->B,b->C));
	ab.B = CADD(CMUL(a->A,b->B),CMUL(a->B,b->D));
	ab.C = CADD(CMUL(a->C,b->A),CMUL(a->D,b->C));
	ab.D = CADD(CMUL(a->C,b->B),CMUL(a->D,b->D));
	return ab;
}


COMPLEX CADD(COMPLEX a, COMPLEX b) {
	COMPLEX result;
	result.x = a.x+b.x;
	result.y = a.y+b.y;
	return result;
}
	
COMPLEX CSUB(COMPLEX a, COMPLEX b) {
	COMPLEX result;
	result.x = a.x-b.x;
	result.y = a.y-b.y;
	return result;
}
	
COMPLEX CMUL(COMPLEX a, COMPLEX b) {
	COMPLEX result;
	result.x = a.x*b.x - a.y*b.y;
	result.y = a.x*b.y + a.y*b.x;
	return result;
}

COMPLEX CDIV(COMPLEX a, COMPLEX b) {
	COMPLEX result;
	TMPREAL magn;
	magn = b.x*b.x+b.y*b.y;
	result.x =  ( a.x*b.x + a.y*b.y) / magn ;
	result.y =  (-a.x*b.y + a.y*b.x) / magn ;
	return result;
}

double CABS(COMPLEX a) {
	return sqrt(a.x*a.x+a.y*a.y);
}

COMPLEX CSQRT(double r) {
	COMPLEX c;
	c.x = (r>0) ? sqrt(r)  : 0;
	c.y = (r<0) ? sqrt(-r) : 0;
#ifdef DEBUG
	printf("CSQRT: %g %g %g\n", r, c.x, c.y);
#endif
	return c;
}

COMPLEX CCSQRT(COMPLEX r) {
	COMPLEX c;
	double r0,theta;
	if (r.y == 0) r.y = 0;
	r0 = sqrt(r.x*r.x+r.y*r.y);
	theta = atan2(r.y, r.x);
	c.x = sqrt(r0)*cos(theta/2);
	c.y = sqrt(r0)*sin(theta/2);
#ifdef DEBUG
	printf("CSQRT: %g %g = %g %g\n", r.x, r.y, c.x, c.y);
#endif
	return c;
}

COMPLEX CPOW(COMPLEX r, double n) {
	COMPLEX c;
	double r0,theta;
	if (r.y == 0) r.y = 0;
	r0 = sqrt(r.x*r.x+r.y*r.y);
	theta = atan2(r.y, r.x);
	c.x = pow(r0,n)*cos(theta*n);
	c.y = pow(r0,n)*sin(theta*n);
#ifdef DEBUG
	printf("CPOW: (%g %g)^%f = %g %g   %g %g\n", r.x, r.y, n, c.x, c.y, r0, theta);
#endif
	return c;
}
