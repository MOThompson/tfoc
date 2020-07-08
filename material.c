/* material.c - database handling for TFOC */

/* ------------------------------ */
/* Feature test macros            */
/* ------------------------------ */

/* ------------------------------ */
/* Standard include files         */
/* ------------------------------ */
#include <stdio.h>
#include <string.h>
#include <math.h>
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

#ifndef PATH_MAX
	#define	PATH_MAX	(260)
#endif

struct _TFOC_MATERIAL_MIX {
	EMA_MODEL EMA_Model;							/* What model for effective medium approximation */
	double fraction[MAX_MIX_TERMS];
	TFOC_MATERIAL *material[MAX_MIX_TERMS];
};

struct _TFOC_MATERIAL {
	char name[MATERIAL_NAME_LENGTH];			/* Material name		*/
	void *n_spline, *k_spline;					/* Spline structures	*/
	TFOC_MATERIAL_MIX *mixed;					/* Non-null ==> mixed phase w/ effective medium */
};

/* ------------------------------- */
/* My external function prototypes */
/* ------------------------------- */
TFOC_MATERIAL *TFOC_FindMaterial(char *name, char *database);
COMPLEX TFOC_FindNK(TFOC_MATERIAL *material, double lambda);
void TFOC_PrintMaterials(void);

/* ------------------------------- */
/* My internal function prototypes */
/* ------------------------------- */
static TFOC_MATERIAL *materials=NULL;
static int num_materials=0;
static int dim_materials=0;

/* ------------------------------- */
/* My usage of other external fncs */
/* ------------------------------- */

/* ------------------------------- */
/* Locally defined global vars     */
/* ------------------------------- */


/* ===========================================================================
-- Routine to scan a string and return the material name.  A material
-- may be defined as either a single continuous non-space string, or any 
-- text between matching [] indicating a complex mixture.  This routine
-- is called by sample.c and below in recursion for interpreting a
-- material.
--
-- Usage:  int TFOC_GetMaterialName(char *str, char *name, size_t namelen, char **endptr);
--
-- Inputs: str     - pointer to the string containing the possible name
--         name    - pointer to space to receive the name
--         namelen - number of characters in name (maximum possible length)
--         endptr  - if not NULL, receives pointer to new value of str past name
--
-- Output: *str    - filled with the name of the material
--         *endptr - filled with pointer to next non-space character in str after name
--
-- Return: 0 - simple name
--         1 - complex name (mixture)
=========================================================================== */
int TFOC_GetMaterialName(char *str, char *name, size_t namelen, char **endptr) {

	int ibrace;
	int rc;

	while (isspace(*str)) str++;						/* Skip leading whitespace */

/* Switch depending on the first character */
	if (*str != '[') {									/* Simple name */
		rc = 0;
		while (*str && ! isspace(*str)) {
			*name = *(str++);
			if (namelen > 1) { name++; namelen--; }
		}
	} else {													/* Compound name, go to matching ] */
		rc = 1;
		*name = *(str++);
		if (namelen > 1) { name++; namelen--; }
		ibrace = +1;										/* Have one brace */
		while (ibrace && *str) {
			if (*str == ']') ibrace--;
			if (*str == '[') ibrace++;
			*name = *(str++);
			if (namelen > 1) { name++; namelen--; }
		}
	}
	*name = '\0';

/* Now, skip space and return endptr if requested */
	if (endptr != NULL) {
		while (isspace(*str)) str++;
		*endptr = str;
	}

/* Return indicating simple name or complex name */
	return rc;
}


/* ===========================================================================
-- Routine to locate a material within the existing list of materials,
-- or try to load it from the database directory.
--
-- Usage: TFOC_MATERIAL *TFOC_FindMaterial(char *name, char *database);
--
-- Inputs: name - name of the material - must resolve to a valid filename
--         database - name of a database directory.  File will be searched
--                    for in this directory.
--
-- Output: Adds an entry to the material database list, allocating space
--         and filling in all required information.
--
-- Return: Pointer to the material database structure for specified name.
--         On error, return is NULL.
=========================================================================== */
#define	NPT_MAX	16384
TFOC_MATERIAL *TFOC_FindMaterial(char *name, char *database) {

	int i, rc;
	double ev[NPT_MAX], n[NPT_MAX], k[NPT_MAX];
	int npt=0;
	char filename[PATH_MAX],line[256],*aptr;
	FILE *funit;
	TFOC_MATERIAL *now;

/* See if it already exists */
	while (isspace(*name)) name++;
	for (i=0; i<num_materials; i++) {
		if (_stricmp(name, materials[i].name) == 0) return(materials+i);
	}

/* ----- Make sure there will be space in the sample structure for this entry */
	if (num_materials >= dim_materials) {
		dim_materials += 10;
		materials = realloc(materials, dim_materials*sizeof(*materials));
	}
	now = materials + num_materials;
	num_materials++;

	memset(now, 0, sizeof(*materials));					/* Zero out the entry */
	strncpy_s(now->name, sizeof(now->name), name, sizeof(now->name));

/* ---------------------------------------------------------------------------
-- Okay, not there.  Try to add instead.  Two options.  First is it is
-- a simple name that we look for in the "database" subdirectory.
-- Second is that it begins with a [ and indicates a mixed layer which
-- will be calculated via an effective medium approximation
--------------------------------------------------------------------------- */
	if (*name == '[') {					/* okay - manage as a mixture */
		int rc,i,imix;
		double fsum;
		char myname[MATERIAL_NAME_LENGTH];
		char subname[MATERIAL_NAME_LENGTH];

		strcpy_s(myname, sizeof(myname), name+1);
		aptr = myname+strlen(myname)-1;
		if (*aptr == ']') *(aptr--) = '\0';
		while (isspace(*aptr) && aptr >= myname) *(aptr--) = '\0';

		now->mixed = calloc(1, sizeof(TFOC_MATERIAL_MIX));
		now->mixed->EMA_Model = SERIES;
		imix = 0;							/* Working on first entry */
		rc = 0;
		aptr = myname;
		while (isspace(*aptr)) aptr++;
		/* First look for a model type specifier */
		if (_strnicmp(aptr, "SERIES", 6) == 0) {
			now->mixed->EMA_Model = SERIES;
			aptr += 6; while (isspace(*aptr)) aptr++;
		} else if (_strnicmp(aptr, "PARALLEL", 8) == 0) {
			now->mixed->EMA_Model = PARALLEL;
			aptr += 8; while (isspace(*aptr)) aptr++;
		} else if (_strnicmp(aptr, "BRUGGEMAN", 9) == 0) {
			now->mixed->EMA_Model = BRUGGEMAN;
			aptr += 9; while (isspace(*aptr)) aptr++;
		} else if (_strnicmp(aptr, "LOOYENGA", 8) == 0) {
			now->mixed->EMA_Model = LOOYENGA;
			aptr += 8; while (isspace(*aptr)) aptr++;
		} else if (_strnicmp(aptr, "MAXWELL-GARNETT", 15) == 0) {
			now->mixed->EMA_Model = MAXWELL_GARNETT;
			aptr += 15; while (isspace(*aptr)) aptr++;
		} else if (isalpha(*aptr)) {
			fprintf(stderr, "ERROR: An unrecognized Effective Medium Approximation model was specified (%s)\n", aptr);
			return NULL;
		}

		while (*aptr && (imix < MAX_MIX_TERMS)) {
			now->mixed->fraction[imix] = strtod(aptr, &aptr);
			TFOC_GetMaterialName(aptr, subname, sizeof(subname), &aptr);
			now->mixed->material[imix] = TFOC_FindMaterial(subname, database);
			if (now->mixed->material[imix] == NULL) { rc = 1; break; }
			imix++;
		}
		for (fsum=0,i=0; i<imix; i++) fsum += now->mixed->fraction[i];
		if (rc != 0 || imix == 0 || fsum == 0) {	/* Nothing or totally invalid */
			fprintf(stderr, "ERROR: Parsing mixed-phase material: %s\n", now->name);
			free(now->mixed);
			num_materials--;
			return NULL;
		}
		for (i=0; i<imix; i++) now->mixed->fraction[i] /= fsum;
		
/* ---------------------------------------------------------------------------
-- Not a mixture.  So just look for a datafile (either in the "database"
-- subdirectory or the one specified on argument list) with the name of the
-- material.  It must contain ev,n,k values which will be loaded, spline
-- fit, and the spline coefficient list maintained.
--------------------------------------------------------------------------- */
	} else {
		strcpy_s(filename, sizeof(filename), database);				/* Start creating name */
		if (*filename != '\0') {					/* If not blank */
			aptr = filename + strlen(filename)-1;
			if (*aptr != '/' && *aptr != '\\') strcat_s(filename, sizeof(filename), "/");
		}
		strcat_s(filename, sizeof(filename), name);
		if ( (rc = fopen_s(&funit, filename, "r")) != 0) {
			fprintf(stderr, "ERROR: Unable to open \"%s\" for material \"%s\" in the directory \"%s\" (rc=%d)\n", filename, name, database, rc);
			num_materials--;
			return NULL;
		}

/* Read in all of the existing data */
		while (npt < NPT_MAX && fgets(line, sizeof(line), funit) != NULL) {
			while ( (aptr = strchr(line, '\n')) != NULL) *aptr = '\0';
			aptr = line;
			while (isspace(*aptr)) aptr++;
			if (*aptr == '\0' || *aptr == '#') continue;			/* Comment */
			if (strncmp(aptr, "/*", 2) == 0) continue;			/* Another comment */
			ev[npt] = strtod(aptr, &aptr);
			n[npt]  = strtod(aptr, &aptr);
			k[npt]  = fabs(strtod(aptr, &aptr));				/* So stored as positive */
			npt++;
		}

/* Fit the spline and fill in the entry */
		now->n_spline = GVFitSpline(NULL, ev, n, npt, 0);
		now->k_spline = GVFitSpline(NULL, ev, k, npt, 0);

		if (now->n_spline == NULL || now->k_spline == NULL) {
			free(now->n_spline); free(now->k_spline);
			num_materials--;
			return(NULL);
		}
	}
	
	return now;
}


/* ===========================================================================
-- Obtain the base n,k values from the material database
--
-- Usage: COMPLEX FindNK(MATERIAL *material, double lambda)
--
-- Inputs: material - database for a given material
--         lambda   - wavelength in nm
--
-- Output: none
--
-- Return: complex index of refraction as interpolated from the
--         database files
=========================================================================== */
COMPLEX TFOC_FindNK(TFOC_MATERIAL *material, double lambda) {
	COMPLEX n, b, e1, e2, f, ni[MAX_MIX_TERMS];
	static COMPLEX one={1.0,0.0}, two={2.0,0.0}, three={3.0,0.0}, four={4.0,0.0}, eight={8.0,0.0};
	int i, first_term, last_term, num_terms;
	EMA_MODEL model;
	
	if (material->mixed == NULL) {
		n.x =  GVEvalSpline(material->n_spline, 1240.0/lambda);
		n.y = -GVEvalSpline(material->k_spline, 1240.0/lambda);
	} else {
		model = material->mixed->EMA_Model;
		num_terms = 0;
		first_term = last_term = -1;
		for (i=0; i<MAX_MIX_TERMS; i++) {
			if (material->mixed->fraction[i] <= 0) continue;
			num_terms++;
			if (first_term < 0) first_term = i;
			last_term = i;
			ni[i] = TFOC_FindNK(material->mixed->material[i], lambda);
		}
		if (num_terms == 0) {
			fprintf(stderr, "ERROR: No element has any fraction - returning air\n");
			n.x = 1; n.y = 0;
		/* Simple one component written in complex form */
		} else if (num_terms == 1) {							/* Defined as mixture but only one there */
			n = ni[first_term];
		/* Series weighted average */
		} else if (model == SERIES) {							/* Weighted fraction summation */
			n.x = n.y = 0;
			for (i=0; i<MAX_MIX_TERMS; i++) {
				if (material->mixed->fraction[i] <= 0) continue;
				e1 = ni[i];
				e1 = CMUL(e1,e1);									/* Go to dielectric constant */
				n.x += material->mixed->fraction[i]*e1.x;
				n.y += material->mixed->fraction[i]*e1.y;
			}
			n = CCSQRT(n); if (n.x < 0) { n.x = -n.x; n.y = -n.y; }
		/* Parallel weighted average */
		} else if (model == PARALLEL) {
			n.x = n.y = 0;
			for (i=0; i<MAX_MIX_TERMS; i++) {
				if (material->mixed->fraction[i] <= 0) continue;
				f.x  = material->mixed->fraction[i]; f.y = 0;		/* Fraction */
				n = CADD(n, CDIV(f, CMUL(ni[i],ni[i])));
			}
			n = CDIV(one, n);										/* Invert */
			n = CCSQRT(n); if (n.x < 0) { n.x = -n.x; n.y = -n.y; }
		/* Looyenga - H. Looyenga, "Dielectric constants of heterogeneous mixtures", Physica 31, 401-406 (1965) */
		} else if (model == LOOYENGA) {
			n.x = n.y = 0;
			for (i=0; i<MAX_MIX_TERMS; i++) {
				if (material->mixed->fraction[i] <= 0) continue;
				f.x  = material->mixed->fraction[i]; f.y = 0;		/* Fraction */
				n = CADD(n, CMUL(f, CPOW(ni[i],2.0/3.0)));			/* epsilon^(1/3) from n */
			}
			n = CPOW(n, 1.5);													/* sqrt(sum^3) */
			if (n.x < 0) { n.x = -n.x; n.y = -n.y; }
		/* Maxwell-Garnett mixing */
		} else if (num_terms == 2 && model == MAXWELL_GARNETT) {
			if (material->mixed->fraction[first_term] < material->mixed->fraction[last_term]) {
				f.x = material->mixed->fraction[first_term]; f.y = 0;
				e1  = ni[first_term]; 		/* Minority phase (inclusion) */
				e2  = ni[last_term];			/* Majority phase (matrix) */
			} else {
				f.x = material->mixed->fraction[last_term]; f.y = 0;
				e1  = ni[last_term];			/* Minority phase (inclusion) */
				e2  = ni[first_term];		/* Majority phase (matrix) */
			}
			e1 = CMUL(e1,e1);	e2 = CMUL(e2,e2);				/* Go to dielectric constant */
			n = CDIV ( CSUB(CMUL(e1,CADD(one,CMUL(two,f))),CMUL(e2,CSUB(CMUL(two,f),two))),
						  CADD(CMUL(e2,CADD(two,f)),CMUL(e1,CSUB(one,f))) );
			n = CMUL(e2, n);
			n = CCSQRT(n); if (n.x < 0) { n.x = -n.x; n.y = -n.y; }
		/* Bruggeman model */
		} else if (num_terms == 2 && model == BRUGGEMAN) {
			f.x = material->mixed->fraction[first_term]; f.y = 0;
			e1 = ni[first_term]; e1 = CMUL(e1,e1);			/* Work in dielectric constant */
			e2 = ni[last_term];	e2 = CMUL(e2,e2);			/* Work in dielectric constant */
			b = CADD(CMUL(e1,CSUB(CMUL(three,f),one)), CMUL(e2,CSUB(two,CMUL(three,f))));	/* Actually -b in quadratic */
			n = CCSQRT(CADD(CMUL(b,b),CMUL(eight,CMUL(e1,e2))));
			n = CDIV(CADD(n,b),four);
			n = CCSQRT(n); if (n.x < 0) { n.x = -n.x; n.y = -n.y; }
		} else {
			fprintf(stderr, "ERROR: Either two many terms for Effective Medium Approximation model or unimplemented model (%d)\n", model);
			n.x = 1; n.y = 0;
		}
	}

#ifdef DEBUG
	printf("FindNK returning: %f+%fi\n", n.x, n.y);
#endif
	return(n);
}

/* ===========================================================================
-- Print the values at common wavelengths
--
-- Usage: void PrintMaterials(void)
--
-- Inputs: none
--
-- Output: Prints to standard output
--
-- Return: void
=========================================================================== */
void TFOC_PrintMaterials(void) {

	int i;
	COMPLEX n;

	for (i=0; i<num_materials; i++) {
		fprintf(stderr, "%-10s:", materials[i].name);
		n = TFOC_FindNK(materials+i, 1064.0); fprintf(stderr, "  %f %f", n.x, n.y);
		n = TFOC_FindNK(materials+i,  632.0); fprintf(stderr, "  %f %f", n.x, n.y);
		n = TFOC_FindNK(materials+i,  532.0); fprintf(stderr, "  %f %f", n.x, n.y);
		n = TFOC_FindNK(materials+i,  308.0); fprintf(stderr, "  %f %f", n.x, n.y);
		fprintf(stderr, "\n");
	}
	return;
}
