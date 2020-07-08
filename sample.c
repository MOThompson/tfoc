/* sample.c - routines for interpreting sample structure to layer configuration */

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

/* ------------------------------- */
/* My external function prototypes */
/* ------------------------------- */
TFOC_SAMPLE *TFOC_LoadSample(char *fname);
void TFOC_MakeLayers(TFOC_SAMPLE *sample, TFOC_LAYER *layers, double T, double lambda);
void TFOC_PrintDetail(TFOC_SAMPLE *sample, TFOC_LAYER *layers);

/* ------------------------------- */
/* My internal function prototypes */
/* ------------------------------- */

/* ------------------------------- */
/* My usage of other external fncs */
/* ------------------------------- */

/* ------------------------------- */
/* Locally defined global vars     */
/* ------------------------------- */

/* ------------------------------- */
/* My share of global vars         */
/* ------------------------------- */
double cpmax = 1E20, cnmax=3E20;				/* Maximum activated concentrations n and p */

#define	LIMIT_DOPING(n)	(((n)>cpmax)?(cpmax):(((n)<-cnmax)?(-cnmax):(n)))

/* ===========================================================================
-- Convert the high level sample description into a number of fundamental
-- layers to run in Fresnel.  The layers need only n and z.  Expand the
-- doping profiles in this routine.
=========================================================================== */
void TFOC_MakeLayers(TFOC_SAMPLE *sample, TFOC_LAYER *layers, double T, double lambda) {
	double a,b,peak,posn,doping, dz, w, temperature;
	int i;
	TFOC_LAYER *lay;
	TFOC_SAMPLE *sam;
	FC_MODE fc_mode = KLAASSEN_MU;

	int ilay;

	for (ilay=0,sam=sample,lay=layers; sam->type!=EOS; ilay++,sam++) {

		if (sam->type == IGNORE_LAYER) continue;
		temperature = (sam->temperature <= 0) ? T : sam->temperature;

		switch (sam->doping_profile) {
			case NO_DOPING:
				lay->layer  = ilay;									/* Defining layer	*/
				lay->type   = SUBLAYER;								/* Type of layer	*/
				lay->n      = sam->n;								/* Index				*/
				lay->z      = sam->z;								/* And thickness	*/
				lay->doping  = 0;
				lay++;
				break;
			case CONSTANT:
				lay->layer  = ilay;									/* Defining layer	*/
				lay->type   = SUBLAYER;								/* Type of layer	*/
				lay->n      = sam->n;								/* Index				*/
				lay->z      = sam->z;								/* And thickness	*/
				lay->doping = LIMIT_DOPING(sam->doping_parms[0]);				/* Doping level	*/
				lay->n.y += -fc_k(lay->doping, 0, temperature, lambda/1000.0, fc_mode);	/* FE Correction	*/
				lay++;
				break;
			case LINEAR:
				dz   = sam->z / sam->doping_layers;				/* Width of each layer	*/
				a    = sam->doping_parms[0];						/* Front level (/cm^3)	*/
				b    = sam->doping_parms[1];						/* Back  level (/cm^3)	*/
				for (i=0; i<sam->doping_layers; i++) {
					posn = (i+0.5)/sam->doping_layers;			/* Relative position [0,1] */
					doping = a*(1-posn) + b*posn;
					lay->layer  = ilay;
					lay->type   = SUBLAYER;
					lay->n      = sam->n;
					lay->z      = dz;
					lay->doping = LIMIT_DOPING(doping);
					lay->n.y += -fc_k(lay->doping, 0, temperature, lambda/1000.0, fc_mode);
					lay++;
				}
				break;
			case LINEAR_IMPLANT:
				dz   = sam->z / sam->doping_layers;				/* Width of each layer	*/
				a    = sam->doping_parms[1];						/* Front level (a.u.)	*/
				b    = sam->doping_parms[2];						/* Back  level (a.u.)	*/
				peak = 2*sam->doping_parms[0]/(sam->z*1E-7) / (a+b);		/* To cm^3	*/
				a = a*peak; b = b*peak;								/* Front/back conc		*/
				for (i=0; i<sam->doping_layers; i++) {
					posn = (i+0.5)/sam->doping_layers;			/* Relative position [0,1] */
					doping = a*(1-posn) + b*posn;
					lay->layer  = ilay; 
					lay->type   = SUBLAYER;
					lay->n      = sam->n;
					lay->z      = dz;
					lay->doping = LIMIT_DOPING(doping);
					lay->n.y += -fc_k(lay->doping, 0, temperature, lambda/1000.0, fc_mode);
					lay++;
				}
				break;
			case EXPONENTIAL:											/* Properly handles integral of exponential */
				dz   = sam->z / sam->doping_layers;				/* Width of each layer	*/
				w    = sam->doping_parms[1];						/* 1/e width				*/
				peak = (sam->doping_parms[0] * (1-exp(-dz/w)) / (1-exp(-sam->z/w))) / (dz*1E-7);
				for (i=0; i<sam->doping_layers; i++) {
					posn = (i*sam->z)/sam->doping_layers;		/* Start position of layer */
					doping = peak*exp(-posn/w);
					lay->layer  = ilay; 
					lay->type   = SUBLAYER;
					lay->n      = sam->n;
					lay->z      = dz;
					lay->doping = LIMIT_DOPING(doping);
					lay->n.y += -fc_k(lay->doping, 0, temperature, lambda/1000.0, fc_mode);
					lay++;
				}
				break;
			default:
				fprintf(stderr, "ERROR: Unrecognized case (%d) in doping_profile (layer=%d)\n", sam->doping_profile, ilay);
				break;
		}
	}

/* Identify the first and last elements */
	layers[0].type = INCIDENT;											/* First will be the incident medium	*/
	(--lay)->type = EOS;													/* Last will be the exit substrate		*/

#ifdef DEBUG
	for (i=0,lay=layers; ; i++,lay++) {
		printf("%d.%d\tz=%g\tn.x=%.5f\tn.y=%.5f\ttype=%d\tdoping=%g\n", lay->layer, i, lay->z, lay->n.x, lay->n.y, lay->type, lay->doping);
		if (lay->type == EOS) break;
	}
#endif
	return;
}

/* ===========================================================================
-- Routine to just jump over optional blank space and = signs
=========================================================================== */
static char *SkipOptEq(char *aptr) {
	while (isspace(*aptr)) aptr++;
	if (*aptr == '=') aptr++;
	while (isspace(*aptr)) aptr++;
	if (strchr("0123456789.+-", *aptr) == NULL) {
		fprintf(stderr, "ERROR: Reading sample structure.  Expected number but got \"%s\"\n", aptr);
		return NULL;
	}
	return aptr;
}

/* ===========================================================================
-- Routine to load the sample structure from a database file.
--
-- Usage: SAMPLE *LoadSample(char *fname)
--
-- Inputs: fname - filename containing text lines describing the sample
--                 If NULL, then input will be taken from <stdin>
--
-- Output: none
--
-- Return: Returns pointer to allocated structure describing sample.  Calling
--         program is ultimately responsible for releasing the structure.
=========================================================================== */
TFOC_SAMPLE *TFOC_LoadSample(char *fname) {

	int i, rc;

	TFOC_SAMPLE *sample=NULL, *sam;
	int num_layers = 0;
	int dim_layers = 0;

	FILE *funit;
	char line[256], *name, *aptr;

	if (fname == NULL) {
		funit = stdin;
	} else if ( (rc = fopen_s(&funit, fname, "r")) != 0) {
		fprintf(stderr, "ERROR: File \"%s\" specifying sample structure does not exist (rc = %d)\n", fname, rc);
		return NULL;
	}

	while (fgets(line, sizeof(line), funit) != NULL) {
		if ((aptr=strchr(line,'\n'))!=NULL) *aptr='\0';			/* Dump any <newline> char */
		aptr = line; while (isspace(*aptr)) aptr++;				/* Skip whitespace			*/
		if (*aptr == '\0' || *aptr == '#' || *aptr == '%' || strncmp(aptr, "/*", 2) == 0) continue;

		if (*aptr == '!') {												/* Internal parameters		*/
			aptr++;															/* Skip over the ! sign		*/
			while (isspace(*aptr)) aptr++;							/* Find start of name		*/
			name = aptr;
			while (*aptr && ! isspace(*aptr) && *aptr != '=') aptr++;				/* Go to the end of name	*/
			if (*aptr != '\0') *(aptr++) = '\0';					/* EOS, or do we fake?		*/

			while (isspace(*aptr)) aptr++;							/* Skip white space			*/
			while (*aptr == '=') aptr++;								/* Skip any = signs			*/
			while (isspace(*aptr)) aptr++;							/* Skip white space			*/

			if (_strnicmp(name, "CMAX", 4) == 0) {					/* Max possible P/N concentration */
				cpmax = cnmax = strtod(aptr, &aptr);
			} else if (_strnicmp(name, "CPMAX", 5) == 0) {
				cpmax = strtod(aptr, &aptr);
			} else if (_strnicmp(name, "CNMAX", 5) == 0) {
				cnmax = strtod(aptr, &aptr);
			} else {
				fprintf(stderr, "ERROR: Unable to process line starting with %s\n", line);
				continue;
			}
			continue;
		}

		/* Create a slot in the sample structure */
		if (num_layers+1 >= dim_layers) {							/* Create space if needed! */
			dim_layers += 20;
			sample = realloc(sample, dim_layers*sizeof(*sample));
		}
		sam = sample+num_layers;										/* This layer					*/
		memset(sam, 0, sizeof(*sam));									/* Zero out all parameters	*/

		/* Fill in the database for the layer - starting with the material name */
		TFOC_GetMaterialName(aptr, sam->name, sizeof(sam->name), &aptr);

		sam->material = NULL;											/* No database loaded		*/
		sam->type     = (num_layers==0)?INCIDENT:SUBLAYER;

		/* Doping profile information */
		sam->doping_profile  = NO_DOPING;							/* No doping at first		*/
		sam->doping_layers   = 1;										/* No sublayers				*/
		sam->temperature     = -1;										/* Temperature undefined	*/
		for (i=0; i<NPARMS_DOPING; i++) sam->doping_parms[i] = 0;

		/* Next value is the thickness in nm (or other units if specified) */
		sam->z = get_nm_value(aptr, &aptr, 0.0);				/* Get the thickness */

/* Check for complex doping profiles */
		while (*aptr != '\0' && *aptr != '\n' && *aptr != '#' && *aptr != '%' && _strnicmp(aptr, "/*", 2) != 0 && _strnicmp(aptr, "//",2) != 0) {
			
			if (_strnicmp(aptr, "temperature", 11) == 0) {			/* Layer temperature */
				if ( (aptr = SkipOptEq(aptr+11)) == NULL) return NULL;
				sam->temperature = strtod(aptr, &aptr);
			} else if (_strnicmp(aptr, "temp", 4) == 0) {				/* Layer temperature */
				if ( (aptr = SkipOptEq(aptr+4)) == NULL) return NULL;
				sam->temperature = strtod(aptr, &aptr);
			} else if (_strnicmp(aptr, "t", 1) == 0) {					/* Layer temperature */
				if ( (aptr = SkipOptEq(aptr+1)) == NULL) return NULL;
				sam->temperature = strtod(aptr, &aptr);

			} else if (_strnicmp(aptr, "doping",6) == 0) {			/* Simple doping spec		*/
				if ( (aptr = SkipOptEq(aptr+6)) == NULL) return NULL;
				sam->doping_profile  = CONSTANT;
				sam->doping_parms[0] = strtod(aptr, &aptr);			/* Doping level (/cm^3)		*/

			} else if (_strnicmp(aptr, "linear_implant",14)==0) {	/* Linear implant profile	*/
				if ( (aptr = SkipOptEq(aptr+14)) == NULL) return NULL;
				sam->doping_profile  = LINEAR_IMPLANT;
				sam->doping_parms[0] = strtod(aptr, &aptr);			/* Dose (/cm^2)				*/
				sam->doping_parms[1] = strtod(aptr, &aptr);			/* Front level					*/
				sam->doping_parms[2] = strtod(aptr, &aptr);			/* Back level					*/
				sam->doping_layers   = strtol(aptr, &aptr, 10);		/* # sublayers					*/
				if (sam->doping_layers <= 1) sam->doping_layers = 10;

			} else if (_strnicmp(aptr, "linear",6)==0) {				/* Linear profile				*/
				if ( (aptr = SkipOptEq(aptr+6)) == NULL) return NULL;
				sam->doping_profile  = LINEAR;
				sam->doping_parms[0] = strtod(aptr, &aptr);			/* Front concentration		*/
				sam->doping_parms[1] = strtod(aptr, &aptr);			/* Back concentration		*/
				sam->doping_layers   = strtol(aptr, &aptr, 10);		/* # sublayers					*/
				if (sam->doping_layers <= 1) sam->doping_layers = 10;

			} else if (_strnicmp(aptr, "exponential", 11) == 0) {	/* Exponential profile */
				if ( (aptr = SkipOptEq(aptr+11)) == NULL) return NULL;
				sam->doping_profile  = EXPONENTIAL;
				sam->doping_parms[0] = strtod(aptr, &aptr);			/* Dose (/cm^2)				*/
				sam->doping_parms[1] = strtod(aptr, &aptr);			/* 1/e width				 	*/
				sam->doping_layers   = strtol(aptr, &aptr, 10);		/* # sublayers					*/
				if (sam->doping_layers <= 1) sam->doping_layers = (int) (5*sam->z/sam->doping_parms[1]+1);

			} else {
				fprintf(stderr, "ERROR: Unrecognized text following layer definition\n\t\"%s\"\n", aptr);
				if (funit != stdin) fclose(funit);
				if (sample != NULL) { free(sample); sample = NULL; }
				return NULL;
			}

			while (isspace(*aptr)) aptr++;								/* Skip any spaces now */
		}

		num_layers++;
	}

	if (num_layers == 0 || sample == NULL) {
		fprintf(stderr, "ERROR: Empty sample structure\n");
		if (sample != NULL) { free(sample); sample = NULL; }
	} else {
		sample[num_layers-1].type = SUBSTRATE;
		sample[num_layers].type   = EOS;						/* End of sample */
	}

	if (funit != stdin) fclose(funit);
	return sample;
}



/* ===========================================================================
-- Print the layer and n/k values used in the calculation.  Essentially
-- current values in the data structures.
--
-- Usage: void TFOC_PrintDetail(TFOC_SAMPLE *sample, TFOC_LAYER *layers);
--
-- Inputs: sample - sample structure
--         layers - layer informaion
--
-- Output: Prints to standard output
--
-- Return: void
=========================================================================== */
void TFOC_PrintDetail(TFOC_SAMPLE *sample, TFOC_LAYER *layers) {

	int i;

	printf("Layer\tnm\tn\tk\tname\n");
	for (i=0; layers->type != EOS; i++,layers++) {
		printf("%d\t%f\t%f\t%f\t%s\n", i, layers->z, layers->n.x, -layers->n.y, sample[layers->layer].name);
	}
	return;
}


/* ===========================================================================
-- Routine to parse a thickness and optionally a dimensional unit
--
-- Usage: double get_nm_value(char *aptr, char **endptr, double dflt);
--
-- Inputs: aptr   - pointer to text
--         endptr - optional pointer to receive next unused character in aptr
--         dflt   - default return value if text contains no value
-- 
-- Output: *endptr - if endptr != NULL, filled with aptr value after scan
--
-- Return: Thickness converted to nm dimensions
--
-- Notes: Allows optional separating whitespace between value and units
--        Units recognized include A, nm, um, mm and cm
=========================================================================== */
double get_nm_value(char *aptr, char **endptr, double dflt) {
	double rc;

	rc = dflt;
	rc = strtod(aptr, &aptr);										/* As much as relevant		*/
	while (isspace(*aptr)) aptr++;								/* Skip whitespace			*/

/* Check for an optional unit specification */
	if (strncmp(aptr, "A", 1) == 0 && (aptr[1] == '\0' || isspace(aptr[1])) ) {
		rc /= 10;														/* A to nm */
		aptr++;
	} else if (strncmp(aptr, "nm", 2) == 0 && (aptr[2] == '\0' || isspace(aptr[2])) ) {
		rc *= 1;															/* nm to nm */
		aptr += 2;
	} else if (strncmp(aptr, "um", 2) == 0 && (aptr[2] == '\0' || isspace(aptr[2])) ) {
		rc *= 1000;														/* um to nm */
		aptr += 2;
	} else if (strncmp(aptr, "mm", 2) == 0 && (aptr[2] == '\0' || isspace(aptr[2])) ) {
		rc *= 1000000;													/* mm to nm */
		aptr += 2;
	} else if (strncmp(aptr, "cm", 2) == 0 && (aptr[2] == '\0' || isspace(aptr[2])) ) {
		rc *= 10000000;												/* cm to nm */
		aptr += 2;
	}
	while (isspace(*aptr)) aptr++;								/* And skip whitespaces again	*/

	if (endptr != NULL) *endptr = aptr;
	return rc;
}
