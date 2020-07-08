/* TFOC - Thin film optical calculator.  Mike Thompson.  1999-2020 */
/* ===========================================================================
-- Initial revisions: Verdant summer 1999
-- May 1, 2000 - Added return of transmission, corrected for angle of
--               propagation and index of substrate material.   Can
--               simply ignore the substrate term.
-- Poor revision control history :-)
--
-- Nov. 2002 - Added structure definitions to allow a layer to be
--             automatically expanded into a sequence with varying
--             dopant concentration.  Needed to handle doping profiles
--             for reflectance calculations.
--
-- May. 2003 - Extracted the free-carrier code to separate routine
--
-- June 2004 - Added ability to vary the maximum activation levels on
--             n and p to determine reflectivity as essentially a function
--             of activation temperature.
--
-- March 2007 - Added ability to set maximum n/p activition via command
--              line.  Also new vlog_p0 variation mode.
--
-- August 2008 - Added individual layer temperature capabitility within
--               the free carrier model.  Required to calculate reflectance
--               from the backside for possible OCT analysis.
--
-- July 2020 - Added unpolarized mode to the TM and TE options (average TM,TE) 
=========================================================================== */

/* ------------------------------ */
/* Standard include files         */
/* ------------------------------ */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

/* ------------------------------ */
/* Local include files            */
/* ------------------------------ */
#define TFOC_CODE
#include "tfoc.h"
#include "gcc_help.h"

/* ------------------------------- */
/* My local typedef's and defines  */
/* ------------------------------- */
#ifndef PATH_MAX
	#define	PATH_MAX	(260)
#endif

typedef struct _NKMOD {						/* Enable ability to change n or k of a layer */
	enum {CHANGE_N, CHANGE_K, CHANGE_DOPING, CHANGE_THICKNESS} type;
	int    ival[4];							/* Integer values */
	double rval[4];							/* Real values */
	struct _NKMOD *last,*next;
} NKMOD;

/* ------------------------------- */
/* My external function prototypes */
/* ------------------------------- */

/* ------------------------------- */
/* My internal function prototypes */
/* ------------------------------- */
static void PrintDetails(void);
static void PrintUsage(void);

/* ------------------------------- */
/* My usage of other external fncs */
/* ------------------------------- */

/* ------------------------------- */
/* Locally defined global vars     */
/* ------------------------------- */


/* ===========================================================================
-- Relatively simple routine to update the N,K values in the sample
-- structure from the materials routines.
=========================================================================== */
int main(int argc, char *argv[]) {

	int rc;
	double lambda=632.8;								/* Laser wavelength (nm)	*/
	double theta=0.0;									/* Incident angle (deg)		*/
	double temperature = 300.0;					/* Temperature (K)			*/			
	POLARIZATION mode=TE;							/* Polarization				*/
	TFOC_SAMPLE *sample=NULL;						/* Sample description		*/
	char *samplefilename=NULL;						/* Sample filename			*/
	char *database=NULL,								/* Directory for database	*/
		  *oname=NULL;									/* Output filename			*/
	FILE *funit=NULL;									/* Output filehandle			*/
	BOOL terse=FALSE;									/* Terse output mode?		*/
	BOOL detail=FALSE;								/* Output layer information? */
	NKMOD *tmp, *tmp2;
	REFL result;

	/* For determining the database directory */
	struct _stat info;
	char env_name[PATH_MAX];

	struct {
		enum {NONE, THICKNESS, DUAL, ANGLE, WAVELENGTH, ENERGY, N, K, EXPLOSIVE, FREE_CARRIER, DOPING_PARM_0, DOPING_PARM_1, DOPING_PARM_2,
				DOPING_PARM_0_LOG, TEMP, CPMAX, CNMAX, CMAX} type;
		int layer;
		double min,max,dx;
	} vary={NONE, 0, 0.0,0.0, 1};

/* Initial the list of parameters to change after parsing options */
	NKMOD *PostLoadChanges=NULL;

/* Fresnel calculation layers and number */
	TFOC_LAYER *layers = NULL;						/* Layers for Fresnel calc	*/
	int nlayers = 0;									/* Number of layers			*/

/* Random variables */
	int i,j,npt;										/* Random variables			*/
	size_t cnt;
	double z;
	char *aptr, *endptr;
	BOOL fatal_error;

/* Process the command line arguments - major work */
	fatal_error = FALSE;
	argc--; argv++;
	while (argc>0 && *argv[0] == '-') {
		aptr = *(argv++); argc--;
		aptr++;
		if (_stricmp(aptr, "help") == 0 || *aptr == '?') {
			PrintUsage();
			return(0);

		} else if (_stricmp(aptr, "terse") == 0) {
			terse = TRUE;
			
		} else if (_stricmp(aptr, "manual") == 0) {
			PrintDetails();
			return(0);

		} else if (_stricmp(aptr, "detail") == 0) {
			detail = TRUE;
			
		} else if (_stricmp(aptr, "cmax") == 0) {			/* Set the maximum n/p-type doping */
			if (argc < 1) goto TooFewArgs;
			cnmax = cpmax = fabs(atof(*argv)); argc--; argv++;

		} else if (_stricmp(aptr, "cpmax") == 0) {			/* Set the maximum p-type doping */
			if (argc < 1) goto TooFewArgs;
			cpmax = fabs(atof(*argv)); argc--; argv++;

		} else if (_stricmp(aptr, "cnmax") == 0) {			/* Set the maximum p-type doping */
			if (argc < 1) goto TooFewArgs;
			cnmax = fabs(atof(*argv)); argc--; argv++;

		} else if (_stricmp(aptr, "vn") == 0) {				/* Vary n value of material */
			if (argc < 4) goto TooFewArgs;
			vary.type = N;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;

		} else if (_stricmp(aptr, "vk") == 0) {				/* Vary k value of material */
			if (argc < 4) goto TooFewArgs;
			vary.type = K;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;

		} else if (_stricmp(aptr, "vt") == 0) {				/* Vary thickness		*/
			if (argc < 4) goto TooFewArgs;
			vary.type = THICKNESS;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			vary.max   = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			vary.dx    = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;

		} else if (_stricmp(aptr, "vfe") == 0) {			/* Vary concentration or dose (logarithmic) */
			if (argc < 4) goto TooFewArgs;
			vary.type = FREE_CARRIER;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;

		} else if (_stricmp(aptr, "vlog_p0") == 0) {		/* Vary concentration or dose (logarithmic) */
			if (argc < 4) goto TooFewArgs;
			vary.type = DOPING_PARM_0_LOG;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;	/* Actually number of steps */

		} else if (_stricmp(aptr, "vp0") == 0) {			/* Vary concentration or dose (linear) */
			if (argc < 4) goto TooFewArgs;
			vary.type = DOPING_PARM_0;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;

		} else if (_stricmp(aptr, "vp1") == 0) {			/* Vary profile parameter 1 */
			if (argc < 4) goto TooFewArgs;
			vary.type = DOPING_PARM_1;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;

		} else if (_stricmp(aptr, "vp2") == 0) {			/* Vary profile parameter 2 */
			if (argc < 4) goto TooFewArgs;
			vary.type = DOPING_PARM_2;
			vary.layer = atoi(*argv);	argc--; argv++;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;
			
		} else if (_stricmp(aptr, "vtemp") == 0) {
			if (argc < 3) goto TooFewArgs;
			vary.type  = TEMP;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;	/* Actually number of steps */
			
		} else if (_stricmp(aptr, "vcpmax") == 0) {		/* Vary maximum P activation level */
			if (argc < 3) goto TooFewArgs;
			vary.type = CPMAX;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;	/* Actually number of steps */

		} else if (_stricmp(aptr, "vcnmax") == 0) {		/* Vary maximum N activation level */
			if (argc < 3) goto TooFewArgs;
			vary.type = CNMAX;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;	/* Actually number of steps */

		} else if (_stricmp(aptr, "vcmax") == 0) {			/* Vary maximum activation level */
			if (argc < 3) goto TooFewArgs;
			vary.type = CMAX;
			vary.min   = atof(*argv);	argc--; argv++;
			vary.max   = atof(*argv);	argc--; argv++;
			vary.dx    = atof(*argv);	argc--; argv++;	/* Actually number of steps */
			
		} else if (_stricmp(aptr, "vd") == 0 || _stricmp(aptr, "vm") == 0) {			/* Vary for melt, sum of two layers constant */
			if (argc < 3) goto TooFewArgs;
			vary.type = DUAL;
			vary.layer = atoi(*argv);	argc--;	argv++;
			vary.min   = 0;
			vary.max   = atof(*argv);	argc--;	argv++;
			vary.dx    = atof(*argv);	argc--;	argv++;

		} else if (_stricmp(aptr, "ex" ) == 0) {				/*for explosive crystllization */
			if (argc < 4) goto TooFewArgs;
			vary.type = EXPLOSIVE;
			vary.layer = atoi(*argv);	argc--;	argv++;		/* melt layer, need +1 & -1 to be a/c respectively */
			vary.min   = atoi(*argv);	argc--;	argv++;		/* ?? melt thickness */
			vary.max   = atof(*argv);	argc--;	argv++;		/* Total travel distance = a+c*/
			vary.dx    = atof(*argv);	argc--;	argv++;		/* step size */
			
		} else if (_stricmp(aptr, "va") == 0) {				/* Vary angle			*/
			if (argc < 3) goto TooFewArgs;
			vary.type = ANGLE;
			vary.min   = atof(*argv);	argc--; argv++; 
			vary.max   = atof(*argv);	argc--; argv++; 
			vary.dx    = atof(*argv);	argc--; argv++; 

		} else if (_stricmp(aptr, "vw") == 0) {				/* Vary wavelength	*/
			if (argc < 3) goto TooFewArgs;
			vary.type = WAVELENGTH;
			vary.min   = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			vary.max   = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			vary.dx    = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;

		} else if (_stricmp(aptr, "ve") == 0) {				/* Vary energy	*/
			if (argc < 3) goto TooFewArgs;
			vary.type = ENERGY;
			vary.min   = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			vary.max   = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			vary.dx    = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;

		} else if (_stricmp(aptr, "w") == 0 || _stricmp(aptr, "wavelength") == 0 || _stricmp(aptr, "lambda") == 0) {
			if (argc < 1) goto TooFewArgs;
			lambda = get_nm_value(*argv, &endptr, 0.0);	
			if (*endptr != '\0') goto TrailingGarbage;
			argc--; argv++; 

		} else if (_stricmp(aptr, "TM") == 0) {
			mode = TM;
		} else if (_stricmp(aptr, "TE") == 0) {
			mode = TE;
		} else if (_stricmp(aptr, "random") == 0 || _strnicmp(aptr, "unpolarized", 5) == 0) {
			mode = UNPOLARIZED;
		} else if (_strnicmp(aptr, "polarization", 5) == 0) {
			if (argc < 1) goto TooFewArgs;
			if (_stricmp(*argv, "s") == 0 || _stricmp(*argv, "TE") == 0) {
				mode = TE;
			} else if (_stricmp(*argv, "p") == 0 || _stricmp(*argv, "TM") == 0) {
				mode = TM;
			} else if (_strnicmp(*argv, "unpolarized",3) == 0 || _stricmp(*argv, "random") == 0) {
				mode = UNPOLARIZED;
			} else {
				fprintf(stderr, "Specified polarization (%s) is not recognized.  Use -h for help\n", *argv);
				fatal_error = TRUE;
			}
			argc--; argv++; 

		} else if (_stricmp(aptr, "a") == 0 || _stricmp(aptr, "angle") == 0 || _stricmp(aptr, "AOI") == 0) {
			if (argc < 1) goto TooFewArgs;
			theta = atof(*argv);
			argc--; argv++; 

		} else if (_stricmp(aptr, "temp") == 0 || _stricmp(aptr, "temperature") == 0) {
			if (argc < 1) goto TooFewArgs;
			temperature = atof(*argv);
			argc--; argv++; 

		} else if (_stricmp(aptr, "o") == 0 || _stricmp(aptr, "output") == 0) {
			if (argc < 1) goto TooFewArgs;
			oname = *argv;
			argc--; argv++; 
		} else if (*aptr == 'o') {
			aptr = (aptr[1]!='\0')?(aptr+1):(argc--,*(argv++));
			oname = aptr;

		} else if (_stricmp(aptr, "d") == 0 || _stricmp(aptr, "database") == 0) {
			if (argc < 1) goto TooFewArgs;
			database = *argv;
			argc--; argv++; 

		} else if (_stricmp(aptr, "s") == 0 || _stricmp(aptr, "sample") == 0) {
			if (argc < 1) goto TooFewArgs;
			aptr = *argv;	argc--; argv++; 
			samplefilename = aptr;
			if ( (sample = TFOC_LoadSample(samplefilename)) == NULL) return(-1);

		} else if (_stricmp(aptr, "t") == 0 || _stricmp(aptr, "thickness") == 0) {
			if (argc < 2) goto TooFewArgs;
			tmp = calloc(sizeof(*tmp), 1);
			tmp->type    = CHANGE_THICKNESS;
			tmp->ival[0] = atoi(*argv); argc--; argv++;			/* Layer number, then thickness */
			tmp->rval[0] = get_nm_value(*argv, &endptr, 0.0);	if (*endptr != '\0') goto TrailingGarbage;	argc--; argv++;
			tmp->last = PostLoadChanges;
			PostLoadChanges = tmp;

		} else if (_stricmp(aptr, "dop") == 0 || _stricmp(aptr, "doping") == 0) {
			if (argc < 2) goto TooFewArgs;
			tmp = calloc(sizeof(*tmp), 1);
			tmp->type    = CHANGE_DOPING;
			tmp->ival[0] = atoi(*argv); argc--; argv++;
			tmp->rval[0] = atof(*argv); argc--; argv++;
			tmp->last = PostLoadChanges;
			PostLoadChanges = tmp;

		} else if (_stricmp(aptr, "n") == 0) {				/* Change n of a layer (after sample load) */
			if (argc < 2) goto TooFewArgs;
			tmp = calloc(sizeof(*tmp), 1);
			tmp->type    = CHANGE_N;
			tmp->ival[0] = atoi(*argv); argc--; argv++;
			tmp->rval[0] = atof(*argv); argc--; argv++;
			tmp->last = PostLoadChanges;
			PostLoadChanges = tmp;

		} else if (_stricmp(aptr, "k") == 0) {				/* Change k of a layer (after sample load) */
			if (argc < 2) goto TooFewArgs;
			tmp = calloc(sizeof(*tmp), 1);
			tmp->type    = CHANGE_K;
			tmp->ival[0] = atoi(*argv); argc--; argv++;
			tmp->rval[0] = atof(*argv); argc--; argv++;
			tmp->last = PostLoadChanges;
			PostLoadChanges = tmp;
			
		} else {
			fprintf(stderr, "Switch unrecognized (-%s).  Try -h for help\n", aptr);
			fatal_error = TRUE;
		}
	}

	/* Abort out on fatal errors (after looking for all) */
	if (fatal_error) return 3;

/* If no sample descriptor file yet, take next command line arg or default */
	if (sample == NULL) {
		if (argc > 0) {
			samplefilename = *argv;	argc--; argv++;
		} else {
			samplefilename = "test.sam";
		}
		if ( (sample = TFOC_LoadSample(samplefilename)) == NULL) return(-1);
	}

/* Determine the appropriate database directory (if not set in options) */
	if (database == NULL) {
		if (getenv_s(&cnt, env_name, sizeof(env_name), "tfocDatabase") == 0 && cnt > 0) {	/* Ignore if required size would be more than PATH_MAX */
			database = env_name;
		} else if ( _stat("./tfocDatabase", &info) == 0 && info.st_mode & S_IFDIR ) {
			database = "./tfocDatabase/";
		} else if ( _stat("c:/tfocDatabase", &info) == 0 && info.st_mode & S_IFDIR ) {
			database = "c:/tfocDatabase/";
		} else if ( _stat("c:/database.nk", &info) == 0 && info.st_mode & S_IFDIR ) {		/* Compatibility with earlier versions */
			database = "c:/database.nk/";
		} else {																									/* Better hope materials are in the same directory */
			database = "./";
		}
	}

/* --------------------------------------------------------------------------------
-- Okay, look up the materials and fill in n,k values for each layer directly from
-- the database.  Although we have temperature data, no ability to modify the n,k
-- values based on this information currently.  Only used with the free carrier
-- modification for IR absorption.
-------------------------------------------------------------------------------- */
	for (i=0; sample[i].type != EOS; i++) {
		if ( (sample[i].material = TFOC_FindMaterial(sample[i].name, database)) == NULL) {
			fprintf(stderr, "ERROR: Unable to locate \"%s\" in the materials database directory\n", sample[i].name);
			return(-1);
		}
		sample[i].n = TFOC_FindNK(sample[i].material, lambda);
	}

/* --------------------------------------------------------------------------------
-- At this point, we've looked up the database N,K -- now possibly modify the values 
-- based on requests given on the command line.
-------------------------------------------------------------------------------- */
	if (PostLoadChanges != NULL) {
		for (tmp=PostLoadChanges; tmp->last!=NULL; tmp=tmp->last) tmp->last->next = tmp;
		while (tmp != NULL) {
			switch (tmp->type) {
				case CHANGE_N:
					i = tmp->ival[0];										/* Layer number */
					sample[i].n.x = tmp->rval[0];						/* New n value  */
					break;
				case CHANGE_K:
					i = tmp->ival[0];										/* Layer number */
					sample[i].n.y = -tmp->rval[0];					/* New k value  */
					break;
				case CHANGE_THICKNESS:
					i = tmp->ival[0];										/* Layer number */
					sample[i].z = tmp->rval[0];						/* New thickness */
					break;
				case CHANGE_DOPING:
					i = tmp->ival[0];										/* Layer number */
					if (sample[i].doping_profile == NO_DOPING) sample[i].doping_profile = CONSTANT;
					sample[i].doping_parms[0] = tmp->rval[0];		/* New doping   */
					break;
			}
			tmp2 = tmp; tmp = tmp->next;
			free(tmp2);
		}
		free(tmp);
	}

/* ----------------------------------------------------------
-- Finally, figure out how big the actual layer array will
-- need to be given expansion of profiles, etc.
---------------------------------------------------------- */
	nlayers = 0;
	for (i=0; sample[i].type != EOS; i++) {
		switch (sample[i].doping_profile) {
			case NO_DOPING:
			case CONSTANT:
				nlayers++;
				break;
			case EXPONENTIAL:
			case LINEAR_IMPLANT:
			case LINEAR:
				nlayers += sample[i].doping_layers;
		}
	}
	layers = calloc(nlayers+1, sizeof(*layers));				/* Allocate space (one extra for safety) for layers */

/*	PrintMaterials(); */

/* Open the output file, or just use stdout */
	if (oname != NULL && (rc = fopen_s(&funit, oname, "w")) != 0) {
		fprintf(stderr, "ERROR: Failed to open file \"%s\" for writing (rc=%d)\n", oname, rc);
		return 3;
	} else {
		funit = stdout;
	}

/* And go! */
	if (vary.type == NONE) {
		TFOC_MakeLayers(sample, layers, temperature, lambda);
		result = TFOC_ReflN(theta, mode, lambda, layers);
		fprintf(funit, "%f %f %f\n", result.R, result.T, 1.0-result.R-result.T);
		if (detail) TFOC_PrintDetail(sample, layers);
	} else {

		if (vary.type == EXPLOSIVE) {						/* Some corrections to this mode */
			if ( (vary.dx = fabs(vary.dx)) == 0) vary.dx = 1;
			sample[vary.layer-1].z = -vary.dx;			/* Will be incremented in loop below */
			sample[vary.layer].z=vary.min;
			sample[vary.layer+1].z = vary.max+vary.dx;
			npt = (int) (vary.max/vary.dx + 1.5);
		} else if (vary.type == CPMAX || vary.type == CNMAX || vary.type == CMAX || vary.type == DOPING_PARM_0_LOG) {
			npt = (int) vary.dx;								/* dx is really number of steps */
			if (npt < 2) npt = 2;							/* Must have at least two steps */
			if (npt > 1048576) npt = 1048576;			/* Limit to 2^20 (1 million) calculations */
		} else {
			vary.dx = fabs(vary.dx);
			if (vary.min > vary.max)  vary.dx = -vary.dx;
			if (vary.dx == 0) vary.dx = (vary.max-vary.min)/200.0;
			if (vary.dx == 0) vary.dx = 1.0;							/* So not to blowup */
			npt = (int) ((vary.max-vary.min)/vary.dx + 1.5);
		}

		if (! terse) {													/* Structure info		*/
			fprintf(funit, "# Thin-film optical calculator [v 2.1]\n");
			fprintf(funit, "# ----------------------------------------------------------------------------\n");
			fprintf(funit, "#  Wavelength:     %.2f\n", lambda);
			fprintf(funit, "#  Polarization:   %s\n", (mode==TM)?"TM (p)":(mode==TE)?"TE (s)":"Unpolarized");
			fprintf(funit, "#  Incident Angle: %.2f\n", theta);
			fprintf(funit, "#  Temperature:    %.2f\n", temperature);
			fprintf(funit, "# Sample structure from %s\n", samplefilename);
			for (i=0; sample[i].type != EOS; i++) {
				if (i == 0) {
					fprintf(funit, "#  %2d INCIDENT  %s", i, sample[i].name);
				} else if (sample[i].type == SUBSTRATE) {
					fprintf(funit, "#  %2d SUBSTRATE %s", i, sample[i].name);
				} else {
					fprintf(funit, "#  %2d %8.2f  %s", i, sample[i].z, sample[i].name);
				}
				switch (sample[i].doping_profile) {
					case NO_DOPING:
						break;
					case CONSTANT:
						fprintf(funit, " constant doping=%G", sample[i].doping_parms[0]);
						break;
					case LINEAR:
						fprintf(funit, " linear doping.  Front=%g  Back = %g  nlayers=%d\n", 
								  sample[i].doping_parms[0], sample[i].doping_parms[1], sample[i].doping_layers);
						break;
					case LINEAR_IMPLANT:
						fprintf(funit, " linear implant.  Dose=%g  Front/back heights = %g/%g  nlayers=%d\n", 
								  sample[i].doping_parms[0], sample[i].doping_parms[1], sample[i].doping_parms[2], sample[i].doping_layers);
						break;
					case EXPONENTIAL:
						fprintf(funit, " exponentail doping.  Dose=%g  1/e width=%g  nlayers=%d\n", 
								  sample[i].doping_parms[0], sample[i].doping_parms[1], sample[i].doping_layers);
						break;
					default:
						fprintf(stderr, "ERROR: Unrecognized doping profile in printout section\n");
				}
				fprintf(funit, "\n");
			}
			switch (vary.type) {
				case NONE:
					break;
				case THICKNESS:
					fprintf(funit, "# Thickness of layer %d varied from %f nm to %f nm in %f nm steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case FREE_CARRIER:
					fprintf(funit, "# Free-carrier density of layer %d varied for log(n) = %f to %f in %f steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case DOPING_PARM_0_LOG:
					fprintf(funit, "# Doping parameter 0 (dose/conc) of layer %d varied logarithmically from %g to %g with %d steps\n", vary.layer, vary.min, vary.max, (int) vary.dx);
					break;
				case DOPING_PARM_0:
					fprintf(funit, "# Doping parameter 0 (dose/conc) of layer %d varied from %f to %f in %f steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case DOPING_PARM_1:
					fprintf(funit, "# Doping parameter 1 of layer %d varied from %f to %f in %f steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case DOPING_PARM_2:
					fprintf(funit, "# Doping parameter 2 of layer %d varied from %f to %f in %f steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case CPMAX:
					fprintf(funit, "# Maximum P activation level varied from %.3g to %.3g with %d intervals\n", vary.min, vary.max, (int) vary.dx);
					break;
				case CNMAX:
					fprintf(funit, "# Maximum N activation level varied from %.3g to %.3g with %d intervals\n", vary.min, vary.max, (int) vary.dx);
					break;
				case CMAX:
					fprintf(funit, "# Maximum N/P activation level varied from %.3g to %.3g with %d intervals\n", vary.min, vary.max, (int) vary.dx);
					break;
				case ANGLE:
					fprintf(funit, "# Angle of incidence varied from %f to %f in %f steps\n", vary.min, vary.max, vary.dx);
					break;
				case TEMP:
					fprintf(funit, "# Temperature varied from %f to %f in %f steps\n", vary.min, vary.max, vary.dx);
					break;
				case N:
					fprintf(funit, "# Real part of index of layer %d varied from %f to %f in %f steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case K:
					fprintf(funit, "# Imaginary part of index of layer %d varied from %f to %f in %f steps\n", vary.layer, vary.min, vary.max, vary.dx);
					break;
				case WAVELENGTH:
					fprintf(funit, "# Incident wavelength varied from %f nm to %f nm in %f nm steps\n", vary.min, vary.max, vary.dx);
					break;
				case ENERGY:
					fprintf(funit, "# Incident photon energy varied from %f eV to %f eV in %f eV steps\n", vary.min, vary.max, vary.dx);
					break;
				case DUAL:
					fprintf(funit, "# Dual layer change (melt).  Layer %d varying from 0 to %f nm, %d from %f nm to 0, in steps of %f\n", vary.layer, vary.max, vary.layer+1, vary.max, vary.dx);
					break;
				case EXPLOSIVE:
					fprintf(funit, "# explosive propogation between layer %d and %d, total of %f nm, with melt thickness of %f nm at %f nm steps\n", vary.layer-1, vary.layer+1,vary.max,vary.min,vary.dx);
					break;
			}
			fprintf(funit, "# ----------------------------------------------------------------------------\n");
			fprintf(funit, "# x\tR\tT (into substrate)\n");
		}

		for (i=0; i<npt; i++) {
			z = vary.min + vary.dx*i;
			switch (vary.type) {
				case NONE:
					break;
				case THICKNESS:
					sample[vary.layer].z = z; break;
				case FREE_CARRIER:
					if (sample[vary.layer].doping_profile == NO_DOPING) sample[vary.layer].doping_profile = CONSTANT;
					sample[vary.layer].doping_parms[0] = ((z<0)?-1:+1)*pow(10.0,fabs(z));
					break;
				case DOPING_PARM_0_LOG:
					z = exp(log(fabs(vary.min))+i/(npt-1.0)*(log(fabs(vary.max))-log(fabs(vary.min))));
					if (vary.min < 0) z = -z;
					sample[vary.layer].doping_parms[0] = z;
					break;
				case DOPING_PARM_0:
					sample[vary.layer].doping_parms[0] = z;
					break;
				case DOPING_PARM_1:
					sample[vary.layer].doping_parms[1] = z;
					break;
				case DOPING_PARM_2:
					sample[vary.layer].doping_parms[2] = z;
					break;
				case CPMAX:
				case CNMAX:
				case CMAX:
					z = exp(log(vary.min)+i/(npt-1.0)*(log(vary.max)-log(vary.min)));
					if (vary.type == CPMAX || vary.type == CMAX) cpmax = z; 
					if (vary.type == CNMAX || vary.type == CMAX) cnmax = z; 
					break;
				case EXPLOSIVE:
					sample[vary.layer-1].z += vary.dx;					/* Pre-set so okay at i==0 */
					sample[vary.layer+1].z -= vary.dx;
					z = sample[vary.layer-1].z;							/* Use different coordinate on output */
					break;
				case DUAL:														/* Change two simultaneously, keeping total constant */
					sample[vary.layer].z   = z;
					sample[vary.layer+1].z = vary.max-z;
					break;
				case N:															/* Change the real part only */
					sample[vary.layer].n.x = z;
					break;
				case K:															/* Change the imaginary part only */
					sample[vary.layer].n.y = -z;
					break;
				case ANGLE:
					theta = z;
					break;
				case TEMP:
					temperature = z;
					break;
				case WAVELENGTH:
					lambda = z;
					for (j=0; sample[j].type != EOS; j++) sample[j].n = TFOC_FindNK(sample[j].material, lambda);
					break;
				case ENERGY:
					lambda = (z > 0) ? 1239.842/z : 0.001 ;
					for (j=0; sample[j].type != EOS; j++) sample[j].n = TFOC_FindNK(sample[j].material, lambda);
					break;
			}
			TFOC_MakeLayers(sample, layers, temperature, lambda);
			result = TFOC_ReflN(theta, mode, lambda, layers);
			fprintf(funit, "%g\t%9.7f\t%9.7f\n", z, result.R, result.T);
		}
	}
	if (funit != stdout) fclose(funit);
	return(0);

TooFewArgs:
	fprintf(stderr, "ERROR: Too few arguments for given -%s option", aptr);
	return 3;

TrailingGarbage:
	fprintf(stderr, "ERROR: Trailing garbage on argument (%s) for option -%s", *argv, aptr);
	return 3;
}

/* ===========================================================================
=========================================================================== */
static void PrintUsage(void) {
	printf(
"\n"
"Thin film optical calculator v. 2.3\n"
"   Michael Thompson - mot1@cornell.edu\n"
"\n"
"Usage: tfoc [options] <sample_descriptor_file>\n"
"\n"
"Options:\n"
"     -?                              This help\n"
"     -manual                         More detailed help\n"
"     -detail                         On single calculation, print n,k per layer\n"
"     -a[ngle]       <theta>          Incident angle (in first medium)\n"
"     -w[avelength]  <lambda>[unit>]  Wavelength w/ optional units (nm default)\n"
"     -lambda        <labmda>[<unit>] Wavelength w/ optional units (nm default)\n"
"     -temp[erature] <K>              Temperature in K\n"
"     -polar[ization] [TE | TM | s | p | unp | random]\n"
"     -TM | -TE | -random | -unpol    Polarization (Transverse Magnetic (p), \n"
"                                     Transverse Electric (s), or unpolarized\n"
"\n"
"     -o[utput]      <filename>       Output filename (written)\n"
"     -d[atabase]    <directory>      Specify material n,k database directory\n"
"                                     Checks tfocDatabase environment variable, and\n"
"                                     for ./tfocDatabase or c:/tfocDatabase\n"
"     -s[ample]      <sample file>    Sample filename - processed immediately\n"
"\n"
"     -vt <layer> <min> <max> <dx>    Vary layer thickness (0=incident media)\n"
"     -vd <layer> <range>     <dx>    Vary 2 layers with constant total\n"
"     -vm <layer> <thickness> <dx>    Vary two layer thicknesses to imitate melt\n"
"                                     (identical to -vd)\n"
"     -ex <layer> <liq_t> <tot> <dx>  Simulate explosive crystallization\n"
"                                     Layer n-1 and n+1 vary between 0 and max\n"
"     -va          <min> <max> <dx>   Vary angle of incidence\n"
"     -vw          <min> <max> <dx>   Vary wavelength (nm)\n"
"     -ve          <min> <max> <dx>   Vary energy (eV)\n"
"     -vtemp <min> <max> <dt>         Vary temperature (K)\n"
"     -vn  <layer> <min> <max> <dx>   Vary real part of index of refraction\n"
"     -vk  <layer> <min> <max> <dx>   Vary imaginary part of index of refraction\n"
"     -vfe <layer> <min> <max> <dx>   Increase k by free electron model for Si.\n"
"                                     min/max are log (19=1E19) with + => p-type\n"
"                                     and - => n-type.  k for free-electron is\n"
"                                     added to the table value of k\n"
"     -vp0 <layer> <min> <max> <dx>   Doping parameters 0,1, or 2 are varied in\n"
"     -vp1 <layer> <min> <max> <dx>   a linear fashion.  0 is the concentration\n"
"     -vp2 <layer> <min> <max> <dx>   or dose, 1 is width for exponential, etc.\n"
"     -vlog_p0 <lay> <min> <max> <n>  Vary parameter 0 in logarithmic fashion from\n"
"                                     min to max with n steps.  Sign of <min> and <max>\n"
"                                     must be the same (eg. -1E15 -1E19)\n"
"     -vcpmax <min> <max> <n>         Vary max p activation level in n log steps\n"
"     -vcnmax <min> <max> <n>         Vary max n activation level in n log steps\n"
"     -vcmax  <min> <max> <n>         Vary max n,p activation in n log steps\n"
"\n"
"     -cmax <max>                     Set the maximum activated n & p dopant concentration\n"
"     -cpmax <max>                    Set the maximum activated p-type dopant concentration\n"
"     -cnmax <max>                    Set the maximum activated n-type dopant concentration\n"
"\n"
"     -terse                          Output 2 column data only (no description)\n"
"\n"
"Following options are executed after all others and the sample structure is\n"
"loaded.  They are applied in order of their appearance on the command line.\n"
"     -t[hickness]    <layer> <nm>    Change layer thickness (units on thickness ok)\n"
"     -dop[ing]       <layer> </cc>   Change layer doping\n"
"     -n              <layer> <nval>  Change real part of the index\n"
"     -k              <layer> <kval>  Change imaginary part of the index\n"			 
"\n"
"In the absence of any of the -vx options, the reflectivity and transmission\n"
"of the stack will be returned as a single number.  The -vx options return three\n"
"columns with the parameter value, the reflectivity and the transmission.\n"
"\n"
"Important revision changes:\n"
"   August 1, 2002: The use of combined option/argument is no longer allowed.\n"
"                   Previous use of -w308 must now be -w 308 or -wavelength 308\n"
"   July 6, 2020: Added random polarization (averaged TM and TE values)\n"
	);
	return;
}

/* ===========================================================================
=========================================================================== */
static void PrintDetails(void) {

	PrintUsage();
	printf(
"\n"
"The sample descripor file contains two columns, a material name and the\n"
"thickness, and potentially material options.  The first line is taken as\n"
"the incident medium (the thickness being obviously ignored), and the last\n"
"entry is the substrate medium (again ignoring thickness).  Any layer with\n"
"a zero or negative thickness is left in the structure, but does not\n"
"contribute to the reflectivity.  Numbering within the structure begins\n"
"with 0 for the incident medium, 1 the first surface layer, etc.\n"
"\n"
"The material name may be complex - a mixture of multiple materials with their\n"
"fractions.  These are specified as \"[ {model} 20 SiO2 30 Si3N4 50 c-Si ]\"\n"
"where the fractions need not sum to 100%% (will be normalized) and an optional\n"
"Effective Medium Approximation (EMA) model may be specified.  Anywhere a\n"
"material is required, the mixture form may be specified, including recursively.\n"
"Models must be one of\n"
"    [ SERIES | PARALLEL | BRUGGEMAN | LOOYENGA | MAXWELL-GARNETT ]\n"
"where the series mixing (simple fractional weight of the dielectric constat)\n"
"is default if no model is specified.  SERIES, PARALLEL and LOOYENGA are\n"
"implemented for an arbitrary number of components while the others are only\n"
"valid for two elements.\n"
"\n"
"The only material options recognized relate to doping and temperature.\n"
"Both invoke a free-carrier model for Si, adding the free-carrier k value to\n"
"that returned from the database for the material.  Positive doping values\n"
"correspond to p-type carriers, negative to n-type.  The doping may also be\n"
"specified at the command line using the -doping <layer> <value> options.\n"
"Temperature will override the command line value but is included to permit\n"
"varying temperature profiles to be properly simulated (in the limit).\n"
"\n"
"WARNING: If doping is added to anything other than Si, or specified for\n"
"         wavelengths where it is inappropriate, tough luck.  The user is\n"
"         always responsible for intelligent use of any Thompson code.\n"
"\n"
"The following example is the sample description file for a layered structure\n"
"consisting of a 10 nm Ti/TiN stack on a 30 nm amorphized silicon wafer.  The\n"
"rear surface of the Si is assumed to be non-reflecting for purposes of the\n"
"interference calculations (left as final medium).  The database directory\n"
"would have to contain entries for air, TiN, Ti, a-Si and c-Si.\n"
"The upper layers are at room temperature while the last layer is at 1273 K.\n"
"      AIR\n"
"      TiN   10\n"
"      [ 0.8 Ti 0.2 TiN ] 10\n"
"      a-Si  30\n"
"      c-Si  50 linear 1E20 1E18 20\n"
"      c-Si  50 linear_implant 1E15 1 0 20\n"
"      c-Si  50 exponential -1E15 10 20\n"
"      c-Si  doping=-1E18 temperature=1273\n"
"\n"
"The doping profiles (valid only for Si) may be specified in a number of forms.\n"
"n-type dopants are specified as negative concentrations with p-type positive.\n"
"The type of doping profile is established by one of the keywords.  Each key\n"
"is followed by 1 to 4 parameters depending on the requirements.  Concentrations\n"
"must be specified in /cm^3 values, while implant doses are in /cm^2.\n"
"  * Simple constant concentration profile through film\n"
"      doping <cnc>\n"
"  * Linearly varying concentration from front to back in n layers.\n"
"      linear <front cnc> <back cnc> [<nlayers>]\n"
"  * Implant with linearly varying concentration.  Heights are a.u.\n"
"      linear_implant <dose> <front hgt> <back hgt> [<nlayers>]\n"
"  * Implant with exponential varying concentration.\n"
"      exponential <dose> <1/e width nm> [<nlayers>]\n"
"  * Layer temperature\n"
"      temperature <K>\n"
"For doping profiles, internally the film is just expanded to n layers with the\n"
"integral of the dose over the layer set established as a constant value.  If\n"
"the number of layers is not specified, a reasonable value will be assumed.  All\n"
"sublayers will be at the same temperature - either as specified on the line or\n"
"as the temperature specified on the command line.\n"
"\n"
"The materials database directory must contain a file for each material\n"
"specified in the sample descriptor.  The file contains three columns giving\n"
"photon energy (in eV), real part of the index (n) and the imaginary part\n"
"of the index (k).  The sign of the imaginary part will be ignored - internal\n"
"convention is n = n-ik.  Given entries in the database will be interpolated for\n"
"all wavelengths using a cubic spline with constant extension at the limits.\n"
"\n"
"Lines beginning with the ! character set simulation parameters.  Currently only\n"
"the maximum activated concentrations for n and p can also be set.\n"
"   ! cmax  = <val>\n"
"   ! cpmax = <val>\n"
"   ! cnmax = <val>\n"
"sets the maximum activated concentrations for both, p and n respectively.\n"
"The default values are 1E20 for p and 3E20 for n type doping.  These values\n"
"may also be overwritten by command line options\n"
"   -cmax <val>\n"
"   -cpmax <val>\n"
"   -cnmax <val>\n"
"The sign of the values does not matter, only the absolute value is used.\n"
          );
	return;
}
