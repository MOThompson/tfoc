/* ===========================================================================
  TFOC_MODULE - Module for linking with other simulation codes
                Main routine reduced to simple subroutine with fixed
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
#define	DATABASE	(NULL)


/* ===========================================================================
-- Routine to determine the reflectance of the sample structure under
-- specified conditions
--
-- Inputs: SampleFile  - name of datafile containing optical description
--         lambda      - wavelenght of light (nm)
--         theta       - angle of incidence
--         mode        - TE, TM, or UNPOLARIZED mode
--         temperature - temperature of the sample
--
-- Output: none
--
-- Return: reflectivity of the sample (absolute)
=========================================================================== */
REFL TFOC_GetRefl(char *SampleFile, double lambda, double theta, POLARIZATION mode, double temperature) {

	int i;
	size_t cnt;
	REFL result;
	struct _stat info;
	char database[PATH_MAX]="", env_name[PATH_MAX];

/* Fresnel calculation layers and number */
/* Variables that are retained after created first time */
	static TFOC_SAMPLE *sample=NULL;				/* Sample description		*/
	static TFOC_LAYER  *layers = NULL;			/* Layers for Fresnel calc	*/
	static int nlayers = 0;							/* Number of layers			*/

/* On first run, load the sample structure for repeated calculations */
	if (sample == NULL) {
		sample = TFOC_LoadSample(SampleFile);
		if (sample == NULL) {
			fprintf(stderr, "ERROR: Input sample description file (%s) is invalid\n", SampleFile);
			abort();
		}

		/* Determine the appropriate database directory (if not set in options) */
		if (*database == '\0') {
			/* Check the environment variables first ... one called tfocDatabase or in the tfoc directory of LocalAppData */
			if (getenv_s(&cnt, env_name, sizeof(env_name), "tfocDatabase") == 0 && cnt > 0) {	/* Ignore if required size would be more than PATH_MAX */
				strcpy_s(database, sizeof(database), env_name);
			} else if (getenv_s(&cnt, env_name, sizeof(env_name), "LocalAppData") == 0 && cnt > 0) {
				sprintf_s(database, sizeof(database), "%s/TFOC/tfocdatabase", env_name);
				if (_stat(database, &info) != 0 || ! (info.st_mode & S_IFDIR) ) {
					sprintf_s(database, sizeof(database), "%s/TFOC/tfocdatabase.nk", env_name);
					if (_stat(database, &info) != 0 || ! (info.st_mode & S_IFDIR) ) {
						sprintf_s(database, sizeof(database), "%s/TFOC/database", env_name);
						if (_stat(database, &info) != 0 || ! (info.st_mode & S_IFDIR) ) {
							sprintf_s(database, sizeof(database), "%s/TFOC/database.nk", env_name);
							if (_stat(database, &info) != 0 || ! (info.st_mode & S_IFDIR) ) *database = '\0';
						}
					}
				}
			}
		}

		/* If not installed formally, check a number of well known names and paths */
		if (*database == '\0') {
			if ( _stat("./tfocDatabase", &info) == 0 && info.st_mode & S_IFDIR ) {
				strcpy_s(database, sizeof(database), "./tfocDatabase");
			} else if ( _stat("./tfocDatabase.nk", &info) == 0 && info.st_mode & S_IFDIR ) {
				strcpy_s(database, sizeof(database), "./tfocDatabase.nk");
			} else if ( _stat("c:/tfocDatabase", &info) == 0 && info.st_mode & S_IFDIR ) {
				strcpy_s(database, sizeof(database), "c:/tfocDatabase");
			} else if ( _stat("c:/database.nk", &info) == 0 && info.st_mode & S_IFDIR ) {		/* Compatibility with earlier versions */
				strcpy_s(database, sizeof(database), "c:/database.nk");										
			} else {																									/* Better hope materials are in the same directory */
				strcpy_s(database, sizeof(database), ".");
			}
		}
		strcat_s(database, sizeof(database), "/");														/* Append trailing path delimiter */

		/* ----------------------------------------------------------
		   -- Pre-process sample structure - don't have temperature yet
			-- Identify the material database information and
			-- at the same time, figure out how big the actual layer
			-- array will need to be given expansion of profiles, etc.
			---------------------------------------------------------- */
		nlayers = 0;
		for (i=0; sample[i].type != EOS; i++) {
			if ( (sample[i].material = TFOC_FindMaterial(sample[i].name, DATABASE)) == NULL) {
				fprintf(stderr, "ERROR: Unable to locate %s in the materials database directory\n", sample[i].name);
				result.R = result.T = -1;
				return result;
			}
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
	}

/* ----------------------------------------------------------
-- Okay, look up the materials and fill in n,k values
-- At the same time, figure out how big the actual layer
-- array will need to be given expansion of profiles, etc.
---------------------------------------------------------- */

/* And go! */
	for (i=0; sample[i].type != EOS; i++) {
		sample[i].n = TFOC_FindNK(sample[i].material, lambda);
	}
	TFOC_MakeLayers(sample, layers, temperature, lambda);
	result = TFOC_ReflN(theta, mode, lambda, layers);
/*	printf("%f %f %f\n", result.R, result.T, 1.0-result.R-result.T); */

	return result;
}


#ifdef TEST

int main(int argc, char *argv[]) {

	printf("%f", TFOC_GetRefl("a.sam", 10600, 75, TE, 300.0));
	return 0;
}

#endif
