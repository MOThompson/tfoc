/* GCC code for VC security library functions (and some string) */
/**
-- VC and Windows have transitioned to secure versions of many string functions 
-- to avoid overflow vulnerabilities.  These routines are a portable implementation
-- of these functions that should compile clean in a simple C environment.
**/

/* ------------------------------ */
/* Feature test macros            */
/* ------------------------------ */

/* ------------------------------ */
/* Standard include files         */
/* ------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <fcntl.h>
#include <ctype.h>
#include <errno.h>

/* ------------------------------ */
/* Local include files            */
/* ------------------------------ */
#include "gcc_help.h"

/* ------------------------------ */
/* My external routines           */
/* ------------------------------ */
/* These are from gcc_help.h */
int	strcpy_s(char *dest, size_t size, const char *src);
int	strncpy_s(char *dest, size_t size, const char *src, size_t count);
int	strcat_s(char *dest, size_t size, const char *src);
int	_stricmp(const char *str1, const char *str2);
int	_strnicmp(const char *str1, const char *str2, size_t count);
int	getenv_s(size_t *requiredSize, char *buffer, size_t bufsize, const char *varname);
int	fopen_s(FILE **pFile, const char *filename, const char *mode);
					
/* ===========================================================================
-- Routine to safely copy a string
-- 
-- Usage:  int strcpy_s(char *dest, size_t size, const char *src);
--
-- Inputs: dest - pointer to destination buffer
--         size - number of elements in dest
--         src  - pointer to source string
--
-- Output: *dest - filled with src up to the defined size
--
-- Return: 0 or error code
--         EINVAL - dest or src NULL
--         ERANGE - size is too small to hold string
--         On error, dest[0] is set to 0 (unless dest is NULL)
=========================================================================== */
int strcpy_s(char *dest, size_t len, const char *src) {
	char *aptr;

	/* Check for errors */
	if (dest == NULL) return EINVAL;
	if (src == NULL || len == 0) { dest[0] = '\0'; return EINVAL; }

	/* Copy over to the EOS */
	aptr = dest;
	while (len-- > 0) {
		if ( (*(aptr++) = *(src++)) == '\0') {
			return 0;
		}
	}
	*dest = '\0';
	return ERANGE;
}


/* ===========================================================================
-- Routine to safely copy up to n characters of a string
-- 
-- Usage:  int strncpy_s(char *dest, size_t size, const char *src, size_t count);
--
-- Inputs: dest  - pointer to destination buffer
--         size  - number of elements in dest
--         src   - pointer to source string
--         count - number of characters to copy
--
-- Output: *dest - filled with src up to the defined size
--
-- Return: 0 or error code
--         EINVAL - dest or src NULL
--         ERANGE - size is too small to hold string
--         On error, dest[0] is set to 0 (unless dest is NULL)
=========================================================================== */
int strncpy_s(char *dest, size_t len, const char *src, size_t count) {
	char *aptr;

	/* Check for errors */
	if (dest == NULL) return EINVAL;
	if (src == NULL || len == 0) { dest[0] = '\0'; return EINVAL; }

	/* Copy over to the EOS or count characters */
	aptr = dest;
	while (len-- > 0 && count-- > 0) {
		if ( (*(aptr++) = *(src++)) == '\0') return 0;
	}

	/* Either there is now space for the terminating EOS, or we have a range error */
	if (len == 0) { *dest = '\0'; return ERANGE; }
	*aptr = '\0';
	return 0;
}


/* ===========================================================================
-- Routine to safely concatenate a string
-- 
-- Usage:  int strcat_s(char *dest, size_t size, const char *src);
--
-- Inputs: dest  - pointer to destination buffer
--         size  - number of elements in dest
--         src   - pointer to source string
--
-- Output: *dest - filled with src up to the defined size
--
-- Return: 0 or error code
--         EINVAL - dest or src NULL
--         ERANGE - size is too small to hold string
--         On error, dest[0] is set to 0 (unless dest is NULL)
=========================================================================== */
int strcat_s(char *dest, size_t len, const char *src) {
	char *aptr;

	/* Check for errors */
	if (dest == NULL) return EINVAL;
	if (src == NULL || len == 0) { dest[0] = '\0'; return EINVAL; }

	/* Walk dest forward either len or to the EOS */
	aptr = dest;
	while (*aptr && len > 0) { aptr++; len--; }

	/* Copy over to the EOS or as much as fits */
	while (len-- > 0) {
		if ( (*(aptr++) = *(src++)) == '\0') return 0;
	}
	*dest = '\0';
	return ERANGE;
}


/* ===========================================================================
--   The stricmp function is a case-insensitive version of ANSI standard
--   strcmp.  The stricmp function compares <string1> and <string2> and returns
--   a value indicating the relationship:
--     < 0       <string1> is less than    <string2>
--     = 0       <string1> is identical to <string2>
--     > 0       <string1> is greater than <string2>
--
-- Usage:   int stricmp(char *string1, char *string2);
--
-- Inputs:  string1 - first  string to be compared
--          string2 - second string to be compared
--
-- Returns: a negative value if <string1> is less than <string2>,
--          0 if <string1> is equal to <string2>, or a positive value
--          if <string1> is greater than <string2>.
=========================================================================== */
int _stricmp(const char *str1, const char *str2) {
	
	register int a,b;
	do {
		if ( (a = tolower(*str1++)) != (b = tolower(*str2++)) ) {
			if (a < b) 
				return(-1);
			else
				return(+1);
		}
	} while ( a != '\0');
	return(0);
}

/* ===========================================================================
--   The strnicmp function is a case-insensitive version of ANSI standard
--   strncmp.  The strnicmp function compares, at most, the first <count>
--   characters of <string1> and <string2> and returns:
--     < 0       <string1> is less than    <string2>
--     = 0       <string1> is identical to <string2>
--     > 0       <string1> is greater than <string2>
--
-- Usage:   int strnicmp(char *string1, char *string2, size_t count);
--
-- Inputs:  string1 - first  string to be compared
--          string2 - second string to be compared
--          count   - maximum number of characters to compare
--
-- Returns: a negative value if <string1> is less than <string2>,
--          0 if <string1> is equal to <string2>, or a positive value
--          if <string1> is greater than <string2>.
=========================================================================== */
int _strnicmp(const char *str1, const char *str2, size_t count) {
	
	register int a,b;
	while (count--) {
		if ( (a = tolower(*str1++)) != (b = tolower(*str2++)) ) {
			if (a < b) 
				return(-1);
			else
				return(+1);
		}
		if ( a == '\0' ) break;
	}
	return(0);
}

/* ===========================================================================
-- Routine to safely copy a string
-- 
-- Usage:  int getenv_s(size_t *requiredSize, char *buffer, size_t bufsize, const char *varname);
--
-- Inputs: requiredSize - pointer to receive needed length to store the environment variable
--         buffer - pointer to destination 
--         bufsize - number of elements allowed in buffer
--         varname - environment variable to be accessed
--
-- Output: *requiredSize - filled with the length of varname or 0 if it does not exit
--         *buffer - filled with the environment variable value if bufsize adequate
--
-- Return: 0 or error code
--         1 if the variable does not exist but otherwise valid
--         EINVAL - requiredSize NULL, varname NULL, or buffer NULL and bufsize > 0
--         ERANGE - buffer is not big enough to store the full string
--         If the variable does not exist, returns 1
=========================================================================== */
int getenv_s(size_t *requiredSize, char *buffer, size_t bufsize, const char *varname) {

	char *aptr;

	if (requiredSize == NULL || (buffer == NULL && bufsize > 0) || varname == NULL) return EINVAL;
	
	if ( (aptr = getenv(varname)) == NULL) {					/* Does not exist */
		*requiredSize = 0;
		if (bufsize > 0) *buffer = '\0';
		return 1;
	}
	if ( (*requiredSize = strlen(aptr)+1) > bufsize) {
		if (bufsize > 0) *buffer = '\0';
		return ERANGE;
	}
	strcpy_s(buffer, bufsize, aptr);
	return 0;
}

/* ===========================================================================
-- Routine to safely open a file
-- 
-- Usage:  int fopen_s(FILE **pFile, const char *filename, const char *mode);
--
-- Inputs: pFile    - pointer to a FILE * variable for the stream descriptor
--         filename - file to be opened
--         mode     - mode file should be opened
--
-- Output: *pFile - filled wit the stream descriptor from fopen()
--
-- Return: 0 on success
--         EINVAL - pFile, filename or mode are NULL
--         errno return value from fopen routine
=========================================================================== */
int fopen_s(FILE **pFile, const char *filename, const char *mode) {

	if (pFile == NULL || filename == NULL || mode == NULL) return EINVAL;

	*pFile = fopen(filename, mode);
	return errno;
}
