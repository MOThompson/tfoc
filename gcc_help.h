#ifndef __gcc_help
   #define __gcc_help

#ifdef __GNUC__ 
		int	strcpy_s(char *dest, size_t size, const char *src);
		int	strncpy_s(char *dest, size_t size, const char *src, size_t count);
		int	strcat_s(char *dest, size_t size, const char *src);
		int	_stricmp(const char *str1, const char *str2);
		int	_strnicmp(const char *str1, const char *str2, size_t count);
		int	getenv_s(size_t *requiredSize, char *buffer, size_t bufsize, const char *varname);
		int	fopen_s(FILE **pFile, const char *filename, const char *mode);
#endif
					
#endif	/* __gcc_help */
