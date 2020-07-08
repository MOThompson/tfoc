# Makefile for laser melting simulation program
# #######################
# IBM PC WIN-NT settings
# #######################

CL =
CC = cl
CFLAGS = /nologo /W3

################################################################
TARGET   = tfoc.exe
LIB_FILE = tfoc.lib
LIB_OBJS = fresnel.obj sample.obj material.obj spline.obj free_carrier.obj
INSTALL_TARGET = ..\tfoc.exe 
################################################################

ALL: $(TARGET) $(LIB_FILE) test_tfoc.exe

INSTALL: $(INSTALL_TARGET)

$(INSTALL_TARGET) : $(TARGET)
	cp $(TARGET) $(INSTALL_TARGET)

CLEAN:
	rm *.obj *.exe *.res *.lib

test_tfoc.exe : test_tfoc.obj $(LIB_FILE)
	$(CC) $(CFLAGS) -Fetest_tfoc.exe test_tfoc.obj $(LIB_FILE)

$(TARGET): tfoc.obj $(LIB_OBJS)
	$(CC) $(CFLAGS) -Fe$(TARGET) tfoc.obj $(LIB_OBJS)

.c.obj:
	$(CC) $(CFLAGS) -c $<

tfoc.obj         : tfoc.h gcc_help.h
fresnel.obj      : tfoc.h gcc_help.h
material.obj     : tfoc.h gcc_help.h
sample.obj       : tfoc.h gcc_help.h
spline.obj       : tfoc.h gcc_help.h
free_carrier.obj : tfoc.h gcc_help.h
tfoc_module.obj  : tfoc.h gcc_help.h

$(LIB_FILE) : tfoc_module.obj $(LIB_OBJS)
	if EXIST $@ del $@
	lib /out:$@ $**
