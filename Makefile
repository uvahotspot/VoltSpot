# Uncomment the following math acceleration flags 
# relevant to your target and set the appropriate
# path and flag options

#Super LU
SuperLUroot	= $(HOME)/Runjie/Temp/SuperLU_4.3
SUPERLULIB 	= $(SuperLUroot)/lib/libsuperlu_4.3.a
BLASLIB    	= -lblas
SLU_HEADER  = $(SuperLUroot)/SRC

INCDIR		= $(SLU_HEADER)
LIBDIR		= 
LIBS  		= -lm $(SUPERLULIB) $(BLASLIB)
EXTRAFLAGS	= -Wno-unused-result

CC 			= gcc
ifeq ($(DEBUG), 1)
OFLAGS		= -O0 -ggdb -Wall
else
ifeq ($(DEBUG), 2)
OFLAGS		= -O3 -pg -ggdb -Wall
else
OFLAGS		= -O3
endif	# DEBUG = 2
endif	# DEBUG = 1

RM			= rm -f
AR			= ar qcv
RANLIB	= ranlib
OEXT		= o
LEXT		= a

# Verbosity level [0-3]
ifndef VERBOSE
VERBOSE	= 0
endif

ifdef INCDIR
INCDIRFLAG = -I$(INCDIR)
endif

ifdef LIBDIR
LIBDIRFLAG = -L$(LIBDIR)
endif

CFLAGS	= $(OFLAGS) $(EXTRAFLAGS) $(INCDIRFLAG) $(LIBDIRFLAG) -DVERBOSE=$(VERBOSE)

# sources, objects, headers and inputs

# Floorplan
FLPSRC	= flp.c
FLPOBJ	= flp.$(OEXT)
FLPHDR	= flp.h

# PDN model
PDNSRC  = PDN_sim.c PDN_analyze.c pad.c matrix.c
PDNOBJ  = PDN_sim.$(OEXT) PDN_analyze.$(OEXT) pad.$(OEXT) matrix.$(OEXT)
PDNHDR	= PDN_sim.h PDN_analyze.h pad.h

# Miscellaneous
MISCSRC = util.c
MISCOBJ = util.$(OEXT)
MISCHDR = util.h

# all objects
OBJ	= $(PDNOBJ) $(FLPOBJ) $(MISCOBJ)

# targets
all:	voltspot lib

voltspot:	voltspot.$(OEXT) $(OBJ)
	$(CC) $(CFLAGS) -o voltspot voltspot.$(OEXT) $(OBJ) $(LIBS)
ifdef LIBDIR
		@echo
		@echo
		@echo "...Done. Do not forget to include $(LIBDIR) in your LD_LIBRARY_PATH"
endif

lib:	voltspot
	$(RM) libvoltspot.$(LEXT)
	$(AR) libvoltspot.$(LEXT) $(OBJ)
	$(RANLIB) libvoltspot.$(LEXT)

#pull in dependency info for existing .o files
-include $(OBJ:.o=.d)

.c.$(OEXT):
	$(CC) $(CFLAGS) -c $*.c
	$(CC) -MM $(CFLAGS) $*.c > $*.d

filelist:
	@echo $(FLPSRC) $(PDNSRC) $(MISCSRC) \
		  $(FLPHDR) $(PDNHDR) $(MISCHDR) \
		  voltspot.h voltspot.c\
		  Makefile
clean:
	$(RM) *.$(OEXT) *.obj *.d core *~ Makefile.bak voltspot libvoltspot.$(LEXT)
	$(RM) $(OBJ)

cleano:
	$(RM) *.$(OEXT) *.obj

