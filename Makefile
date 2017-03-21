# conex version
VER := 5.40
# conex extension
# activate this, if you want to manipulate cross-sections and 
# multiplicities etc., or you want to write/read particle lists
# during the simulation
CONEX_EXTENSIONS+=-DCONEX_EXTENSIONS  # 
# system type and directory layout
SYSTEM := $(shell uname -s)
ifdef CONEXBINROOT
  CONEX_PREFIX := $(CONEXBINROOT)
else
  CONEX_PREFIX := $(PWD)
endif
ifndef CONEXTABROOT
  CONEXTABROOT := $(CONEX_PREFIX)/tabs
endif
OBJDIR = $(CONEX_PREFIX)/obj/$(SYSTEM)
LIBDIR = $(CONEX_PREFIX)/lib/$(SYSTEM)
BINDIR = $(CONEX_PREFIX)/bin/
CFGDIR = $(CONEX_PREFIX)/cfg
TABDIR = $(CONEXTABROOT)
# general compiling options
MC3D   = 
#MC3D   = -D__MC3D__ 
LOWMEMORY = -D__SORT_FOR_ENERGY__
#LOWMEMORY= -D__SAVEMEMO__  -D__SORT_FOR_ENERGY__ -D__CXDEBUG__
CPPFLAGS = -D__URQMD__ -D__CXSUB__ -D__PRESHOW__  $(LOWMEMORY) $(CONEX_EXTENSIONS) $(MC3D) -fPIC
# compiler selection
COMP := `which gfortran 2>&1 | grep -v which`
COMPIL := "$(shell echo $(COMP))"
#COMPIL = ""
ifeq ($(COMPIL),"")
 FC  = g77  # for g77
 G2C = g2c  # for g77
else
 FC = gfortran  # for gfortran
 G2C = gfortran  # for gfortran
 FFLAGS += -std=legacy -fPIC
endif
ifndef CXX
  CXX = g++
endif
ifndef CC
  CC = gcc
endif
LD = $(CXX)
CCFLAGS+= $(ROOTCFLAGS)

ifeq ("$(SYSTEM)", "Linux")
   TRAPFPE =  src/trapfpe.c 
   SOFLAGS = -fPIC -ggdb3 -Wall -shared -l$(G2C)
   SOEXT = so
else
ifeq ("$(SYSTEM)", "Darwin")
   TRAPFPE =  src/trapfpe_darwin.c 
   SOFLAGS = -fPIC -ggdb3 -Wall -dynamiclib -l$(G2C) $(ROOTLDFLAGS)
   SOEXT = dylib
endif
endif

CXXFLAGS += $(CONEX_EXTENSIONS)
FFLAGS += -w -fno-second-underscore -fno-automatic -Isrc/ -Isrc/urqmd13/
CXXFLAGS += -g -fPIC $(MC3D) -Isrc 
CCFLAGS += -g -fPIC $(MC3D) 

DEBUG = -ggdb3
OPTIMZE = -O0
FFLAGS += $(OPTIMZE) $(DEBUG)
CFLAGS += $(OPTIMZE) $(DEBUG)
CXXFLAGS += $(OPTIMZE) $(DEBUG)

ROOTCXXFLAGS := -I$(ROOTSYS)/include $(shell root-config --cflags)
ROOTLDFLAGS := -Wl,--no-as-needed $(shell root-config --libs)
ifeq ("$(SYSTEM)", "Darwin")
   LDFLAGS += -bind_at_load -dynamic
endif

OBJS := $(OBJDIR)/conex_f.o $(OBJDIR)/preshw.o $(OBJDIR)/utils.o $(OBJDIR)/veto.o
ifdef CONEX_EXTENSIONS
 OBJS += $(OBJDIR)/particle.o $(OBJDIR)/particleDict.o
 OBJS += $(OBJDIR)/conexResampling.o $(OBJDIR)/conexExtensionsF.o $(OBJDIR)/conexExtensions.o 
 OBJS += $(OBJDIR)/Modifier.o
 OBJS += $(LIBDIR)/libresample.a
 CXXFLAGS += -Iresampler
 LIBRESAMPLE_SOURCES := $(wildcard resampler/*.cc) $(wildcard resampler/*.h)
endif

ifeq ("$(SYSTEM)", "Darwin")
	OBJS += $(OBJDIR)/CxRoot.o $(OBJDIR)/ConexDynamicInterface.o
endif

PARAMS := $(patsubst src/%.paramin, $(CFGDIR)/%.param, $(wildcard src/*.paramin))


.PHONY: std epos qII qgsjetII qgsjet sibyll param_files config_auto \
        commons rebuild svnrevision dirs help first_target all

first_target: help
all: epos sibyll qgsjet qgsjetII

epos: commons binary $(LIBDIR)/libCONEXepos.$(SOEXT)
qgsjet: commons binary $(LIBDIR)/libCONEXqgsjet.$(SOEXT)
qgsjetII: commons binary $(LIBDIR)/libCONEXqgsjetII.$(SOEXT)
qII: commons binary $(LIBDIR)/libCONEXqgsjetII.$(SOEXT)
sibyll: commons binary $(LIBDIR)/libCONEXsibyll.$(SOEXT)

crossSection:
	$(MAKE) -f Makefile.crossSection

commons: svnrevision dirs param_files config_auto
param_files: $(PARAMS)
binary: $(BINDIR)/conex2r

config_auto: 
ifdef CONEX_EXTENSIONS
	echo "#ifndef _include_conexconfigauto_h__" > src/conexConfigAuto.h
	echo "#define _include_conexconfigauto_h__" >> src/conexConfigAuto.h
	echo "#define CONEX_EXTENSIONS 1" >> src/conexConfigAuto.h
	echo "#endif" >> src/conexConfigAuto.h
else
	echo "#ifndef _include_conexconfigauto_h__" > src/conexConfigAuto.h
	echo "#define _include_conexconfigauto_h__" >> src/conexConfigAuto.h
	echo "#endif"	 >> src/conexConfigAuto.h
endif

rebuild: clean all

$(CFGDIR)/%.param: src/%.paramin
	@(echo "==[make]==> generating $@ from : $^")
	@cat $^ | sed s%@CONEXTABROOT@%$(TABDIR)%g > $@ ;


$(OBJDIR)/ConexDynamicInterface.o: src/ConexDynamicInterface.cc src/ConexDynamicInterface.h src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS)  -D_CONEX_PREFIX=\"$(CONEX_PREFIX)\" -D_CONEX_SYSTEM=\"$(SYSTEM)\"  src/ConexDynamicInterface.cc -c -o $@

$(OBJDIR)/CxRoot.o: src/CxRoot.cc src/CxRoot.h src/conexConfig.h src/leadingInteractionsData.cc src/leadingInteractionsData.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS)  -D_CONEX_PREFIX=\"$(CONEX_PREFIX)\" -D_CONEX_SYSTEM=\"$(SYSTEM)\" src/CxRoot.cc -c -o $@

$(BINDIR)/conex2r: src/conex.cc $(OBJDIR)/ConexDynamicInterface.o $(OBJDIR)/CxRoot.o $(OBJDIR)/leadingInteractionsData.o $(OBJDIR)/trapfpe.o src/conexConfig.h
	@(echo "==[make]==> compiling CONEX-ROOT interface")
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) $(ROOTLDFLAGS) $(LDFLAGS) -D_CONEX_PREFIX=\"$(CONEX_PREFIX)\" -D_CONEX_SYSTEM=\"$(SYSTEM)\" $(OBJDIR)/ConexDynamicInterface.o $(OBJDIR)/CxRoot.o $(OBJDIR)/leadingInteractionsData.o -l$(G2C) src/conex.cc -o $@

$(OBJDIR)/conex_f.o: src/conex2r.F src/conexConfig.h
	$(FC) $(CPPFLAGS) $(FFLAGS) -c src/conex2r.F -o $(OBJDIR)/conex_f.o

$(OBJDIR)/Modifier.o: src/Modifier.h src/Modifier.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/Modifier.cc -o $(OBJDIR)/Modifier.o

$(OBJDIR)/leadingInteractionsData.o: src/leadingInteractionsData.h src/leadingInteractionsData.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/leadingInteractionsData.cc -o $(OBJDIR)/leadingInteractionsData.o

$(OBJDIR)/conexExtensions.o: src/conexExtensions.h src/conexExtensions.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/conexExtensions.cc -o $(OBJDIR)/conexExtensions.o

$(OBJDIR)/conexExtensionsF.o: src/conexExtensionsF.h src/conexExtensionsF.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/conexExtensionsF.cc -o $(OBJDIR)/conexExtensionsF.o

$(OBJDIR)/conexResampling.o: src/conexResampling.h src/conexResampling.cc src/conexConfig.h
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c src/conexResampling.cc -o $(OBJDIR)/conexResampling.o

$(OBJDIR)/particleDict.cc: src/particle.h src/particleLinkDef.h
	rootcint -f $@ -c $(CXXFLAGS) $(ROOTCXXFLAGS) $^ 

$(OBJDIR)/particleDict.o: $(OBJDIR)/particleDict.cc
	$(CXX) $(CXXFLAGS) -I. $(ROOTCXXFLAGS) -c $^ -o $@

$(OBJDIR)/particle.o: src/particle.cc
	$(CXX) $(CXXFLAGS) $(ROOTCXXFLAGS) -c $^ -o $@

$(OBJDIR)/preshw.o: src/preshw.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/utils.o: src/utils.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/veto.o: src/veto.c
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/trapfpe.o: $(TRAPFPE)
	$(CC) $(CCFLAGS)  -D_GNU_SOURCE -c $^ -o $@

$(OBJDIR)/conex_sub_sibyll.o: src/conex_sub.F  src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for SIBYLL (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__SIBYLL21__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_epos.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for EPOS (this may take a while ...)")
	$(FC) $(FFLAGS)  $(CPPFLAGS) -D__EPOS__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_qgsjetII.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJETII (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJETII__ -c src/conex_sub.F -o $@
$(OBJDIR)/conex_sub_qgsjet.o: src/conex_sub.F src/conexConfig.h
	@(echo "==[make]==> compiling main CONEX routines for QGSJET01 (this may take a while ...)")
	$(FC) $(FFLAGS) $(CPPFLAGS) -D__QGSJET__ -c src/conex_sub.F -o $@


$(LIBDIR)/libresample.a: $(LIBRESAMPLE_SOURCES)
	@(echo "==[make]==> compiling particle resampling routines ...")
	$(MAKE) -C resampler
	ln -fs ../../resampler/libresample.a $@

clean:
	@rm -rf $(CONEX_PREFIX)/obj
	@rm -rf $(CONEX_PREFIX)/lib
	@rm -rf $(BINDIR)
	@rm -rf $(CFGDIR)
	$(MAKE) clean -C resampler

dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(CFGDIR)


svnrevision:
	@if [ -d .svn ] ; then \
	echo "==[svn revision number]==> $(shell svnversion .)"; \
	sed 's+@SVNREVISION@+"$(shell svnversion .)"+' src/SvnRevisionNumber.h.in > src/SvnRevisionNumber.h; \
	else \
	if [ -e src/SvnRevisionNumber.h.in ] ; then \
	echo "==[svn revision number]==> unknown"; \
	sed 's+@SVNREVISION@+"unknown"+' src/SvnRevisionNumber.h.in > src/SvnRevisionNumber.h; \
	fi \
	fi

$(LIBDIR)/libCONEXqgsjetII.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_qgsjetII.o
	@(echo "")
	@(echo "==[make]==> creating QGSJETII shared library")
	$(LD) $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXqgsjet.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_qgsjet.o
	@(echo "")
	@(echo "==[make]==> creating QGSJET01 shared library")
	$(LD)  $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXsibyll.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_sibyll.o
	@(echo "")
	@(echo "==[make]==> creating SIBYLL shared library")
	$(LD) $^ -o $@ $(SOFLAGS)
$(LIBDIR)/libCONEXepos.$(SOEXT): $(OBJS) $(OBJDIR)/conex_sub_epos.o
	@(echo "")
	@(echo "==[make]==> creating EPOS shared library")
	$(LD) $^ -o $@ $(SOFLAGS)

help:
	@echo ""
	@echo "  CONEX installation Makefile "
	@echo " "
	@echo "  Select interaction models you want to compile "
	@echo "  Options are: 'epos', 'sibyll', 'qgsjet', 'qgsjetII', "
	@echo "               'all' selects all of them"
	@echo " "
	@echo "  Edit 'src/conexConfig.h' for special options"
	@echo " "

dist: svnrevision
	@echo "==[preparing files]==> this will take a while ...";
	@rm -fr conex2r$(VER)
	@mkdir conex2r$(VER)
	@cp -r Makefile README macros/ src/ tabs/ resampler/ eas_disection/ conex2r$(VER)/
	@rm -rf conex2r$(VER)/src/*.h.in
	@find conex2r$(VER) -name "CVS" -type d -print | xargs rm -rf
	@find conex2r$(VER) -name ".svn" -type d -print | xargs rm -rf
	@find conex2r$(VER) -name "\#*" -type f -print | xargs rm -f
	@find conex2r$(VER) -name "*~" -type f -print | xargs rm -f
	@find conex2r$(VER) -name "*.o" -type f -print | xargs rm -f 
	@find conex2r$(VER) -name "*.so" -type f -print | xargs rm -f 
#	@rm -fr conex2r$(VER)/tabs/qgsdat-II-03* # too large
	@rm -fr conex2r$(VER)/tabs/*.*flu* # not needed publicly
	@rm -fr conex2r$(VER)/tabs/*.p2* # not needed publicly
	@echo "==[packing files]==> this will take another while ...";
	@tar -zcf conex2r$(VER).tgz conex2r$(VER)
	@echo "==[finished]==>";
