# $Id: GNUmakefile 68058 2013-03-13 14:47:43Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := serc19_cmsmodel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

CPPFLAGS += $(shell root-config --cflags)
EXTRALIBS := $(shell root-config --libs) $(shell root-config --glibs) -lconfig++  -ltbb -lGeom
EXTRALIBS +=-L/usr/lib64 -lQtCore -lQtGui -lQtOpenGL -lXi


include $(G4INSTALL)/config/architecture.gmk

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

