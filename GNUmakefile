# GNUmakefile, v1.4 21/10/2008e
# GNUmakefile for MEPED module. Karl Yando, 16/07/2006

name := MExED_v1.4
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

# Added 01/22/2007 from N03
visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
