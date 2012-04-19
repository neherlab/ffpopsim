#############################################################################
#
# Licence:	
# Author:	Richard Neher, Boris Shraiman, Fabio Zanini
# Date:		2012/04/19
#
# Description:
# ------------
# Makefile for the PopGenLib library.
#
##==========================================================================
DIRS := doc src
CLEANDIRS = $(DIRS:%=clean-%)

##==========================================================================
# Rules for compiling and generating documentation.
all: $(DIRS)
$(DIRS):
	$(MAKE) -C $@

clean: $(CLEANDIRS)
$(CLEANDIRS): 
	$(MAKE) -C $(@:clean-%=%) clean


.PHONY: $(DIRS)
.PHONY: $(CLEANDIRS)
.PHONY: all clean
#############################################################################
