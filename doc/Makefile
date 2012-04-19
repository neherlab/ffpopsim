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
# Documentation generator
DOXY       = doxygen
DOXYFILE   = Doxyfile

##==========================================================================
# Rules for generating the documentation.
doc:
	$(VERBOSE)$(DOXY) $(DOXYFILE)

doc-%:
	$(VERBOSE)$(DOXY) $(DOXYFILE).$*

##==========================================================================
# Rules for cleaning.
clean:
	rm -rf latex html

#############################################################################
