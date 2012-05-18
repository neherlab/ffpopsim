/**
 * @file hivgene.cpp
 * @brief Implementation of an HIV gene
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-05-18
 */
#include "hivpopulation.h"

/**
 * @brief Default constructor for an HIV gene
 *
 * @param start_in first site of the gene in the HIV genome
 * @param end_in last site + 1 of the gene in the HIV genome
 */
hivgene::hivgene(unsigned int start_in, unsigned int end_in) {
		if((start_in>=HIVGENOME) || (end_in>HIVGENOME) || (end_in<=start_in))
			throw (int)HIVPOP_BADARG;
		else {
			start = start_in;
			end = end_in;
		}
}
