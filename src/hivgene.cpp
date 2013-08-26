/**
 * @file hivgene.cpp
 * @brief Implementation of an HIV gene
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-05-18
 *
 * Copyright (c) 2012-2013, Richard Neher, Fabio Zanini
 * All rights reserved.
 *
 * This file is part of FFPopSim.
 *
 * FFPopSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FFPopSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FFPopSim. If not, see <http://www.gnu.org/licenses/>.
 */
#include "hivpopulation.h"

/**
 * @brief Default constructor for an HIV gene
 *
 * @param start_in first site of the gene in the HIV genome
 * @param end_in last site + 1 of the gene in the HIV genome
 */
hivgene::hivgene(unsigned int start_in, unsigned int end_in, unsigned int second_start_in, unsigned int second_end_in) {
		if((start_in>=HIVGENOME) || (end_in>HIVGENOME) || (end_in<=start_in) ||
		   (second_start_in>=HIVGENOME) || (second_end_in>HIVGENOME) || (second_end_in<second_start_in) ||
		   ((second_end_in==second_start_in) && (second_start_in!=0)))
			throw (int)HIVPOP_BADARG;
		else {
			start = start_in;
			end = end_in;
			second_start = second_start_in;
			second_end = second_end_in;
		}
}
