/**
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

/* licence */
%pythoncode{
LICENSE = '''FFPopSim is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. FFPopSim is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with FFPopSim. If not, see <http://www.gnu.org/licenses/>.'''
}

/* NUMPY */
%pythoncode {
import numpy as _np
}

/* renames and ignores */
%ignore SAMPLE_ERROR;
%ignore sample;


/*****************************************************************************/
/* INDEX_VALUE_PAIR_T                                                        */
/*****************************************************************************/
%feature("autodoc", "Pair of an index and a value") index_value_pair_t;
%rename(index_value_pair) index_value_pair_t;
%extend index_value_pair_t {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"index: %u, val: %.2e", (unsigned int)$self->index,
                                               $self->val);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"(%u, %.2e)", (unsigned int)$self->index, $self->val);
        return &buffer[0];
}

%feature("autodoc", "Index") index;
%feature("autodoc", "Value") val;
} /* extend index_value_pair_t */
/*****************************************************************************/

/*****************************************************************************/
/* GENOTYPE_VALUE_PAIR_T                                                     */
/*****************************************************************************/
%feature("autodoc", "Pair of a genotype and a value") genotype_value_pair_t;
%rename(genotype_value_pair) genotype_value_pair_t;
%extend genotype_value_pair_t {
/* string representations */
const char* __str__() {
        static char buffer[255];
        unsigned long L = ($self->genotype).size();
        if(L > 2)
                sprintf(buffer,"genotype: %d...%d, val: %.2e",
                               (int)($self->genotype)[0],
                               (int)($self->genotype)[L-1],
                               $self->val);
        else if(L == 2)
                sprintf(buffer,"genotype: %d%d, val: %.2e",
                               (int)($self->genotype)[0],
                               (int)($self->genotype)[1],
                               $self->val);
        else if(L == 1)
                sprintf(buffer,"genotype: %d, val: %.2e",
                               (int)($self->genotype)[0],
                               $self->val);
        else
                sprintf(buffer,"genotype: (empty), val: %.2e",
                               $self->val);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        unsigned long L = ($self->genotype).size();
        if(L > 2)
                sprintf(buffer,"([%d, ..., %d], %.2e)",
                               (int)($self->genotype)[0],
                               (int)($self->genotype)[L-1],
                               $self->val);
        else if(L == 2)
                sprintf(buffer,"([%d, %d], %.2e)",
                               (int)($self->genotype)[0],
                               (int)($self->genotype)[1],
                               $self->val);
        else if(L == 1)
                sprintf(buffer,"([%d], %.2e)",
                               (int)($self->genotype)[0],
                               $self->val);
        else
                sprintf(buffer,"([], %.2e)", $self->val);
        return &buffer[0];
}

%feature("autodoc", "Value") val;
%feature("autodoc", "Genotype") genotype;
} /* extend genotype_value_pair_t */
/*****************************************************************************/

/*****************************************************************************/
/* STAT_T                                                                    */
/*****************************************************************************/
%feature("autodoc", "Mean and variance of a statistical distribution") stat_t;
%rename(stat) stat_t;
%extend stat_t {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"mean: %e, variance: %e", $self->mean, $self->variance);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"(%e, %e)", $self->mean, $self->variance);
        return &buffer[0];
}

%feature("autodoc", "Mean") mean;
%feature("autodoc", "Variance") variance;
} /* extend stat_t */
/*****************************************************************************/
