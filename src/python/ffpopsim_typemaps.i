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

/*******************************************************************/
/* Typemaps for FFPopSim                                           */
/*******************************************************************/
/* convert between Python bool array and boost::dynamic_bitset */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  boost::dynamic_bitset<>
{
        $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in) boost::dynamic_bitset<> (boost::dynamic_bitset<> temp)
{
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Create boost::dynamic_bitset from Python list */
        temp.resize(L);
        long tmplong;
        for(size_t i=0; i < L; i++) {
                tmplong = PyInt_AsLong(PySequence_Fast_GET_ITEM(tmplist, i));
                if(tmplong < 0) {
                        PyErr_SetString(PyExc_ValueError, "Expecting an array of bool.");
                        SWIG_fail;
                }
                temp[i] = (bool)tmplong; 
        }      
        $1 = temp;
}
%typemap(out) boost::dynamic_bitset<>
{
        unsigned long L = $1.size();
        npy_intp dims[1] = {(npy_intp) L};
        PyObject *array = PyArray_ZEROS(1, dims, NPY_BOOL, 0);
        if (!array) SWIG_fail;

        /* no checks on memory alignments, since we create a new array */
        bool *ptr = (bool *)PyArray_DATA(array);
        for(int i=0; i < L; i++, ptr++)
                if($1.test(i))
                        *ptr = true;
        $result = array;
}

/* STL::Vector (FIX) */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  std::vector<double>
{
        $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in) std::vector<double> (std::vector<double> temp)
{
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Create boost::dynamic_bitset from Python list */
        temp.reserve(L);
        for(size_t i=0; i < L; i++) {
                temp.push_back(PyFloat_AsDouble(PySequence_Fast_GET_ITEM(tmplist, i)));
                if(PyErr_Occurred()) {
                        PyErr_SetString(PyExc_ValueError, "Expecting an array of double.");
                        SWIG_fail;
                }
        }      
        $1 = temp;


}
%typemap(out) std::vector<double>
{
        unsigned long L = $1.size();
        npy_intp dims[1] = {(npy_intp) L};
        PyObject *array = PyArray_ZEROS(1, dims, NPY_FLOAT, 0);
        if (!array) SWIG_fail;

        /* no checks on memory alignments, since we create a new array */
        double *ptr = (double *)PyArray_DATA(array);
        for(std::vector<double>::iterator iter = $1.begin(); iter != $1.end(); iter++, ptr++)
                *ptr = *iter;
        $result = array;
}

/* array of length L */
%typemap(in,
         fragment="NumPy_Fragments")
  double* IN_ARRAY1_L
  (PyArrayObject* array=NULL, int is_new_object=0)
{
  npy_intp size[1] = { -1 };
  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,
                                                   &is_new_object);
  if (!array || !require_dimensions(array, 1) ||
      !require_size(array, size, 1)) SWIG_fail;
  $1 = (double*) array_data(array);
}
%typemap(freearg)
  double* IN_ARRAY1_L
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

/* streams and strings */
%typemap(in) istream &model (std::ifstream temp) {
        if (!PyString_Check($input)) {
                PyErr_SetString(PyExc_ValueError, "Expecting a string");
                return NULL;
        }
        temp.open(PyString_AsString($input));
        $1 = &temp;
}
%typemap(in) ostream &out_genotypes (std::ofstream temp) {
        if (!PyString_Check($input)) {
                PyErr_SetString(PyExc_ValueError, "Expecting a string");
                return NULL;
        }
        temp.open(PyString_AsString($input));
        $1 = &temp;
}


/* LOWD */
/* recombination rates */
%typemap(in) double *rec_rates {
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Get circular and L properties from the class (we are in the Python world here) */
        bool circular = (bool)PyInt_AsLong(PyObject_GetAttrString($self, "circular"));
        unsigned long Lint = PyInt_AsLong(PyObject_GetAttrString($self, "L"));

        /* Check lengths */
        if((!(circular)) && (L != Lint - 1)) {
                PyErr_SetString(PyExc_ValueError, "Expecting an array of length L-1.");
                SWIG_fail;
        }        
        if((circular) && (L != Lint)) {
                PyErr_SetString(PyExc_ValueError, "Expecting an array of length L.");
                SWIG_fail;
        } 

        /* Create C array from Python list */
        $1 = new double[L];
        double tmpdouble;
        for(size_t i=0; i < L; i++) {
                tmpdouble = (double)PyFloat_AsDouble(PySequence_Fast_GET_ITEM(tmplist, i));
                if (tmpdouble < 0) {
                        PyErr_SetString(PyExc_ValueError,"Expecting a sequence of positive floats");
                        SWIG_fail;
                }
                $1[i] = tmpdouble;
        }
}
%typemap(freearg) double *rec_rates {
  if($1) delete[] $1;
}


/* HIGHD */
/* vectors of multilocus coefficients */
%typemap(in) vector<int> loci (std::vector<int> temp) {
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Create std::vector from Python list */
        temp.reserve(L);
        long tmplong;
        for(size_t i=0; i < L; i++) {
                tmplong = PyInt_AsLong(PySequence_Fast_GET_ITEM(tmplist, i));
                if(tmplong < 0) {
                        PyErr_SetString(PyExc_ValueError, "Expecting an array of positive integers (the loci).");
                        SWIG_fail;
                }
                temp.push_back((int)tmplong); 
        }      
        $1 = temp;
}
%typemap(out) vector<coeff_t> {
        PyObject *tmpitem, *tmploci;
        double tmpval;
        int tmporder;
        coeff_t *tmpcoeff;

        PyObject *tmplist = PyList_New(0);
        for(int i=0; i < $1.size(); i++) {
                /* set the i-th coefficient */
                tmpcoeff = &($1.at(i));
                tmpval = tmpcoeff->value;
                tmporder = tmpcoeff->order; 
                tmpitem = PyTuple_New(2);
                /* set the value */
                PyTuple_SET_ITEM(tmpitem, 0, PyFloat_FromDouble(tmpval));
                /* set the loci */
                tmploci = PyTuple_New(tmporder);
                for(int j=0; j < tmporder; j++)
                        PyTuple_SET_ITEM(tmploci, j, PyInt_FromLong((long)(tmpcoeff->loci[j])));
                PyTuple_SET_ITEM(tmpitem, 1, tmploci);
                /* set the coeff */
                PyList_Append(tmplist, tmpitem);
        }
        $result = tmplist;
}


/* typemaps extensions: CHECK HERE FOR INPUT ARGUMENT SIGNATURES! */
%apply double* IN_ARRAY1_L {double* frequencies};
