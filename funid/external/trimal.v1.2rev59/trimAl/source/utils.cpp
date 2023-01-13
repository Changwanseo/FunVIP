/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.2: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses.

    readAl v1.2: a tool for automated alignment conversion among different
                 formats.

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
                       [scapella, tgabaldon]@crg.es

    This file is part of trimAl/readAl.

    trimAl/readAl are free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl/readAl are distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl/readAl. If not, see <http://www.gnu.org/licenses/>.

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "utils.h"
#include "values.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++
| void utils::initVect(int *, int, int)          |
|      This method is used to initializate all   |
|      positions of a vector with a given value. |
++++++++++++++++++++++++++++++++++++++++++++++++*/

void utils::initlVect(int *vector, int tam, int valor) {

  for(int i = 0; i < tam; i++) vector[i] = valor;

}

void utils::initlVect(float *vector, int tam, float valor) {

  for(int i = 0; i < tam; i++) vector[i] = valor;

}


/*+++++++++++++++++++++++++++++++++++++++++++++
| void utils::copyVect(int *, int *, int)     |
|      This method copies integer vector 1 to |
|      integer vector 2.                      |
+++++++++++++++++++++++++++++++++++++++++++++*/

void utils::copyVect(int *vect1, int *vect2, int tam) {

  for(int i = 0; i < tam; i++) vect2[i] = vect1[i];

}


/*+++++++++++++++++++++++++++++++++++++++++++++++
| void utils::copyVect(float *, float *, float) |
|      This method copies float vector 1 to     |
|      float vector 2.                          |
+++++++++++++++++++++++++++++++++++++++++++++++*/

void utils::copyVect(float *vect1, float *vect2, int tam) {

  for(int i = 0; i < tam; i++) vect2[i] = vect1[i];

}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::roundToInf(double)           |
|      This method rounds a double number |
|      to the inferior integer.           |
+++++++++++++++++++++++++++++++++++++++++*/

int utils::roundToInf(double number) {

  return ((int) number);
}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::roundInt(double)             |
|      This method rounds a double number |
|      to a integer.                      | 
+++++++++++++++++++++++++++++++++++++++++*/

int utils::roundInt(double number) {    
 
  return ((int) ((double) number + 0.5));
}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::roundToSup(double)           |
|      This method rounds a double number |
|      to the greater integer.            | 
+++++++++++++++++++++++++++++++++++++++++*/

int utils::roundToSup(double number) {
  
  return ((int) ((double) number + 1.0));
}


/*+++++++++++++++++++++++++++++++++++++++++
| int utils::max(int, int)                |
|      This method returns the maximum    |
|      value of the two given arguments.  |
+++++++++++++++++++++++++++++++++++++++++*/

int utils::max(int x, int y) { 
  
  if(x > y) return x;
  else      return y;
}

float utils::max(float x, float y) { 
  
  if(x > y) return x;
  else      return y;
}

double utils::max(double x, double y) { 
  
  if(x > y) return x;
  else      return y;
}




/*+++++++++++++++++++++++++++++++++++++++++++
| bool utils::isNumber(char *)              |
|      This method checks if the given      |
|      string is a float number.            |
+++++++++++++++++++++++++++++++++++++++++++*/

bool utils::isNumber(char *num){
  
  int tam = strlen(num);
  int i, flt = 1, expn = 1, sgn = 1;

  for(i = 0; i < tam; i++) {
    if(num[i] == '.' && flt) 
      flt = 0;

    else if(((num[i] == 'e') ||(num[i] == 'E')) && expn)
      expn = 0;

    else if(((num[i] == '+') ||(num[i] == '-')) && sgn) {
      if(!expn) sgn = 0;
    }
    else if(num[i] > '9' || num[i] < '0')
      return false;
  }
  
  return true;
  
}


/*+++++++++++++++++++++++++++++++++++++++++++
| bool utils::compare(char *, char *)       |
|      This method compares the two strings |
|      given, and returns true if the two   |
|      strings are equal.                   |
+++++++++++++++++++++++++++++++++++++++++++*/

bool utils::compare(char *a, char *b){

  return(!strcmp(a,b));
}


/*++++++++++++++++++++++++++++++++++++++++++
| void utils::removeSpaces(char *, char *) |
|      This method removes spaces in the   |  
|      input string and put the result in  |
|      the output string.                  |
++++++++++++++++++++++++++++++++++++++++++*/

void utils::removeSpaces(char *in, char *out){

  unsigned int i, j = 0;

  for(i = 0; i < strlen(in); i++){

    if(in[i] != ' ' && in[i] != '\t'){
      out[j] = in[i];
      j++;
    }
  }
  out[j] = '\0';
}


/*++++++++++++++++++++++++++++++++++++++++++
| void utils::quicksort(float *, int, int) |
|      This method sorts the vector using  |
|      the quicksort method.               |
++++++++++++++++++++++++++++++++++++++++++*/

void utils::quicksort(float *vect, int ini, int fin) {

  float elem_div;
  int i, j;

  if ((ini >= fin) || (fin < 0))
    return;

  elem_div = vect[fin];
  i = ini - 1;
  j = fin;

  while (1) {

    while (vect[++i] < elem_div)
      if(i == fin)
        break;

    while (vect[--j] > elem_div)
      if(j == 0)
        break;

    if(i < j)
      swap(&vect[i], &vect[j]);
    else 
      break;
  }

  swap(&vect[i], &vect[fin]);

  quicksort(vect, ini, i - 1);
  quicksort(vect, i + 1, fin);
}


/*++++++++++++++++++++++++++++++++++++++++
| void utils::swap(float *, float *)     |
|      This method swaps the values in a |
|      and b.                            |
++++++++++++++++++++++++++++++++++++++++*/

void utils::swap(float *a, float *b){

  float temp;

  temp = *a;
  *a = *b;
  *b = temp;

}

/*++++++++++++++++++++++++++++++++++++++++++
| void utils::quicksort(float *, int, int) |
|      This method sorts the vector using  |
|      the quicksort method.               |
++++++++++++++++++++++++++++++++++++++++++*/

void utils::quicksort(int *vect, int ini, int fin) {

  int i, j, elem_div;

  if ((ini >= fin) || (fin < 0))
    return;

  elem_div = vect[fin];
  i = ini - 1;
  j = fin;

  while (1) {

    while (vect[++i] < elem_div)
      if(i == fin)
        break;

    while (vect[--j] > elem_div)
      if(j == 0)
        break;

    if(i < j)
      swap(&vect[i], &vect[j]);
    else 
      break;
  }

  swap(&vect[i], &vect[fin]);

  quicksort(vect, ini, i - 1);
  quicksort(vect, i + 1, fin);
}


/*++++++++++++++++++++++++++++++++++++++++
| void utils::swap(float *, float *)     |
|      This method swaps the values in a |
|      and b.                            |
++++++++++++++++++++++++++++++++++++++++*/

void utils::swap(int *a, int *b) {

  int temp;

  temp = *a;
  *a = *b;
  *b = temp;

}

/* ****************************************************************************************************************** */
/* ****************************************************************************************************************** */
