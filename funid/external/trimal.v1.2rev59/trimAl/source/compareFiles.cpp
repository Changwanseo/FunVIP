/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
   ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 

    trimAl v1.2: a tool for automated alignment trimming in large-scale 
                 phylogenetics analyses 

    Copyright (C) 2009 Capella-Gutierrez S. and Gabaldon, T.
                       [scapella, tgabaldon]@crg.es

    This file is part of trimAl.

    trimAl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, the last available version.

    trimAl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with trimAl. If not, see <http://www.gnu.org/licenses/>.

 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include "compareFiles.h"
#include "alignment.h"

#define LONG 80

void compareFiles::printStatisticsFileColumns(int numAminos, float *compareVect) {

  int i;

  cout << "| Residue\tConsistency |" << endl;
  cout << "| Number \t   Value    |" << endl;
  cout << "+---------------------------+" << endl;
  cout.precision(10);

  for(i = 0; i < numAminos; i++) 
    cout << "  " << setw(5) << i+1 << "\t\t" << compareVect[i] << endl;
}

void compareFiles::printStatisticsFileAcl(int numAminos, float *compareVect) {

  float refer, *vectAux = new float[numAminos];
  int i, num;

  utils::copyVect(compareVect, vectAux, numAminos);
  utils::quicksort(vectAux, 0, numAminos-1);


  /* We set the output precision and print the header. */
  cout << "| Number of\t        \t|\t Cumulative \t% Cumulative\t|   Consistency   |" << endl;
  cout << "| Residues \t% Length\t|\tNumberResid.\t   Length   \t|      Value      |" << endl;
  cout << "+-------------------------------+---------------------------------------+-----------------+" << endl;
  cout.precision(10);

  refer = vectAux[0];
  num = 1;

  for(i = 1; i < numAminos; i++) {

    if(refer != vectAux[i]) {
      cout << "  " << num << "\t\t" << setw(10) << ((float) num/numAminos * 100.0)
           << "\t\t" << i << "\t\t" << setw(10) << ((float) i/numAminos * 100.0) << "\t"
           << setw(15) << refer << endl;
      refer = vectAux[i]; 
      num = 1;
    }
    else num++;
  }

  cout << "  " << num << "\t\t" << setw(10) << ((float) num/numAminos * 100.0)
       << "\t\t" << i << "\t\t" << setw(10) << ((float) i/numAminos * 100.0) << "\t"
       << setw(15) << refer << endl;

  delete [] vectAux;
}

bool compareFiles::applyWindow(int columns, int halfWindow, float *columnsValue) {

  int i, j, window;
  float *vectAux;

  /* If halfWindow value is greater than 1/4 of alignment length, trimAl fails */
  if(halfWindow > columns/4)
     return false;
  else
    window = (2 * halfWindow + 1);

  /* Allocate local memory. Copy the array values to auxiliar memory. */
  vectAux = new float[columns];
  utils::copyVect(columnsValue, vectAux, columns);

  /* For each column from the alignment selected, compute a average value
     for its consistency value */
  for(i = 0; i < columns; i++) {
    /* This average is computed from halfWindow positions before to halfWindow
       positions after. */
    for(j = i - halfWindow, columnsValue[i] = 0; j <= i + halfWindow; j++) {
      if(j < 0)
        columnsValue[i] += vectAux[-j];
      else if(j >= columns)
        columnsValue[i] += vectAux[((2 * columns - j) - 2)];
      else
        columnsValue[i] += vectAux[j];
    }

    /* Finally, the column value is divided by the window size in order to
       compute the average score. */
    columnsValue[i] /= window;
  }

  delete [] vectAux;
  return true;
}



int compareFiles::algorithm(alignment **vectAlignments, char **fileNames, float *columnsValue, int numAlignments, bool verbosity) {

  int i, j, k, l, m, pairSeq = 0, acert = 0, alig = 0;
  int *numSpeciesAlig, *numAminosAlig;
  int *columnSeqMatrix = NULL, *columnSeqMatrixAux = NULL;
  float max = 0, value = 0, **vectAciertos;

  numSpeciesAlig = new int[numAlignments];
  numAminosAlig  = new int[numAlignments];
  vectAciertos   = new float*[numAlignments];

  for(i = 0; i < numAlignments; i++) {
    numSpeciesAlig[i] = vectAlignments[i] -> getNumSpecies();
    numAminosAlig[i] =  vectAlignments[i] -> getNumAminos();
  }

  for(i = 0; i < numAlignments; i++, value = 0) {

    if(verbosity)
      cout << endl;

    vectAciertos[i] = new float[numAminosAlig[i]];
    utils::initlVect(vectAciertos[i], numAminosAlig[i], 0);

    for(j = 0, k = 0; j < numAminosAlig[i]; j++, pairSeq = 0, acert = 0) {

      delete[] columnSeqMatrix;	
      columnSeqMatrix = new int[numSpeciesAlig[i]];
      utils::initlVect(columnSeqMatrix, numSpeciesAlig[i], 0);

      vectAlignments[i] -> getColumnSeqMatrix(j, columnSeqMatrix);

      for(k = 0; k < numSpeciesAlig[i]; k++) {

        if(columnSeqMatrix[k] != 0) {

          for(l = 0; l < numAlignments; l++) {
            if((l != i) && (k < numSpeciesAlig[l])) {

              delete[] columnSeqMatrixAux;	
              columnSeqMatrixAux = new int[numSpeciesAlig[l]];
              utils::initlVect(columnSeqMatrixAux, numSpeciesAlig[l], 0);

              vectAlignments[l] -> getColumnSeqMatrix(columnSeqMatrix[k], k, columnSeqMatrixAux);

              for(m = k+1; (m < numSpeciesAlig[i]) && (m < numSpeciesAlig[l]); m++)
                if(columnSeqMatrix[m] != 0) {
                  if(columnSeqMatrix[m] == columnSeqMatrixAux[m]) 
                    acert++;
                  pairSeq++;
                }
            }
          }
        } 
      }

      if(pairSeq != 0) {
         vectAciertos[i][j] += ((1.0 * acert)/pairSeq);
         value += vectAciertos[i][j];
      }
    }		

    if(verbosity)
      cout << "File:\t\t" << fileNames[i] << endl << "Values:\t\tSequences: " 
           << numSpeciesAlig[i] << "\tResidues: " << numAminosAlig[i] <<  "\tPond. Hits: " << setw(8) 
           << value << "\t%Consistency: " << value/numAminosAlig[i] << endl;

    if((value/numAminosAlig[i]) > max) {
      alig = i;
      max = value/numAminosAlig[i];
    } 
  }

  if(verbosity) {
    cout << "\t\t\t\t\t--------------" << endl;
    cout << endl << "File Selected:\t" << fileNames[alig] << endl << "Value:\t\t" << max << endl << endl;
  }

  if(columnsValue != NULL)
    for(i = 0; i < numAminosAlig[alig]; i++) 
      columnsValue[i] = vectAciertos[alig][i];

  for(i = 0; i < numAlignments; i++)
    delete [] vectAciertos[i];

  delete[] columnSeqMatrix;
  delete[] columnSeqMatrixAux;

  delete [] numSpeciesAlig;
  delete [] numAminosAlig;
  delete [] vectAciertos;

  return alig;
}

void compareFiles::compare(alignment **vectAlignments, char **fileNames, int numAlignments) {
  int i, j, k, l, ll, lng, l3ng, l2ng = 0, inc = 1, tamNames = 0, *OthLongs, **mapOthAligs;
  char *str = NULL, *s2tr = NULL, **vectAcrt = NULL, line[256];

  lng = vectAlignments[0] -> getNumSpecies() + 1;
  str  = new char[lng]; s2tr = new char[lng];

  vectAcrt = new char*[numAlignments];
  lng = vectAlignments[0] -> getNumAminos();

  for(i = 0; i < numAlignments; i++)
    vectAcrt[i] = new char[lng];

  for(i = 0; i < numAlignments; i++)
    tamNames = (tamNames > (int) strlen(fileNames[i])) ? tamNames : (int) strlen(fileNames[i]);

  OthLongs = new int[numAlignments];
  mapOthAligs = new int*[numAlignments];

  mapOthAligs[0] = NULL;

  for(i = 1; i < numAlignments; i++) {
    OthLongs[i] = vectAlignments[i] -> getNumAminos();
    mapOthAligs[i] = new int[OthLongs[i]];
    utils::initlVect(mapOthAligs[i], OthLongs[i], -1);
  }

  for(i = 0; i < lng; i += LONG, inc = 0) {
    l2ng = vectAlignments[0] -> alignmentToFile(cout, inc, i, LONG, tamNames);

    for(j = i; j < (i + LONG) && j < lng; j++) {
      vectAlignments[0] -> getColumn(j, str);

      for(k = 1; k < numAlignments; k++) {
        vectAcrt[k][j] = '-';

        for(l = 0; l < OthLongs[k]; l++) {
          vectAlignments[k] -> getColumn(l, s2tr);
          if((!strcmp(str, s2tr)) && (mapOthAligs[k][l] == -1)) {
            vectAcrt[k][j] = '*';
            mapOthAligs[k][l] = j;
            break; 
          }
        }

        for(ll = l+1; ll < OthLongs[k]; ll++)
          if(mapOthAligs[k][ll] < j) { 
            vectAcrt[k][mapOthAligs[k][ll]] = '-';
            mapOthAligs[k][ll] = -1;
          }
      }
    }

    for(k = 1; k < numAlignments; k++) {
      strncpy(line, fileNames[k], l2ng);
      line[l2ng] = '\0'; cout << setw(l2ng+2) << line;

      if((i + LONG) > lng) l3ng = lng - i;
      else l3ng = LONG;

      strncpy(line, &vectAcrt[k][i], l3ng);
      line[l3ng] = '\0'; cout << line << endl;
    }
    cout << endl;
  }

  cout << endl << "Summary: " << endl;
  cout << "Reference Sequence:\t" << fileNames[0] << "\thas " << lng << " positions" << endl << endl;
  cout << "Compare Sequence/s:\t" << endl;

  for(i = 1, k = 0; i < numAlignments; i++, k = 0) {
    for(j = 0; j < OthLongs[i]; j++)
      if(mapOthAligs[i][j] != -1) k++;

    strncpy(line, fileNames[i], l2ng);
    line[l2ng] = '\0'; cout << "\t\t\t" << setw(l2ng+2) << line;
    cout << "\thas " << OthLongs[i] << " positions\tFound: " << k << " positions" << endl;
  }
  cout << endl;

}
