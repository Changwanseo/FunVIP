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

#include "alignment.h"

/* ********************************************************************************************************* */
/* ********************************************************************************************************* */
/* ********** 					NEW CODE 					  ********** */
/* ********************************************************************************************************* */
/* ********************************************************************************************************* */

#define GAPPYOUT 1
#define STRICT 2

int alignment::calculateSeqIdentity(void) {

  char indet;
  int i, j, k, hit, dst;
  float mx, avg, maxSeq = 0, avgSeq = 0, **values;

  values = new float*[sequenNumber];

  if(getTypeAlignment() == AAType)
    indet = 'X';
  else
    indet = 'N';

  for(i = 0; i < sequenNumber; i++) {
    values[i] = new float[sequenNumber];

    for(k = 0; k < i; k++) 
      values[i][k] = values[k][i];
    values[i][i] = 0;

    for(k = i + 1; k < sequenNumber; k++) {
      for(j = 0, hit = 0, dst = 0; j < residNumber; j++) {
        if(((alignmentMatrix[i][j] != indet) && (alignmentMatrix[i][j] != '-')) ||
           ((alignmentMatrix[k][j] != indet) && (alignmentMatrix[k][j] != '-'))) {
          dst++;
          if(alignmentMatrix[i][j] == alignmentMatrix[k][j]) 
            hit++;
        }
      }
      values[i][k] = (float) hit/dst;
    }

    for(k = 0, mx = 0, avg = 0; k < sequenNumber; k++) {
      if(i != k) {
        mx  = mx < values[i][k] ? values[i][k] : mx;
        avg += values[i][k];
      }
    }

    avgSeq += avg/(sequenNumber - 1);
    maxSeq += mx;
  }

  avgSeq = avgSeq/sequenNumber;
  maxSeq = maxSeq/sequenNumber;

  if(avgSeq >= 0.55)
    return GAPPYOUT;

  else if(avgSeq <= 0.38)
    return STRICT;

  else {

    if(sequenNumber <= 20)
      return GAPPYOUT;

    else {

      if((maxSeq >= 0.5) && (maxSeq <= 0.65))
        return GAPPYOUT;

      else
        return STRICT;
    }
  }
}

void alignment::printSeqIdentity(void) {

  char indet;
  int i, j, k, hit, dst, pos, maxLongName;
  float mx, avg, maxSeq = 0, avgSeq = 0, **maxs, **avgs;

  avgs = new float*[sequenNumber];
  maxs = new float*[sequenNumber];

  if(getTypeAlignment() == AAType)
    indet = 'X';
  else
    indet = 'N';

  for(i = 0; i < sequenNumber; i++) {
    avgs[i] = new float[sequenNumber];
    maxs[i] = new float[2];

    for(k = 0; k < i; k++) 
      avgs[i][k] = avgs[k][i];
    avgs[i][i] = 0;

    for(k = i + 1; k < sequenNumber; k++) {
      for(j = 0, hit = 0, dst = 0; j < residNumber; j++) {
        if(((alignmentMatrix[i][j] != indet) && (alignmentMatrix[i][j] != '-')) ||
           ((alignmentMatrix[k][j] != indet) && (alignmentMatrix[k][j] != '-'))) {
          dst++;
          if(alignmentMatrix[i][j] == alignmentMatrix[k][j])
            hit++;
        }
      }
      avgs[i][k] = (float) hit/dst;
    }

    for(k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
      if(i != k) {
        if(mx < avgs[i][k]) {
          mx  = avgs[i][k];
          pos = k;
        }
        avg += avgs[i][k];
      }
    }

    avgSeq += avg/(sequenNumber - 1);
    maxSeq += mx;

    maxs[i][0] = mx;
    maxs[i][1] = pos;
  }

  avgSeq = avgSeq/sequenNumber;
  maxSeq = maxSeq/sequenNumber;

  for(i = 0, maxLongName = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, strlen(sequenNames[i]));

  cout.precision(6);
  cout << endl << "#Mean Percentage of identity:                            " << avgSeq;
  cout << endl << "#Mean Percentage of identity with most similar sequence: " << maxSeq;

  cout << endl << endl << "#Percentage of identity matrix:";
  for(i = 0; i < sequenNumber; i++) {
    cout << endl << setw(maxLongName + 2) << left << sequenNames[i] << "\t";
    for(j = 0; j < sequenNumber; j++)
      cout << setiosflags(ios::left) << setw(10) << avgs[i][j] * 100 << "\t";
  }

  cout << endl << endl << "#Percentage of identity with most similar sequence:" << endl;

  for(i = 0; i < sequenNumber; i++)
    cout << setw(maxLongName + 2) << left << sequenNames[i] << "\t" << setiosflags(ios::left) << setw(5) << maxs[i][0] * 100
         << "\t\t" << sequenNames[(int) maxs[i][1]] << endl;

  cout << endl;

}
