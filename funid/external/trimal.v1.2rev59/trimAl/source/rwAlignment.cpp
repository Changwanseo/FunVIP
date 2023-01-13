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

#include "alignment.h"
#include "utils.h"

extern int errno;
#include <errno.h>
#include <ctype.h>
#include <string>

#define DELIMITERS     " \t\n"
#define OTHDELIMITERS  " \t\n,:*"
#define OTH2DELIMITERS " \t\n,:;*"
#define PHYLIPDISTANCE 10

/* ********************************************************************************************************* */
/* ********************************************************************************************************* */
/* ********** 					NEW CODE 					  ********** */
/* ********************************************************************************************************* */
/* ********************************************************************************************************* */

  using namespace std;

bool alignment::fillMatrices(string *Sequences, bool aligned) {
  int i, j;

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  residuesNumber  = new   int[sequenNumber];
  alignmentMatrix = new char*[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++) {
    residuesNumber[i] = Sequences[i].size();
    alignmentMatrix[i] = new char[residuesNumber[i] + 1];
    strcpy(alignmentMatrix[i], Sequences[i].c_str());
  }

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if(sequenNumber == 1) isAligned = true;
  for(i = 1; i < sequenNumber; i++) {
    isAligned = false;
    if(residuesNumber[i] != residuesNumber[i-1]) break;
    else isAligned = true;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** Fill some info ***** ***** ***** */
  if(residNumber == 0)
    residNumber = residuesNumber[0];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** Check content ***** ***** ***** */
  if(aligned) {
    for(i = 0; i < sequenNumber; i++) {
      if(residuesNumber[i] != residNumber) {
        cerr << endl << "ERROR: The sequence \"" << sequenNames[i] << "\" (" << residuesNumber[i] 
             << ") does not have the same number of residues fixed by the alignment (" << residNumber << ").";
        return false;
      }
    }
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  if((aligned) || (isAligned)) {

    /* Asign its position to each column. */
    saveResidues = new int[residNumber];
    for(i = 0; i < residNumber; i++) saveResidues[i] = i;

    /* Asign its position to each sequence. */
    saveSequences = new int[sequenNumber];
    for(i = 0; i < sequenNumber; i++) saveSequences[i] = i;
  }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  for(i = 0; i < sequenNumber; i++)
    for(j = 0; j < residuesNumber[i]; j++)
      if((!isalpha(alignmentMatrix[i][j])) && (!ispunct(alignmentMatrix[i][j]))) {
        cerr << endl << "ERROR: The sequence \"" << sequenNames[i] << "\" has an unknow (" << alignmentMatrix[i][j]
             <<") character.";
        return false;
      }
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  return true;
}

int alignment::formatInputAlignment(char *alignmentFile) {

  char c, *firstWord , *line = NULL;
  int state, format = 0, blocks = 0;
  long begin, end;
  ifstream file;
  string nline;

  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return -1;

  /* ***** Check the alignment content ***** */
  begin = file.tellg();
  file.seekg (0, ios::end);
  end = file.tellg();
  if(!(end - begin)) return -1;
  file.seekg (0, ios::beg);
  /* ***** ***** ***** ***** ***** ***** ***** */

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;

  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  delete [] line;
  line = new char [nline.size()+1];
  strcpy(line, nline.c_str());
  /* ******************************************************************** */
  firstWord = strtok(line, OTHDELIMITERS);

  /* Clustal Format */
  if((!strcmp(firstWord, "CLUSTAL")) || (!strcmp(firstWord, "clustal")))
    format = 1;

  /* NBRF/PIR Format */
  else if(firstWord[0] == '>' && firstWord[3] == ';')
    format = 3;

  /* Fasta Format */
  else if(firstWord[0] == '>')
    format = 8;

  /* Nexus Format */
  else if((!strcmp(firstWord, "#NEXUS")) || (!strcmp(firstWord, "#nexus")))
    format = 17;

  /* Mega Format */
  else if((!strcmp(firstWord, "#MEGA")) || (!strcmp(firstWord, "#mega"))) {

    blocks = 0;
    do{ file.read(&c, 1); } while((c != '#') && (!file.eof()));

    do {
      while((c != '\n') && (!file.eof())) file.read(&c, 1);
      file.read(&c, 1); if(c == '#') blocks++;
    } while((c != '\n') && (!file.eof()));

    /* Mega NonInterleaved */
    if(!blocks) format = 22;
    /* Mega Interleaved */
    else format = 21;
  }

  /* Phylip Format */
  else {

    sequenNumber = atoi(firstWord);
    firstWord = strtok(NULL, DELIMITERS);

    if(firstWord != NULL) residNumber = atoi(firstWord);

    if((sequenNumber != 0) && (residNumber != 0)) {

    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    delete [] line;
    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    firstWord = strtok(line, DELIMITERS);
    if(firstWord != NULL) blocks = 1;

    while(1) {
      firstWord = strtok(NULL, DELIMITERS);
      if(firstWord != NULL) blocks++;
      else break;
    }

    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    delete [] line;
    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    firstWord = strtok(line, DELIMITERS);
    if(firstWord != NULL) blocks--;

    while(1) {
      firstWord = strtok(NULL, DELIMITERS);
      if(firstWord != NULL) blocks--;
      else break;
    }

    /* Phylip Format */
    if(!blocks) format = 12;

    /* Phylip3.2 Format */
    else format = 11;
    }
  }

  delete [] line;
  file.close();
  return format;
}


bool alignment::loadPhylipAlignment(char *alignmentFile) {

  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;
  int i, state;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL)
    freeAlignment();
  /* ******************************************************************** */

  /* ******************************************************************** */
  do{ file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;

  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  delete [] line;
  line = new char [nline.size()+1];
  strcpy(line, nline.c_str());
  /* ******************************************************************** */

  /* ******************************************************************** */
  str = strtok(line, DELIMITERS);
  if(str != NULL) sequenNumber = atoi(str);
  else return false;

  str = strtok(NULL, DELIMITERS);
  if(str != NULL) residNumber = atoi(str);
  else return false;

  /* ******************************************************************** */
  /* Check correct read parameters */
  if((sequenNumber == 0) || (residNumber == 0)) 
    return false;

  /* Allocate memory for alignment matrix and sequenNumber names vector */
  sequenNames = new  char*[sequenNumber];
  seqs         = new string[sequenNumber];

  /* Read the first block of lines that contains two tokens */
  for(i = 0; i < sequenNumber; i++) {

    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    /* First token: sequenNumber name */
    str = strtok(line, DELIMITERS);
    sequenNames[i] = new char[strlen(str) + 1];
    strcpy(sequenNames[i], str);

    /* Next token: first N aminoacid block */
    str = strtok(NULL, DELIMITERS);
    seqs[i].append(str, strlen(str));

    /* we clean the line of spaces and tabulations */
    /* and store it in the alignment matrix        */
    while(1) {
      str = strtok(NULL, DELIMITERS);
      if(str != NULL) seqs[i].append(str, strlen(str));
      else break;
    }
  }

  /* Read the next aminoacid blocks */
  while(!file.eof()) {
    for(i = 0; i < sequenNumber; i++) {

      /* ******************************************************************** */
      do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
      for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
        nline += c;

      if(file.eof()) break;
      delete [] line;

      state = nline.find("\r", 0);
      while(state != (int) string::npos) {
        nline.erase(state, 1);
        state = nline.find("\r", (state + 1));
      }

      line = new char [nline.size() + 1];
      strcpy(line, nline.c_str());
      /* ******************************************************************** */

      str = strtok(line, DELIMITERS);
      seqs[i].append(str, strlen(str));

      while(1) {
        str = strtok(NULL, DELIMITERS);
        if(str != NULL) seqs[i].append(str, strlen(str));
        else break;
      }
    }
  }

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, true);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

bool alignment::loadFastaAlignment(char *alignmentFile) {

  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;
  int i, state;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL)
    freeAlignment();
  /* ******************************************************************** */

  while(1) {
    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); ((c != '\n') && (!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());

    /* ******************************************************************** */
    str = strtok(line, DELIMITERS);
    if(str[0] == '>') sequenNumber++;
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* Allocate memory for sequenNumber names vector */
  sequenNames =  new char*[sequenNumber];
  seqs        = new string[sequenNumber];

  for(i = -1; i < sequenNumber; ) {
    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); ((c != '\n') && (!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    /* First token: sequenNumber name */
    str = strtok(line, OTHDELIMITERS);
    if(str[0] == '>') {
      i++;
      if(strlen(str) == 1) str = strtok(NULL, OTHDELIMITERS);
      else str = str + 1;

      sequenNames[i] = new char[strlen(str) + 1];
      strcpy(sequenNames[i], str);
    }
    else {
      while(1) {
        if(str == NULL) break;
        else seqs[i].append(str, strlen(str));
        str = strtok(NULL, DELIMITERS);
      }
    }
  }

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, false);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}


bool alignment::loadClustalAlignment(char *alignmentFile) {

  int i, state, length = 0, firstBlock = true;
  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL)
   freeAlignment();
  /* ******************************************************************** */

  /* ***** ***** ***** ***** Title ***** ***** ***** ***** */
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  /* ***** ***** ***** ***** ***** ***** ***** ***** ***** */

  while(1) {
    /* ******************************************************************** */
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    for(i = 0, length = 0; i < (int) strlen(line); i++)
      if((isalpha(line[i])) || (line[i] == '-')) length++;

    /* ******************************************************************** */
    if(!length) break;
    sequenNumber++;
    file.read(&c, 1);
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  sequenNames = new  char*[sequenNumber];
  seqs        = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;

  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  delete [] line;
  line = new char [nline.size()+1];
  strcpy(line, nline.c_str());
  /* ******************************************************************** */
  alignmentInfo = new char[strlen(line) + 1];
  strcpy(alignmentInfo, line);

  /* ******************************************************************** */
  firstBlock = true; i = 0;

  while(1) {
     if(i >= sequenNumber) i = 0;
    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */
    state = i;
    for(i = 0, length = 0; i < (int) strlen(line); i++)
      if((isalpha(line[i])) || (line[i] == '-')) length++;
    i = state;
    /* ******************************************************************** */
    if(length) {
      str = strtok(line, OTHDELIMITERS);

      if(str != NULL) {
        if(firstBlock) {
          sequenNames[i] = new char[strlen(str) + 1];
          strcpy(sequenNames[i], str);
        }

        str = strtok(NULL, OTHDELIMITERS);
        if(str != NULL) seqs[i].append(str, strlen(str));
        i++;
      }
    } else
        firstBlock = false;
  }

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, true);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

bool alignment::loadNexusAlignment(char *alignmentFile) {

  char c, *frag, *str, *line = NULL;
  int i, state = false, firstBlock = true;
  string nline, auxline, *seqs;
  ifstream file;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL) 
    freeAlignment();
  /* ******************************************************************** */

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  while((c != '\n') && (!file.eof())) { file.read(&c, 1); }
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  /* ******************************************************************** */

  /* ******************************************************************** */
  do {

    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    delete [] line;
    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());

    str = strtok(line, DELIMITERS);
    /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
    if(str != NULL) {

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      for(i = 0; i < (int) strlen(str); i++)
        str[i] = toupper(str[i]);

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      if(!strcmp(str, "MATRIX")) break;

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      if(!strcmp(str, "BEGIN")) state = true;

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      if((!strcmp(str, "DIMENSIONS")) && true) {
        str = strtok(NULL, DELIMITERS);
        frag = strtok(NULL, DELIMITERS);

        str = strtok(str, "=;");
        sequenNumber = atoi(strtok(NULL, "=;"));

        frag = strtok(frag, "=;");
        residNumber = atoi(strtok(NULL, "=;"));

        if(auxline.size() > 0) auxline += ' ' + nline;
        else auxline = nline;
      }

      /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
      if((!strcmp(str, "FORMAT")) && true) {
        if(auxline.size() > 0) auxline += ' ' + nline;
        else auxline = nline;
      }
    }
    /* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));

  } while((c != '\n') && (!file.eof()));

  /* ******************************************************************** */
  /* Check correct read parameters */
  if(strcmp(str, "MATRIX") || (sequenNumber == 0) || (residNumber == 0))
    return false;
  /* ******************************************************************** */

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  sequenNames =  new char*[sequenNumber];
  seqs        = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  i = 0;
  while(1) {
    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size()+1];
    strcpy(line, nline.c_str());

    if((!nline.compare(0,4,"end;")) || (!nline.compare(0,4,"END;"))) break;
    /* ******************************************************************** */
    str = strtok(line, OTH2DELIMITERS);
    if(str != NULL) {
      if(firstBlock) {
        sequenNames[i] = new char[strlen(str) + 1];
        strcpy(sequenNames[i], str);
      }

      while(1) {
        str = strtok(NULL, OTH2DELIMITERS);
        if(str != NULL) seqs[i].append(str, strlen(str));
        else break;
      }
      i++;
    }
    if(i >= sequenNumber) {
      firstBlock = false;
      i = 0;
    }
  }

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, true);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

bool alignment::loadMegaInterleavedAlignment(char *alignmentFile) {

  int i, pos, npos, next, state, firstBlock = true;
  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL) 
    freeAlignment();
  /* ******************************************************************** */

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  /* ******************************************************************** */
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;

  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  line = new char [nline.size()+1];
  strcpy(line, nline.c_str());
  /* ******************************************************************** */
  str = strtok(line, "!: ");
  for(i = 0; i < (int) strlen(str); i++)
    str[i] = toupper(str[i]);
  /* ******************************************************************** */
  if(!strcmp(str, "TITLE")) {
    delete [] inputFileName;

    if(nline[0] == '!') {
      inputFileName = new char [nline.size() + 1];
      inputFileName[nline.size()] = '\0';
      strcpy(inputFileName, nline.c_str());
    }
    else {
      inputFileName = new char [nline.size() + 2];
      inputFileName[nline.size() + 1] = '\0';
      strcpy(inputFileName, "!");
      strcat(inputFileName, nline.c_str());
    }
    nline.clear();
  }

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  nline.clear(); nline += c;

  while((c != '\n') && (!file.eof())) {
   for(file.read(&c, 1); (c != '\n') && (!file.eof()) && (c != '#'); file.read(&c, 1))
     nline += c;

   if(file.eof() || (c == '#')) break;
   nline += c;
   file.read(&c, 1);
  }
  /* ******************************************************************** */
  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  alignmentInfo = new char [nline.size() + 1];
  strcpy(alignmentInfo, nline.c_str());
  /* ******************************************************************** */
  while((c != '#') && (!file.eof())) { file.read(&c, 1); }
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;
  /* ******************************************************************** */

  while(1) {
    /* ******************************************************************** */
    if(nline[0] != '!') {
      next = nline.size() + 1; pos  = -1;

      /* ***** ***** We need to do this because it could be some ***** ***** 
         ***** ***** comment in the middle of the sequence ***** *****  */
      while(1) {
        pos  = nline.find("\"", (pos + 1));
        npos = nline.find("\"", (pos + 1));
        next = nline.rfind("\"", next - 1);

        if(pos != (int) string::npos) {
          if(npos == next) {
            nline.erase(pos, (next - pos + 1));
            pos = -1; next = nline.size() + 1;
          }
        } else break;
      }
      do {
        pos = -1; npos = pos;
        while((pos = nline.find("[", (pos + 1))) != (int) string::npos) npos = pos;
        if((int) nline.find("]", (npos + 1)) != (int) string::npos)
          nline.erase(npos, (nline.find("]", (npos + 1)) - npos + 1));
      } while(npos != -1);
      /* ******************************************************************** */
      state = nline.find("\r", 0);
      while(state != (int) string::npos) {
        nline.erase(state, 1);
        state = nline.find("\r", (state + 1));
      }

      delete [] line;
      line = new char [nline.size()+1];
      strcpy(line, nline.c_str());

      if(strcmp(line, "")) {
        if(line[0] == '#') sequenNumber++;
      } else break;
    }

    /* ******************************************************************** */
    for(nline.clear(), file.read(&c, 1); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;
    if(file.eof()) break;
    /* ******************************************************************** */
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  sequenNames =  new char*[sequenNumber];
  seqs        = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  do { file.read(&c, 1); } while((c != '#') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;
  /* ******************************************************************** */

  i = 0;
  while(true) {
    if(nline[0] != '!') {
      next = nline.size() + 1; pos  = -1;

      /* ***** ***** We need to do this because it could be some ***** ***** 
         ***** ***** comment in the middle of the sequence ***** *****  */
      while(true) {
        pos  = nline.find("\"", (pos + 1));
        npos = nline.find("\"", (pos + 1));
        next = nline.rfind("\"", next - 1);

        if(pos != (int) string::npos) {
          if(npos == next) {
            nline.erase(pos, (next - pos + 1));
            pos = -1; next = nline.size() + 1;
          }
        } else break;
      }

      do {
        pos = -1; npos = pos;
        while((pos = nline.find("[", (pos + 1))) != (int) string::npos) npos = pos;

        if((int) nline.find("]", (npos + 1)) != (int) string::npos)
          nline.erase(npos, (nline.find("]", (npos + 1)) - npos + 1));
      } while(npos != -1);

      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      state = nline.find("\r", 0);
      while(state != (int) string::npos) {
        nline.erase(state, 1);
        state = nline.find("\r", (state + 1));
      }

      delete [] line;
      line = new char [nline.size() + 1];
      strcpy(line, nline.c_str());
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      str = strtok(line, " #\n");
      if(firstBlock) {
        sequenNames[i] = new char[strlen(str) + 1];
        strcpy(sequenNames[i], str);
      }

      while(1) {
        str = strtok(NULL, " \n");
        if(str == NULL) break;
        else seqs[i].append(str, strlen(str));
      }
      i++;
    }
    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    if(i >= sequenNumber) {
      i = 0; firstBlock = false; 
    }
    /* ******************************************************************** */
  }

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, true);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

bool alignment::loadMegaNonInterleavedAlignment(char *alignmentFile) {

  int i, pos, npos, next, state, firstLine = true;
  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL) 
    freeAlignment();

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  /* ******************************************************************** */
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;

  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  line = new char [nline.size()+1];
  strcpy(line, nline.c_str());
  /* ******************************************************************** */
  str = strtok(line, "!: ");
  for(i = 0; i < (int) strlen(str); i++)
    str[i] = toupper(str[i]);
  /* ******************************************************************** */
  if(!strcmp(str, "TITLE")) {
    delete [] inputFileName;

    if(nline[0] == '!') {
      inputFileName = new char [nline.size() + 1];
      inputFileName[nline.size()] = '\0';
      strcpy(inputFileName, nline.c_str());
    }
    else {
      inputFileName = new char [nline.size() + 2];
      inputFileName[nline.size() + 1] = '\0';
      strcpy(inputFileName, "!");
      strcat(inputFileName, nline.c_str());
    }
    nline.clear();
  }

  /* ******************************************************************** */
  do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  nline.clear(); nline += c;

  while((c != '\n') && (!file.eof())) {
   for(file.read(&c, 1); (c != '\n') && (!file.eof()) && (c != '#'); file.read(&c, 1))
     nline += c;

   if(file.eof() || (c == '#')) break;
   nline += c;
   file.read(&c, 1);
  }
  /* ******************************************************************** */
  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  alignmentInfo = new char [nline.size() + 1];
  strcpy(alignmentInfo, nline.c_str());
  /* ******************************************************************** */
  while((c != '#') && (!file.eof())) { file.read(&c, 1); }
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;
  /* ******************************************************************** */

  while(1) {
    /* ******************************************************************** */
    if(nline[0] != '!') {
      next = nline.size() + 1; pos  = -1;

      /* ***** ***** We need to do this because it could be some ***** ***** 
         ***** ***** comment in the middle of the sequence ***** *****  */
      while(1) {
        pos  = nline.find("\"", (pos + 1));
        npos = nline.find("\"", (pos + 1));
        next = nline.rfind("\"", next - 1);

        if(pos != (int) string::npos) {
          if(npos == next) {
            nline.erase(pos, (next - pos + 1));
            pos = -1; next = nline.size() + 1;
          }
        } else break;
      }
      do {
        pos = -1; npos = pos;
        while((pos = nline.find("[", (pos + 1))) != (int) string::npos) npos = pos;
        if((int) nline.find("]", (npos + 1)) != (int) string::npos)
          nline.erase(npos, (nline.find("]", (npos + 1)) - npos + 1));
      } while(npos != -1);

      /* ******************************************************************** */
      state = nline.find("\r", 0);
      while(state != (int) string::npos) {
        nline.erase(state, 1);
        state = nline.find("\r", (state + 1));
      }

      delete [] line;
      line = new char [nline.size()+1];
      strcpy(line, nline.c_str());
      /* ******************************************************************** */

      if(line[0] == '#') sequenNumber++;
    }

    /* ******************************************************************** */
    for(nline.clear(), file.read(&c, 1); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;
    if(file.eof()) break;
    /* ******************************************************************** */
  }

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  sequenNames =  new char*[sequenNumber];
  seqs        = new string[sequenNumber];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */

  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  do { file.read(&c, 1); } while((c != '\n') && (!file.eof()));
  do { file.read(&c, 1); } while((c != '#') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;
  /* ******************************************************************** */

  i = 0;
  while(true) {

    if(nline[0] != '!') {
      next = nline.size() + 1; pos  = -1;

      /* ***** ***** We need to do this because it could be some ***** ***** 
         ***** ***** comment in the middle of the sequence ***** *****  */
      while(true) {
        pos  = nline.find("\"", (pos + 1));
        npos = nline.find("\"", (pos + 1));
        next = nline.rfind("\"", next - 1);

        if(pos != (int) string::npos) {
          if(npos == next) {
            nline.erase(pos, (next - pos + 1));
            pos = -1; next = nline.size() + 1;
          }
        } else break;
      }

      do {
        pos = -1; npos = pos;
        while((pos = nline.find("[", (pos + 1))) != (int) string::npos) npos = pos;

        if((int) nline.find("]", (npos + 1)) != (int) string::npos)
          nline.erase(npos, (nline.find("]", (npos + 1)) - npos + 1));
      } while(npos != -1);

      /* ***** ***** ***** ***** ***** ***** ***** ***** */
      state = nline.find("\r", 0);
      while(state != (int) string::npos) {
        nline.erase(state, 1);
        state = nline.find("\r", (state + 1));
      }

      delete [] line;
      line = new char [nline.size() + 1];
      strcpy(line, nline.c_str());
      /* ***** ***** ***** ***** ***** ***** ***** ***** */

      str = strtok(line, " #\n");
      if(firstLine) {
        sequenNames[i] = new char[strlen(str) + 1];
        strcpy(sequenNames[i], str);
        firstLine = false;
        str = strtok(NULL, " \n");
      }

      while(1) {
        if(str == NULL) break;
        else seqs[i].append(str, strlen(str));
        str = strtok(NULL, " \n");
      }
    }
    /* ******************************************************************** */
    do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    if(nline[0] == '#') {
      i++; firstLine = true;
    }
    /* ******************************************************************** */
  }

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, true);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

bool alignment::loadNBRF_PirAlignment(char *alignmentFile) {

  int i, state, firstLine = true, seqLines = false;
  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  if(alignmentMatrix != NULL) 
    freeAlignment();
  /* ******************************************************************** */

  do {

    /* ******************************************************************** */
    do{ file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size() + 1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    if((line[0] == '>') && (line[3] == ';')) sequenNumber++;

  } while(1);

  /* Finish to preprocess the input file. */
  file.clear();
  file.seekg(0);

  /* ***** ***** ***** Allocate memory ***** ***** ***** */
  sequenNames =  new char*[sequenNumber];
  seqs        = new string[sequenNumber];

  seqInfo     = new  char*[sequenNumber];
  for(i = 0; i < sequenNumber; i++) seqInfo[i] = new char[3];
  /* ***** ***** ***** ***** ***** ***** ***** ***** */
  i = 0;
  while(1) {
    /* ******************************************************************** */

    if(firstLine) do { file.read(&c, 1); } while((c == '\n') && (!file.eof()));
    else file.read(&c, 1);

    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size() + 1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    if((line[0] == '>') && (line[3] == ';') && (firstLine)) {
      firstLine = false;
      str = strtok(line, ">;"); 
      strcpy(seqInfo[i], str);
      str = strtok(NULL, ">;");
      sequenNames[i] = new char[strlen(str) + 1];
      strcpy(sequenNames[i], str);
    }
    /* ******************************************************************** */
    else if((!firstLine) && (!seqLines)) seqLines = true;
    /* ******************************************************************** */
    else if(seqLines) {

      if(line[strlen(line)-1] == '*') {
        seqLines = false; firstLine = true;
      }

      str = strtok(line, OTHDELIMITERS);
      while(1) {
        if(str == NULL) break;
        else seqs[i].append(str, strlen(str));
        str = strtok(NULL, OTHDELIMITERS);
      }
    }
    /* ******************************************************************** */
    if((firstLine) && (!seqLines)) i++;
    /* ******************************************************************** */
  }
  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, false);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

bool alignment::loadPhylip3_2Alignment(char *alignmentFile) {

  int i, state, firstLine = true;
  char c, *str, *line = NULL;
  string nline, *seqs;
  ifstream file;

  /* ******************************************************************** */
  /* Open alignment file and open check */
  file.open(alignmentFile, ifstream::in);
  if(!file) return false;

  if(alignmentMatrix != NULL)
    freeAlignment();

  inputFileName = new char[strlen(alignmentFile) + 9];
  inputFileName[strlen(alignmentFile) + 8] = '\0';
  strcpy(inputFileName, "!Title ");
  strcat(inputFileName, alignmentFile);
  strcat(inputFileName, ";");

  /* ******************************************************************** */
  do{ file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
    nline += c;

  state = nline.find("\r", 0);
  while(state != (int) string::npos) {
    nline.erase(state, 1);
    state = nline.find("\r", (state + 1));
  }

  delete [] line; 
  line = new char [nline.size()+1];
  strcpy(line, nline.c_str());
  /* ******************************************************************** */

  /* ******************************************************************** */
  str = strtok(line, DELIMITERS);
  if(str != NULL) sequenNumber = atoi(str);
  else return false;

  str = strtok(NULL, DELIMITERS);
  if(str != NULL) residNumber = atoi(str);
  else return false;

  /* Check correct read parameters */
  if((sequenNumber == 0) || (residNumber == 0)) 
    return false;

  /* Allocate memory for alignment matrix and sequenNumber names vector */
  sequenNames = new  char*[sequenNumber];
  seqs        = new string[sequenNumber];

  /* ******************************************************************** */
  do{ file.read(&c, 1); } while((c == '\n') && (!file.eof()));
  i = 0;
/* ******************************************************************** */

  do {
    /* ******************************************************************** */
    for(nline.clear(); (c != '\n') && ((!file.eof())); file.read(&c, 1))
      nline += c;

    if(file.eof()) break;
    delete [] line;

    state = nline.find("\r", 0);
    while(state != (int) string::npos) {
      nline.erase(state, 1);
      state = nline.find("\r", (state + 1));
    }

    line = new char [nline.size() + 1];
    strcpy(line, nline.c_str());
    /* ******************************************************************** */

    str = strtok(line, OTHDELIMITERS);
    if(firstLine) {
      firstLine = false;
      sequenNames[i] = new char[strlen(str) + 1];
      strcpy(sequenNames[i], str);
      str = strtok(NULL, OTHDELIMITERS);
    }

    while(1) {
      if(str == NULL) break;
      else seqs[i].append(str, strlen(str));
      str = strtok(NULL, OTHDELIMITERS);
    }

    file.read(&c, 1); 
    if(c == '\n') {
      i++; firstLine = true;
    }
    while((c == '\n') && (!file.eof())){ file.read(&c, 1); }
    /* ******************************************************************** */
  } while(!file.eof());

  /* Close the alignment file */
  file.close();

  /* ***** ***** ***** Fill our data structure.  ***** ***** ***** */
  state = fillMatrices(seqs, true);

  /* ***** ***** ***** Free the dinamic memory.  ***** ***** ***** */
  delete [] line;
  delete [] seqs;

  /* ***** ***** ***** Return the load operation state  ***** ***** ***** */
  return state;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
| void alignment::alignmentToFile(ostream &file)	|
|  	Private method that put the alignment on the    |
|  	parameter stream in PHYLIP2               |
+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void alignment::alignmentPhylipToFile(ostream &file) {
  
  int i, j, maxLongName = 0;
  char str[61], *name;
  
  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not have the sequences aligned." << endl << endl;
    return;
  }

  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;
  
  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, strlen(sequenNames[i])); 
  
  /* It's necessary to do this in order to have not problems with phylip's programs */
  maxLongName = maxLongName > 10 ? 10 : maxLongName;

  /* Create a string with a maxLongName's lenght and concatenate '\0' at the end of the string. */
  name = new char[maxLongName+1];
  name[maxLongName] = '\0';

  /* Concatenate '\0' to the final of the str string */
  str[60] = '\0';
  
  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
    strncpy(name, sequenNames[i], maxLongName);
    file << setw(PHYLIPDISTANCE + 3) << left << name;
    strncpy(str, alignmentMatrix[i], 60);
    file << str << endl;
  }
  file << endl;
  
  /* Put on the stream the rest of the blocks */
  for(i = 60; i < residNumber; i += 60) {
    for(j = 0; j < sequenNumber; j++) {
      strncpy(str, &alignmentMatrix[j][i], 60);
      file << str << endl;
    }
    file << endl;
  }
  file << endl;

  delete [] name;
}

void alignment::alignmentClustalToFile(ostream &file) {

  int i, j, maxLongName = 0;
  char str[61]; str[60] = '\0';

  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not have the sequences aligned." << endl << endl;
    return;
  }

  if((alignmentInfo != NULL)  && (iformat == oformat))
    file << alignmentInfo << endl << endl;
  else 
    file << "CLUSTAL W (1.8) multiple sequence alignment" << endl << endl;

  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, strlen(sequenNames[i])); 

  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(j = 0; j < residNumber; j += 60) {
    file << endl;
    for(i = 0; i < sequenNumber; i++) {
      file << setw(maxLongName+5) << left << sequenNames[i];
      strncpy(str, &alignmentMatrix[i][j], 60); 
      file << str << endl;
    }
    file << endl;
  }
}


void alignment::alignmentNBRF_PirToFile(ostream &file) {

  int i, j, k; 
  char str[11];

  str[10] = '\0';

  for(i = 0; i < sequenNumber; i++) {
    file << ">";

    if((seqInfo != NULL) && (iformat == oformat)) file << seqInfo[i];
    else {
      getTypeAlignment();
      switch(dataType) {
        case DNAType: file << "DL"; break;
        case RNAType: file << "RL"; break;
        case AAType:  file << "P1"; break;
      }
    }

    file << ";" << sequenNames[i] << endl;
    file << sequenNames[i] << " " << residuesNumber[i] << " bases" << endl;

    for(j = 0; j < residuesNumber[i]; j += 50) {
      for(k = j; k < residuesNumber[i] && k < j + 50; k += 10) {
        strncpy(str, &alignmentMatrix[i][k], 10); 
        file << " " << str;
      }

      if((j + 50) >= residNumber)
        file << "*"; 
      file << endl;
    }
    file << endl;
  }
}

void alignment::alignmentFastaToFile(ostream &file) {

  int i, j;
  char str[61];

  str[60] = '\0';

  for(i = 0; i < sequenNumber; i++) {
    file << ">" << sequenNames[i] << " " << dec << residuesNumber[i] << " bp" << endl;

    for(j = 0; j < residuesNumber[i]; j += 60) {
      strncpy(str, &alignmentMatrix[i][j], 60);
      file << str << endl;
    }
  }

}

void alignment::alignmentPhylip3_2ToFile(ostream &file) {

  int i, j, k, maxLongName = 0;
  char str[11];

  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not have the sequences aligned." << endl << endl;
    return;
  }

  /* Include in the first line the sequenNumber and the aminoacids of the alignment */
  file << " " << sequenNumber << " " << residNumber << endl;

  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, strlen(sequenNames[i]));

  /* Concatenate '\0' to the final of the str string */
  str[10] = '\0';

  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
    file << setw(maxLongName+3) << left << sequenNames[i];

    for(j = 0; j < residNumber; j += 50) {
      for(k = j; k < residNumber && k < j + 50; k += 10) {
        strncpy(str, &alignmentMatrix[i][k], 10); 
        file << str;

        if(k+10 < residNumber && k+10 < j + 50)
          file << " ";
      }

      file << endl;
      if(j + 50 < residNumber) 
        file << setw(maxLongName+3) << " ";
    }

    file << endl;
  }
}

void alignment::alignmentNexusToFile(ostream &file) {

  int i, j, k, maxLongName = 0;
  char line[11], *str = NULL;

  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not have the sequences aligned." << endl << endl;
    return;
  }

  line[10] = '\0';
  file << "#NEXUS" << endl;

  if((alignmentInfo != NULL) && (iformat == oformat)) {
    str = strtok(alignmentInfo, "|");

    if(strcmp(str, "\n"))
      file << str << endl << endl;
    else file << endl << endl;

    str = strtok(NULL, "|");

    for(i = 0; i < (int) strlen(str); i++) 
      str[i] = toupper(str[i]);

  } else  file << endl << endl;

  file << "BEGIN DATA;" << endl;
  file << " DIMENSIONS NTAX=" << sequenNumber << " NCHAR=" << residNumber <<";" << endl;

  if((str != NULL) && (iformat == oformat)) file << str << endl;

  else {
    getTypeAlignment();

    file << " FORMAT";
    switch(dataType) {
      case DNAType: file << " DATATYPE=DNA "; break;
      case RNAType: file << " DATATYPE=RNA "; break;
      case AAType:  file << " DATATYPE=PROTEIN "; break;
    }
    file << "INTERLEAVE=yes GAP=-;" << endl;
  }

  for(i = 0; i < sequenNumber; i++)
    file << "[Name: " << sequenNames[i] << "\t\t\tLen: " << residNumber << " Check: 0]" << endl;

  file << endl << "MATRIX" << endl;

  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, strlen(sequenNames[i])); 

  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(j = 0; j < residNumber; j += 50) {
    for(i = 0; i < sequenNumber; i++) {
      file << setw(maxLongName+4) << left << sequenNames[i];
      for(k = j; k < (j + 50) && k < residNumber; k += 10) {
        strncpy(line, &alignmentMatrix[i][k], 10); file << " " << line;
      }
      file << endl;
    }
    file << endl;
  }
  file << endl << ";" << endl << "END;" << endl;
}

void alignment::alignmentMegaToFile(ostream &file) {

  int i, j, k;
  char str[11];

  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not have the sequences aligned." << endl << endl;
    return;
  }
  str[10] = '\0';

  /* Phylemon Webserver Version */
  /* file << "#MEGA" << endl << "Alignment file" << endl; */

  /* Standard Version */
  file << "#MEGA" << endl << inputFileName << endl;

  getTypeAlignment();
  switch(dataType) {
    case DNAType:
      file << "!Format DataType=DNA ";
      break;
    case RNAType:
      file << "!Format DataType=RNA ";
      break;
    case AAType:
      file << "!Format DataType=protein ";
      break;
  }
  file << "NSeqs=" << sequenNumber << " Nsites=" << residNumber << " indel=- CodeTable=Standard;" << endl << endl;

  for(i = 0; i < sequenNumber; i++) {
    file << "#" << sequenNames[i] << endl;

    for(j = 0; j < residNumber; j += 50) {
      for(k = j; ((k < residNumber) && (k < j + 50)); k += 10) {
        strncpy(str, &alignmentMatrix[i][k], 10);
        file << str;

        if(((k + 10) < residNumber) && ((k + 10) < (j + 50)))
          file << " ";
        else
          file << endl;
      }
    }
    file << endl;
  }
}

int alignment::alignmentToFile(ostream &file, bool header, int init, int lenght, int tamNames) {
  int i, maxLongName = 0;
  char *str;

  if(!isAligned) {
    if(header) file << "ERROR: Imposible to generate header. The input file does not have the sequences "
                    << "aligned." << endl;
 } else
    if(header) file << sequenNumber << " " << residNumber << endl;

  for(i = 0; i < sequenNumber; i++) maxLongName = utils::max(maxLongName, strlen(sequenNames[i])); 
  maxLongName = (maxLongName > tamNames) ? maxLongName : tamNames;

  str = new char[lenght+1];
  str[lenght] = '\0'; 

  /* Put on the stream the sequenNumber names and the first 60 aminoacids block */
  for(i = 0; i < sequenNumber; i++) {
    file << setw(maxLongName+2) << left << sequenNames[i];
    strncpy(str, &alignmentMatrix[i][init], lenght);
    file << str << endl;
  }

  file << endl;

  delete [] str;
  return maxLongName;
}

bool alignment::alignmentSummaryHTML(char *destFile, int residues, int sequences, int *selectedRes, int *selectedSeq) {

  int i, j, k, maxLongName = 0;
  bool *res, *seq;
  ofstream file;
  char *seqname;

  if(!isAligned) {
    cerr << endl << "ERROR: The input file does not have the sequences aligned." << endl << endl;
    return false;
  }

  /* Open file and check the operations */
  file.open(destFile);
  if(!file) return false;

  /* Look for the maximum sequenNumber name size */
  for(i = 0; i < sequenNumber; i++)
    maxLongName = utils::max(maxLongName, strlen(sequenNames[i])); 
  seqname = new char [(maxLongName + 10)];

  res = new bool[residNumber];
  for(i = 0; i < residNumber; i++) res[i] = false;
  for(i = 0; i < residues; i++) res[selectedRes[i]] = true;

  seq = new bool[sequenNumber];
  for(i = 0; i < sequenNumber; i++) seq[i] = false;
  for(i = 0; i < sequences; i++) seq[selectedSeq[i]] = true;

  file << "<!DOCTYPE html>" << endl << "<html><head>" << endl;
  file << "    <meta http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />" << endl;
  file << "    <title>trimAl v1.2 Summary</title>" << endl;

  file << "    <style type=\"text/css\">" << endl;
  file << "    .sel  { background-color: #C9C9C9; }\n";
  file << "    </style>\n  </head>\n\n";
  file << "  <body>\n";

  file << "\t<pre>" << endl;
  file << "  <span class=sel>Selected Residue / Sequence</span>" << endl;
  file << "  Deleted Residue / Sequence";

  for(j = 0; j < residNumber; j += 60) {

    file << endl;
    for(i = j + 1; ((i <= residNumber) && (i <= (j + 10))); i++)
      if(!((i + 1) % 10)) file << setw(maxLongName + 19) << right << (i + 1);
    for(i = j + 11; ((i <= residNumber) && (i <= (j + 60))); i++)
      if(!((i + 1) % 10)) file << setw(10) << right << (i + 1);

    file << endl << setw(maxLongName + 10);
    for(i = j + 1; ((i <= residNumber) && (i <= (j + 60))); i++)
    if(!(i % 10)) file << "+";
    else file << "=";

    file << endl;
    for(i = 0; i < sequenNumber; i++) {

      strcpy(seqname, sequenNames[i]); strcat(seqname, "</span>");
      if(seq[i]) file << "    <span class=sel>" << setw(maxLongName + 12) << left << seqname;
      else       file << "    <span>" << setw(maxLongName + 12) << left << seqname;

      for(k = j; ((k < residNumber) && (k < (j + 60))); k++) {
        if((seq[i]) && (res[k])) file << "<span class=sel>" << alignmentMatrix[i][k] << "</span>";
        else                     file << alignmentMatrix[i][k];
      }
      file << endl;
    }
  }
  file << "\t  </pre>" << endl;
  file << "  </body>" << endl << "</html>" << endl;

  file.close();
  delete [] seq;
  delete [] res;
  delete [] seqname;

  return true;
}

