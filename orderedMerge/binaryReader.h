// This code is implemented as part of the paper "Multicore Triangle
// Computations Without Tuning" by Julian Shun and Kanat Tangwongsan in
// Proceedings of the IEEE International Conference on Data Engineering
// (ICDE), 2015.
//
// Copyright (c) 2015 Julian Shun and Kanat Tangwongsan
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef _BINARY_READER_
#define _BINARY_READER_

#include <iostream>
#include "gettime.h"
#include "utils.h"
#include "graph.h"
#include "parallel.h"
#include "IO.h"
#include "graphIO.h"
#include "binaryReader.h"

graphC<uintT,uint> readGraphCFromBinary(char* iFile) {
  char* config = (char*) ".config";
  char* adj = (char*) ".adj";
  char* idx = (char*) ".idx";
  char configFile[strlen(iFile)+7];
  char adjFile[strlen(iFile)+4];
  char idxFile[strlen(iFile)+4];
  strcpy(configFile,iFile);
  strcpy(adjFile,iFile);
  strcpy(idxFile,iFile);
  strcat(configFile,config);
  strcat(adjFile,adj);
  strcat(idxFile,idx);

  ifstream in(configFile, ifstream::in);
  long n;

  in >> n;
  in.close();

  ifstream in2(adjFile,ifstream::in | ios::binary); //stored as uints
  in2.seekg(0, ios::end);
  long size = in2.tellg();
  in2.seekg(0);
  long m = size/sizeof(uint);
  //cout << "adj size = " << size << endl;
  char* s = (char *) malloc(size);
  in2.read(s,size);
  in2.close();
  //cout << "n = "<<n<<" m = "<<m << endl;
  uint* edges = (uint*) s;
  ifstream in3(idxFile,ifstream::in | ios::binary); //stored as longs
  in3.seekg(0, ios::end);
  size = in3.tellg();
  in3.seekg(0);
  if(n != size/sizeof(intT)) { cout << "File size wrong\n"; abort(); }

  //make room for one more intT
  char* t = (char *) malloc(size + sizeof(intT)); 
  in3.read(t,size);
  in3.close();
  uintT* offsets = (uintT*) t;
  offsets[n] = m; //use last int in offsets array for edge case
  return graphC<uintT,uint>(offsets,edges,n,m);
}

#endif
