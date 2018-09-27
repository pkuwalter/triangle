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
#include "sequence.h"
#include "graph.h"
#include "utils.h"
#include "parallel.h"
#include "sampleSort.h"
using namespace std;

struct isSameColor{
  isSameColor(uint mycolor, uint *colors) : colors_(colors), me_(mycolor) {};
  bool operator () (uint v) {return me_ == colors_[v];};
  uint me_;
  uint *colors_;
};

// **************************************************************
// * using the degree heuristic to order the vertices
// **************************************************************

struct nodeLT {
  nodeLT(uintT* Degrees_) : Degrees(Degrees_) {};
  bool operator() (uint a, uint b) {
    return Degrees[a] < Degrees[b];
  };  
  uintT* Degrees;
};

void rankNodes(uintT* Degrees, uint *r, uint *o, uintT n)
{
  parallel_for (long i=0;i<n;i++) o[i] = i;
  compSort(o, n, nodeLT(Degrees));
  parallel_for (long i=0;i<n;i++) r[o[i]] = i;  
}

struct isFwd{
  isFwd(uint myrank, uint *r) : r_(r), me_(myrank) {};
  bool operator () (uint v) {return r_[v] > me_;};
  uint me_;
  uint *r_;
};

struct intLT {
  bool operator () (uintT a, uintT b) { return a < b; };
};

long countCommon(uint *A, uintT nA, uint *B, uintT nB)
{
  uintT i=0,j=0;
  long ans=0;
  while (i < nA && j < nB) {
    if (A[i]==B[j]) i++, j++, ans++;
    else if (A[i] < B[j]) i++;
    else j++;
  }
  return ans;
}

struct countFromA {  
  countFromA (uint *_e, uintT *_start)
    : edges(_e), start(_start) {} ;
  long operator() (uintT _a) {
    long count = 0;
    uintT a = _a;
    uintT tw = 0;
    uintT sza = start[a+1]-start[a];
    for (uintT bi=0;bi<sza;bi++) {
      uintT b=edges[start[a]+bi];
      uintT szb = start[b+1]-start[b];
      tw += min<uintT>(sza, szb);
    }
    if (tw > 10000) {
      long *cc = newA(long, sza);
      parallel_for (uintT bi=0;bi<sza;bi++) {
        uintT b=edges[start[a]+bi];
	uintT szb = start[b+1]-start[b];
        cc[bi] = countCommon(edges + start[a], sza,
            edges + start[b], szb);
      }
      count = sequence::plusReduce(cc, sza);
      free(cc);
    }
    else {
      for (uintT bi=0;bi<sza;bi++) {
	uintT b=edges[start[a]+bi];
	uintT szb = start[b+1]-start[b];
	count += countCommon(edges + start[a], sza,
			     edges + start[b], szb);
      }
    }  
   return count;
  }
  uint* edges;
  uintT *start;
};

long countTriangle(graphC<uintT,uint> GG, double p, long seed)
{
  long n = GG.n, m = GG.m;
  //sample
  uintT numColors = max<uintT>(1,1/p);
  uint* colors = newA(uint,n);
  //assign colors
  parallel_for(long i=0;i<n;i++) colors[i] = utils::hash(seed+i) % numColors;

  uint* edgesSparse = newA(uint,m);
  //keep only neighbors matching own color
  uintT* Degrees = newA(uintT,n+1); Degrees[n] = 0;
  parallel_for(long i=0;i<n;i++) {
    uintT start = GG.offsets[i];
    uintT d = GG.offsets[i+1]-start;
    uint color = colors[i];
    if(d > 10000) {
      Degrees[i] = sequence::filter(GG.edges+start, edgesSparse+start, 
    			       d, isSameColor(color, colors));
    } else {
      uintT k = 0;
      uint* newNghs = edgesSparse + start;
      uint* Nghs = GG.edges + start;
      for(uintT j=0;j<d;j++) {
	uintT ngh = Nghs[j];
	if(colors[ngh] == color) newNghs[k++] = ngh;
      }
      Degrees[i] = k;
    }
  }
  free(colors);

  uint *rank = newA(uint, n);
  uint *order = newA(uint, n);

  rankNodes(Degrees, rank, order, n);

  free(order);

  long new_m = sequence::plusScan(Degrees,Degrees,n+1);
  m = new_m;

  uintT *sz = newA(uintT,n+1); sz[n] = 0;

  //creating subgraph and filtering non-forward edges at the same time
  uint* new_edges = newA(uint,m);
  parallel_for(long i=0;i<n;i++) {
    uintT start = GG.offsets[i];
    uint* Nghs = edgesSparse + start;
    uintT o = Degrees[i];
    uintT d = Degrees[i+1]-o;

    if(d > 10000) {
      sz[i] = sequence::filter(Nghs, new_edges+o,
    			       d, isFwd(rank[i], rank));
    } else {
      uintT k=0;
      uint* newNgh = new_edges+o;
      for(uintT j=0;j<d;j++) { 
	uintT ngh = Nghs[j];
	if(rank[i] < rank[ngh]) newNgh[k++] = ngh;
      }
      sz[i] = k;
    }
    compSort(new_edges+o, sz[i], intLT());
  }
  free(edgesSparse);
  free(rank); 

  //pack edges down to array of half the size
  uint* edges2 = newA(uint,m/2);
  uintT* sz2 = newA(uintT,n+1); 
  sequence::plusScan(sz,sz2,n+1);

  parallel_for(long i=0;i<n;i++) {
    uintT o = Degrees[i];
    uintT d = sz[i];
    uintT start = sz2[i];
    if(d > 10000)
      parallel_for(long j=0;j<d;j++) edges2[start+j] = new_edges[o+j];
    else 
      for(long j=0;j<d;j++) edges2[start+j] = new_edges[o+j];
  }
  free(sz); free(new_edges); free(Degrees);
  new_edges = edges2;

  // start counting
  long count = 
    sequence::reduce<long>((long) 0, (long) n, 
			   utils::addF<long>(),
			   countFromA(new_edges,sz2));
  free(new_edges); free(sz2); 
  long triCount = (p == 1.0) ? count : (long) ((double)count/(p*p));
  cout << "tri. count = " << triCount << endl;
  return triCount;
}
