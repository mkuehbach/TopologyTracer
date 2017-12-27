//-------------------------------------------------------------------------/
//Copyright (C) 2003, 2004, 2005, ALL RIGHTS RESERVED.
//Centre for Sys. Eng. & App. Mech.           FEMAGSoft S.A.
//Universite Cathalique de Louvain            4, Avenue Albert Einstein
//Batiment Euler, Avenue Georges Lemaitre, 4  B-1348 Louvain-la-Neuve
//B-1348, Louvain-la-Neuve                    Belgium
//Belgium
//-------------------------------------------------------------------------/
//Name:         defs.h
//Purpose:      some definition for mesh generation
//Author:       Wu Liang
//Created:      03/2000
//-------------------------------------------------------------------------/

#ifndef POLYTRI_DEFS_H
#define POLYTRI_DEFS_H

#include <math.h>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <queue>
#include <set>
#include <stdlib.h>

#include <iomanip>
#include <fstream>
#include <limits>

#include <stack>


using namespace std;


typedef double precision;
struct tri {
	precision p1x;		//points are ordered clockwise!
	precision p1y;
	precision p2x;
	precision p2y;
	precision p3x;
	precision p3y;
	tri() : p1x(0.0), p1y(0.0), p2x(0.0), p2y(0.0), p3x(0.0), p3y(0.0) {}
};

#define sqr(t)  ((t)*(t))

//#define INFINITY	(1.0e38);
//#define PI			(3.141592653589793238462643383279502884197169399375105820974944592308)
enum  Type			{ UNKNOWN, INPUT, INSERT, START, END, MERGE, SPLIT, REGULAR_UP, REGULAR_DOWN};

class   Pointbase;
class   Linebase;

template <class T, class KeyType>      class    SplayTree;
typedef map<unsigned int, Pointbase*>           PointbaseMap;
typedef map<unsigned int, Linebase*>            LineMap;
typedef priority_queue<Pointbase>               PQueue;
typedef SplayTree<Linebase*, double>            EdgeBST;
typedef list<unsigned int>                      Monopoly;
typedef list<Monopoly>                          Monopolys;
typedef vector<unsigned int>                    Triangle;
typedef list<Triangle>                          Triangles;
typedef map<unsigned int, set<unsigned int> >   AdjEdgeMap;
#endif
