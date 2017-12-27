//-------------------------------------------------------------------------/
//Copyright (C) 2003, 2004, 2005, ALL RIGHTS RESERVED.
//Centre for Sys. Eng. & App. Mech.           FEMAGSoft S.A.
//Universite Cathalique de Louvain            4, Avenue Albert Einstein
//Batiment Euler, Avenue Georges Lemaitre, 4  B-1348 Louvain-la-Neuve
//B-1348, Louvain-la-Neuve                    Belgium
//Belgium
//-------------------------------------------------------------------------/
//
//Name:         geometry.h (all geometry premitives related polygon triang-
//              ulation by sweep line algorithm)
//Author:       Liang, Wu (wu@mema.ucl.ac.be, wuliang@femagsoft.com)
//Created:      03/2001
//Modified:     10/2005. Modified and simplified only for polygon triangul-
//              ation purpose.
//-------------------------------------------------------------------------/

#ifndef POLYTRI_GEOMETRY_H
#define POLYTRI_GEOMETRY_H

//##includes now in defs

#include "PolyTri_Defs.h"
#include "PolyTri_Splay.h"
#include "../../TopologyTracer2D3D_Datatypes.h"


//base class for points;
class Pointbase
{
public:
	Pointbase() {}
	Pointbase(const Pointbase& pb);
	Pointbase(double xx, double yy) :id(0), x(xx), y(yy), type(UNKNOWN) {}
	Pointbase(int idd, double xx, double yy) :id(idd), x(xx), y(yy), type(UNKNOWN) {}
	Pointbase(double xx, double yy, Type ttype) :id(0), x(xx), y(yy), type(ttype) { }
	Pointbase(int idd, double xx, double yy, Type ttype) :id(idd),x(xx), y(yy), type(ttype) { }

	//operator overloading
	friend  bool operator==(const Pointbase&, const Pointbase&); 
	friend  bool operator>(const Pointbase&, const Pointbase&);
	friend  bool operator<(const Pointbase&, const Pointbase&);
	friend  bool operator!=(const Pointbase&, const Pointbase&);
	friend  ostream &operator<<(ostream &os, const Pointbase& point);

	//public data
	unsigned int	id;				//id of point;
	double			x, y;			//coordinates;
	Type			type;			//type of points;
	bool			left;			//left chain or not;
};


//base class for polygon boundary
//Linebase class is a directed line segment with start/end point
class Linebase
{
	public:
	//constructors and destructor
	Linebase();
	Linebase(Pointbase* ep1, Pointbase* ep2, Type type, unsigned int idd);
	Linebase(const Linebase& line);
	~Linebase() {};

	unsigned int id() const { return _id; }

	//two end points 
	Pointbase*   endPoint(int i) const { return _endp[i]; }
	Type         type() const { return _type; }
	double       keyValue() const { return _key; }
	void         setKeyValue(double y); 
	//slightly increased the key to avoid duplicated key for searching tree. 
	void         increaseKeyValue(const double diff) { _key+=diff; }
	//reverse a directed line segment; reversable only for inserted diagonals
	void         reverse(); 

	//set and return helper of a directed line segment;
	void         setHelper(unsigned int i) { _helper=i; }
	unsigned int helper() { return _helper; }

	//operator overloading             
	friend ostream &operator<<(ostream &os, const Linebase& line);

	protected:
	unsigned int _id;           //id of a line segment;
	Pointbase*   _endp[2];      //two end points, ##MK::points memory freed elsewhere???

	Type         _type;         //type of a line segement, input/insert
	double       _key;          //key of a line segment for splay tree searching
	unsigned int _helper;       //helper of a line segemnt
};


class Polygon
{
public:
	Polygon(vector<point2d>* points);
	Polygon(vector<double>* points);
	~Polygon();

	//##MK::public?
	PointbaseMap& points() { return _points; }
	LineMap& edges() { return _edges; }

	// main member function for polygon triangulation
	bool partition2Monotone();
	bool searchMonotones();
	bool triangulation();
	void result();
	vector<triangle>* giveTriangles();
	void verboseTriangles();

	//return all triangles
	Triangles triangles() { return _triangles; }

	unsigned int get_lid() { return l_id; }
	void incr_lid() { l_id = l_id + 1; }
	void set_lid(unsigned int i ) { l_id = i; }
	bool is_healthy(){ return myhealth; }

private:
	void initialize();
	bool does_point_exist_already(double x, double y);

	//prev or next point/edge id for a given ith point/edge
	unsigned int prev(unsigned int i);
	unsigned int next(unsigned int i);

	//handle event vertext according to vertex type
	void handleStartVertex(unsigned int);
	void handleEndVertex(unsigned int);
	void handleSplitVertex(unsigned int);
	void handleMergeVertex(unsigned int);
	void handleRegularVertexUp(unsigned int);
	void handleRegularVertexDown(unsigned int);

	//add diagonal between two vertices
	void addDiagonal(unsigned int i, unsigned int j);

	//angle ABC for three given points, for monotone polygon searching purpose
	double angleCosb(double *A, double *B, double *C);
	//find the next edge, for monotone polygon searching purpose
	unsigned int selectNextEdge(Linebase* edge);

	//triangulate a monotone polygon piece
	bool triangulateMonotone(Monopoly& mpoly);

//private data memers
	PQueue      			_qpoints;				//priority queue for event points
	EdgeBST     			_edgebst;				//edge binary searching tree (splaytree)
	Monopolys   			_mpolys;				//all monotone polygon piece list;
	Triangles   			_triangles;				//all triangle list;

	//data for monotone piece searching purpose;
	AdjEdgeMap				_startAdjEdgeMap;		//all edges starting from given points (map)
	LineMap					_diagonals;				//added diagonals to partition polygon to
													//monotone pieces, not all diagonals of given poly
	vector<unsigned int>	_nVertices;
	PointbaseMap			_points;				//all vertices
	LineMap					_edges;					//all edges 
	double					xmin;
	double					xmax;
	double					ymin;
	double					ymax;					//AABB for polygon
	unsigned int			l_id;
	bool					myhealth;
};

#endif