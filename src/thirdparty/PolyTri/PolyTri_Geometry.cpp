//-------------------------------------------------------------------------/
//Copyright (C) 2003, 2004, 2005, ALL RIGHTS RESERVED.
//Centre for Sys. Eng. & App. Mech.           FEMAGSoft S.A.
//Universite Cathalique de Louvain            4, Avenue Albert Einstein
//Batiment Euler, Avenue Georges Lemaitre, 4  B-1348 Louvain-la-Neuve
//B-1348, Louvain-la-Neuve                    Belgium
//Belgium
//-------------------------------------------------------------------------/
//
//Name:         geometry.cc (all geometry premitive implementations related 
//              to polygon triangulation by sweep line algorithm)
//Author:       Liang, Wu (wu@mema.ucl.ac.be, wuliang@femagsoft.com)
//Created:      03/2001
//Modified:     10/2005. Modified and simplified only for polygon triangul-
//              ation purpose.
//-------------------------------------------------------------------------/

//#include <stack>
//#include <limits>
//#include "thirdparty/PolyTri/
#include "PolyTri_Geometry.h"


//Jonathan Schewchuk's exact arithmetic code, see predicates.cc for details; 
extern double orient2d(double* pa, double* pb, double* pc);

//----------------------------------------------------------------------------
//square of the distance of two points;
//----------------------------------------------------------------------------
double dist_sqr(const Pointbase& sp, const Pointbase& ep)
{
	return sqr(sp.x-ep.x)+sqr(sp.y-ep.y);
}

//----------------------------------------------------------------------------
//square of the distance of two points;
//----------------------------------------------------------------------------
double dist_sqr(double *pa, double *pb)
{
	return sqr(pa[0]-pb[0])+sqr(pa[1]-pb[1]);
}

void UpdateKey(BTreeNode<Linebase*,double>* node, double y)
{
     node->data()->setKeyValue(y);
}

//----------------------------------------------------------------------------
//copy constructor
//----------------------------------------------------------------------------
Pointbase::Pointbase(const Pointbase& pb)
{
	this->id=pb.id; 
	this->x=pb.x;
	this->y=pb.y;
	this->type=pb.type;
	this->left=pb.left;
}

//----------------------------------------------------------------------------
//operator ( ==, >, < and != ) overloading for pointbase class
//----------------------------------------------------------------------------
bool operator==(const Pointbase& pa, const Pointbase& pb)
{
   	return (pa.x==pb.x && pa.y==pb.y);
}

//----------------------------------------------------------------------------
bool operator>(const Pointbase& pa, const Pointbase& pb)
{
   	return( (pa.y > pb.y) || ( (pa.y==pb.y) && (pa.x < pb.x)) );
}

//----------------------------------------------------------------------------
bool operator<(const Pointbase& pa, const Pointbase& pb)
{
	return( (pa.y < pb.y) || ( (pa.y==pb.y) && (pa.x > pb.x)) );
}

//----------------------------------------------------------------------------
bool operator!=(const Pointbase& pa, const Pointbase& pb)
{
   	return !(pa.x==pb.x && pa.y==pb.y);
}

//----------------------------------------------------------------------------
//operator for debugging
//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const  Pointbase& point)
{
    	os<<point.id<<" "<<setw(35)<<point.x<<setw(35)<<point.y<<'\n';
    	return os;
}

//----------------------------------------------------------------------------
//operator for debugging
//----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const Linebase &line)
{
	os<< "Linebase:(" << line._id << ")" << '\n';
	os<< *(line._endp[0]) << *(line._endp[1]);
	os<< " Type=" << line._type << " key value:" << line.keyValue() << '\n';
	return os;
}

//----------------------------------------------------------------------------
//Linebase construct
//----------------------------------------------------------------------------
Linebase::Linebase():_type(UNKNOWN), _id(0)
{
	_endp[0] = NULL;
	_endp[1] = NULL;
}

//-----------------------------------------------------------------------------
//Linebase construct
//-----------------------------------------------------------------------------
Linebase::Linebase(Pointbase* sp, Pointbase* ep, Type type, unsigned int idd) :_type(type), _id(idd)
{
	_endp[0] = sp;
	_endp[1] = ep;
	//_key=_endp[0]->x < _endp[1]->x ? _endp[0]->x:_endp[1]->x;
}

//----------------------------------------------------------------------------
//copy constructor
//----------------------------------------------------------------------------
Linebase::Linebase(const Linebase& line)
{
	this->_id=line._id;
	this->_endp[0]=line._endp[0];
	this->_endp[1]=line._endp[1];
	this->_key=line._key;
	this->_helper=line._helper;
}


//----------------------------------------------------------------------------
//reverse a directed line segment, reverseable only for insert diagonals
//----------------------------------------------------------------------------
void Linebase::reverse() 
{ 
	assert(_type==INSERT); 
	Pointbase* tmp=_endp[0]; 
	_endp[0]=_endp[1];
	_endp[1]=tmp;
}

void Linebase::setKeyValue(double y)
{
        if( _endp[1]->y==_endp[0]->y )  
		_key=_endp[0]->x < _endp[1]->x ? _endp[0]->x:_endp[1]->x;
	else    _key=( y - _endp[0]->y ) * ( _endp[1]->x - _endp[0]->x ) / (_endp[1]->y - _endp[0]->y ) + _endp[0]->x; 
}






//----------------------------------------------------------------------------
//polygon class constructor
//----------------------------------------------------------------------------
Polygon::Polygon( vector<point2d>* points )
{
	//set up variables
	//vector and maps deallocate themselves...
	myhealth = false;
	xmin = std::numeric_limits<double>:: max();
	xmax = std::numeric_limits<double>:: lowest();
	ymin = std::numeric_limits<double>:: max();
	ymax = std::numeric_limits<double>:: lowest();
	this->set_lid(0);

	double x,y;
	Type type;
	type = INPUT;

	//unsigned int sz = points->size();
	_nVertices.push_back( points->size() ) ; // - 1; //assuming closed contour without holes!
	//unsigned int first = 1;
	//unsigned int last = first + _nVertices.at(0) - 1;

	unsigned int i = 1;
	for (unsigned int j = 0; j < _nVertices.at(0); j++, i++) {
		x = points->at(j).x;
		y = points->at(j).y;

		if( does_point_exist_already(x, y) ) {
			cout << "Error! I found the duplicate point " << j << "-->" << x << ";" << y << endl;
			return; //exit(1);
		} 
		else {
			Pointbase* point = new Pointbase( i, x, y, type );
			if(x > xmax) xmax=x;
			if(x < xmin) xmin=x;
			if(y > ymax) ymax=y;
			if(y < ymin) ymin=y;
			_points[i] = point; //polyTri counts point indices on interval [1,size()]
		}
	}

	unsigned int sid;
	unsigned int eid;
	type = INPUT;
	//run over all points connect to lines starting at sid ending at eid
	for( i = 1; i <= _nVertices.at(0); i++ ) {
		sid = 0 + i;
		//eid = 0 + i + 1;
		eid= (i == _nVertices.at(0) ) ? (0+1) : (0+i+1);
		this->incr_lid();
		Linebase* line = new Linebase(_points[sid], _points[eid], type, this->get_lid() );
		_edges[this->get_lid()] = line;
	}
	//fuse last point into a line ending at the first point, to avoid _nVertices ifs per npoints
	/*sid = 0 + _nVertices.at(0);
	eid = 0 + 1;
	type = INPUT;
	this->incr_lid();
	Linebase* line = new Linebase(_points[sid], _points[eid], type, this->get_lid() );
	_edges[_nVertices.at(0)] = line;*/
	
	//link this only one contour into the datastructure and algorithm of polyTri
	this->initialize();
	myhealth = true;
}


Polygon::Polygon( vector<double>* points )
{
	//set up variables
	//vector and maps deallocate themselves...
	myhealth = false;
	xmin = std::numeric_limits<double>:: max();
	xmax = std::numeric_limits<double>:: lowest();
	ymin = std::numeric_limits<double>:: max();
	ymax = std::numeric_limits<double>:: lowest();
	this->set_lid(0);

	double x,y;
	Type type;
	type = INPUT;

	_nVertices.push_back( points->size() / 2) ; // - 1; //assuming closed contour without holes!
	//unsigned int first = 1;
	//unsigned int last = first + _nVertices.at(0) - 1;

	unsigned int i = 1;
	for (unsigned int j = 0; j < _nVertices.at(0); j++, i++) {
		x = points->at(2*j+0);
		y = points->at(2*j+1);

		if( does_point_exist_already(x, y) ) {
			cout << "Error! I found the duplicate point " << j << "-->" << x << ";" << y << endl;
			return;
			//exit(1);
		} 
		else {
			Pointbase* point = new Pointbase( i, x, y, type );
			if(x > xmax) xmax=x;
			if(x < xmin) xmin=x;
			if(y > ymax) ymax=y;
			if(y < ymin) ymin=y;
			_points[i] = point; //polyTri counts point indices on interval [1,size()]
		}
	}

	unsigned int sid;
	unsigned int eid;
	type = INPUT;
	//run over all points connect to lines starting at sid ending at eid
	for( i = 1; i <= _nVertices.at(0); i++ ) {
		sid = 0 + i;
		//eid = 0 + i + 1;
		eid= (i == _nVertices.at(0) ) ? (0+1) : (0+i+1);
		this->incr_lid();
		Linebase* line = new Linebase(_points[sid], _points[eid], type, this->get_lid() );
		_edges[this->get_lid()] = line;
	}
	//fuse last point into a line ending at the first point, to avoid _nVertices ifs per npoints
	/*sid = 0 + _nVertices.at(0);
	eid = 0 + 1;
	type = INPUT;
	this->incr_lid();
	Linebase* line = new Linebase(_points[sid], _points[eid], type, this->get_lid() );
	_edges[_nVertices.at(0)] = line;*/
	
	//link this only one contour into the datastructure and algorithm of polyTri
	this->initialize();
	myhealth = true;
}

//----------------------------------------------------------------------------
//polygon destructor 
//----------------------------------------------------------------------------
Polygon::~Polygon()
{
	PointbaseMap::iterator itp=_points.begin();
	for(; itp!=_points.end(); itp++)
		delete itp->second;
	
	LineMap::iterator itl=_edges.begin();
	for(; itl!=_edges.end(); itl++)
		delete itl->second;  
}


//---------------------------------------------------------------------------
//check there exists duplicated points or not;
//---------------------------------------------------------------------------
bool Polygon::does_point_exist_already(double x, double y)
{
	PointbaseMap::iterator it=_points.begin();
	for( ; it!=_points.end(); ++it)
		if( it->second->x == x && it->second->y == y )
			return true;
	
	return false;
}

//----------------------------------------------------------------------------
//return the previous point (or edge) id for a given ith point (or edge);
//----------------------------------------------------------------------------
unsigned int Polygon::prev(unsigned int i)
{
	unsigned int j(0),prevLoop(0),currentLoop(0);

	while ( i > _nVertices[currentLoop] ) {
		prevLoop=currentLoop;
		currentLoop++;
	}
 
	if( i==1 || (i==_nVertices[prevLoop]+1) ) j=_nVertices[currentLoop];
	else if( i <= _nVertices[currentLoop] ) j=i-1;

	return j;
}

//----------------------------------------------------------------------------
//return the next point (or edge) id for a given ith point (or edge);
//----------------------------------------------------------------------------
unsigned int Polygon::next(unsigned int i)
{
   unsigned int j(0),prevLoop(0),currentLoop(0);
  
   while ( i > _nVertices[currentLoop] ) 
   {
     prevLoop=currentLoop;	   
     currentLoop++;  
   }
 
   if( i < _nVertices[currentLoop] ) j=i+1;
   else if ( i==_nVertices[currentLoop] ) 
   {
	   if( currentLoop==0) j=1;
	   else j=_nVertices[prevLoop]+1;
   }
   
   return j;
}


//----------------------------------------------------------------------------
//polygon initialization;
//to find types of all polygon vertices;
//create a priority queue for all vertices;
//construct an edge set for each vertex (the set holds all edges starting from 
//the vertex, only for loop searching purpose). 
//----------------------------------------------------------------------------
void Polygon::initialize()
{
	PointbaseMap::iterator it=_points.begin();
	for(; it!=_points.end(); it++) {
		int id=it->first; 
		int idp=prev(id);
		int idn=next(id);
		Pointbase p = *_points[id], pnext = *_points[idn], pprev = *_points[idp];

		if( p > pnext && pprev > p )
			_points[id]->type=REGULAR_DOWN;
		else if (p > pprev && pnext > p)
			_points[id]->type=REGULAR_UP;
		else
		{
			double pa[2], pb[2], pc[2];

			pa[0]=_points[idp]->x;
			pa[1]=_points[idp]->y;

			pb[0]=_points[id]->x;
			pb[1]=_points[id]->y;

			pc[0]=_points[idn]->x;
			pc[1]=_points[idn]->y;

			double area=orient2d(pa,pb,pc);

			if( pprev > p && pnext > p ) _points[id]->type=(area > 0.0) ? END: MERGE ;
			if( pprev < p && pnext < p ) _points[id]->type=(area > 0.0) ? START : SPLIT;
	     }
	    
	    _qpoints.push(*(it->second));

	    _startAdjEdgeMap[id].insert(id);
	}
}

//----------------------------------------------------------------------------
//Add a diagonal from point id i to j
//----------------------------------------------------------------------------  
void Polygon::addDiagonal(unsigned int i, unsigned int j)
{
	Type type=INSERT;
	this->incr_lid();
	Linebase* diag = new Linebase(_points[i], _points[j], type, this->get_lid() );
	_edges[diag->id()]=diag;

	_startAdjEdgeMap[i].insert(diag->id());
	_startAdjEdgeMap[j].insert(diag->id());

	_diagonals[diag->id()]=diag;

	//if(_debug) _logfile<<"Add Diagonal from "<<i<<" to "<<j<<'\n';     
}

//----------------------------------------------------------------------------
//Handle start vertex
//----------------------------------------------------------------------------
void Polygon::handleStartVertex(unsigned int i)
{
     double y=_points[i]->y;	
     _edgebst.InOrder(UpdateKey, y);
     
     _edges[i]->setHelper(i);
     _edges[i]->setKeyValue(y);      
     _edgebst.Insert(_edges[i]);

     /*if(_debug)
     {
     	_logfile<<"set e"<<i<<" helper to "<<i<<'\n';
     	_logfile<<"Insert e"<<i<<" to splay tree\n";
     	_logfile<<"key:"<<_edges[i]->keyValue()<<'\n';
     }*/
}

//----------------------------------------------------------------------------
//Handle end vertex
//----------------------------------------------------------------------------
void Polygon::handleEndVertex(unsigned int i)
{
     double y=_points[i]->y;
     _edgebst.InOrder(UpdateKey, y);
     
     unsigned int previ=prev(i);
     Linebase* edge=_edges[previ];
     unsigned int helper=_edges[previ]->helper();
    
     
     if(_points[helper]->type==MERGE) addDiagonal(i, helper);
     _edgebst.Delete(edge->keyValue());

    /*if(_debug)
    {
     	_logfile<<"Remove e"<<previ<<" from splay tree\n";
     	_logfile<<"key:"<<edge->keyValue()<<'\n';
    }*/
}

//----------------------------------------------------------------------------
//Handle split vertex
//----------------------------------------------------------------------------
void Polygon::handleSplitVertex(unsigned int i)
{
     double x=_points[i]->x, y=_points[i]->y;
     _edgebst.InOrder(UpdateKey, y);
     
     BTreeNode<Linebase*, double>*  leftnode;
     _edgebst.FindMaxSmallerThan(x, leftnode);
     Linebase* leftedge=leftnode->data();
     
     unsigned int helper=leftedge->helper();
     addDiagonal(i, helper);

     /*if(_debug)
     {
     	_logfile<<"Search key:"<<x<<" edge key:"<<leftedge->keyValue()<<'\n';
     	_logfile<<"e"<<leftedge->id()<<" is directly left to v"<<i<<'\n';  
     	_logfile<<"Set e"<<leftedge->id()<<" helper to "<<i<<'\n';
     	_logfile<<"set e"<<i<<" helper to "<<i<<'\n';
     	_logfile<<"Insert e"<<i<<" to splay tree\n"; 
     	_logfile<<"Insert key:"<<_edges[i]->keyValue()<<'\n';
     }*/
    
     leftedge->setHelper(i);
     _edges[i]->setHelper(i);
     _edges[i]->setKeyValue(y);
     _edgebst.Insert(_edges[i]);  
}


//----------------------------------------------------------------------------
//Handle merge vertex
//----------------------------------------------------------------------------
void Polygon::handleMergeVertex(unsigned int i)
{
     double x=_points[i]->x, y=_points[i]->y;
     _edgebst.InOrder(UpdateKey, y);
    
     unsigned int previ=prev(i);
     unsigned int helper=_edges[previ]->helper();
     if (_points[helper]->type==MERGE) addDiagonal(i, helper);
     _edgebst.Delete(_edges[previ]->keyValue());
     /*if(_debug)
     {
	     _logfile<<"e"<<previ<<" helper is "<<helper<<'\n';
	     _logfile<<"Remove e"<<previ<<" from splay tree.\n";
     }*/

     BTreeNode<Linebase*, double>*  leftnode;
     _edgebst.FindMaxSmallerThan(x, leftnode);
     Linebase* leftedge=leftnode->data();
            
     helper=leftedge->helper();
     if(_points[helper]->type==MERGE) addDiagonal(i, helper);
     
     leftedge->setHelper(i);

     /*if(_debug)
     {
     	_logfile<<"Search key:"<<x<<" found:"<<leftedge->keyValue()<<'\n';
     	_logfile<<"e"<<leftedge->id()<<" is directly left to v"<<i<<'\n';
     	_logfile<<"Set e"<<leftedge->id()<<" helper to "<<i<<'\n';
     }*/
}

//----------------------------------------------------------------------------
//Handle regular down vertex
//----------------------------------------------------------------------------
void Polygon::handleRegularVertexDown(unsigned int i)
{
     double y=_points[i]->y;
     _edgebst.InOrder(UpdateKey, y);
     
     unsigned int previ=prev(i);
     unsigned int helper=_edges[previ]->helper();
     if(_points[helper]->type==MERGE) addDiagonal(i, helper); 
	  
     _edgebst.Delete(_edges[previ]->keyValue());
     _edges[i]->setHelper(i);
     _edges[i]->setKeyValue(y);
     _edgebst.Insert(_edges[i]);

     /*if(_debug)
     {
	   _logfile<<"e"<<previ<<" helper is "<<helper<<'\n';
	   _logfile<<"Remove e"<<previ<<" from splay tree.\n";
	   _logfile<<"Set e"<<i<<" helper to "<<i<<'\n'; 
	   _logfile<<"Insert e"<<i<<" to splay tree\n"; 
	   _logfile<<"Insert key:"<<_edges[i]->keyValue()<<'\n';
     }*/
}


//----------------------------------------------------------------------------
////Handle regular up vertex
//----------------------------------------------------------------------------
void Polygon::handleRegularVertexUp(unsigned int i)
{	
     double x=_points[i]->x, y=_points[i]->y;
     _edgebst.InOrder(UpdateKey, y);
     
     BTreeNode<Linebase*, double>*  leftnode;
     _edgebst.FindMaxSmallerThan(x, leftnode);
     
     Linebase* leftedge=leftnode->data();
		     
     unsigned int helper=leftedge->helper();
     if(_points[helper]->type==MERGE) addDiagonal(i, helper);
     leftedge->setHelper(i);

     /*if(_debug)
     {
     	_logfile<<"Search key:"<<x<<" found:"<<leftedge->keyValue()<<'\n';
     	_logfile<<"e"<<leftedge->id()<<" is directly left to v"<<i<<" and its helper is:"<<helper<<'\n';
     	_logfile<<"Set e"<<leftedge->id()<<" helper to "<<i<<'\n';
     }*/
}

//----------------------------------------------------------------------------
//partition polygon to monotone polygon pieces
//----------------------------------------------------------------------------
bool Polygon::partition2Monotone()
{
	if(_qpoints.top().type!=START) {
		cout<<"Partition::Please check your input polygon:\n1)orientations?\n2)duplicated points?\n" <<" poly2tri stopped.\n";
		return false;
		//exit(1);
	}

	while(!_qpoints.empty()) {
		Pointbase vertex=_qpoints.top();
		_qpoints.pop();
		unsigned int id=vertex.id;

	  /*if(_debug)
	  {
		 string stype;
	 	 switch (vertex.type)
	  	{
	     		case START:        stype="START";       break;
			case END:          stype="END";         break;
			case MERGE:        stype="MERGE";       break;
			case SPLIT:        stype="SPLIT";       break;
			case REGULAR_UP:   stype="REGULAR_UP";  break;
			case REGULAR_DOWN: stype="REGULAR_DOWN";break;
			default: 
					   cout<<"No duplicated points please! poly2tri stopped\n";
					   exit(1); break;
		}
				    
	  	_logfile<<"\n\nHandle vertex:"<<vertex.id<<" type:"<<stype<<'\n';
	  }*/


		switch(vertex.type) {
			case START:        handleStartVertex(id);       break;
			case END:          handleEndVertex(id);         break;
			case MERGE:        handleMergeVertex(id);       break;
			case SPLIT:        handleSplitVertex(id);       break;
			case REGULAR_UP:   handleRegularVertexUp(id);   break;
			case REGULAR_DOWN: handleRegularVertexDown(id); break;
			default:
				cout<<"Partition::No duplicated points please! poly2tri stopped\n";
				return false; //exit(1); break;                         
		}
	}
	return true;
}


//----------------------------------------------------------------------------
//two Auxiliary functions to find monotone polygon pieces
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//calculate angle B for A, B, C three given points
//----------------------------------------------------------------------------
double Polygon::angleCosb(double *pa, double *pb, double *pc)
{
  double dxab = pa[0] - pb[0];
  double dyab = pa[1] - pb[1];

  double dxcb = pc[0] - pb[0];
  double dycb = pc[1] - pb[1];

  double dxab2 = dxab * dxab;
  double dyab2 = dyab * dyab;
  double dxcb2 = dxcb * dxcb;
  double dycb2 = dycb * dycb;
  double ab = dxab2 + dyab2;
  double cb = dxcb2 + dycb2;

  double cosb = dxab * dxcb + dyab * dycb;
  double denom = sqrt( ab * cb);
  
  cosb/=denom;
  
  return cosb;
}

//----------------------------------------------------------------------------
//for any given edge, find the next edge we should choose when searching for
//monotone polygon pieces; 
//----------------------------------------------------------------------------
unsigned int Polygon::selectNextEdge(Linebase* edge)
{

    unsigned int eid= edge->endPoint(1)->id;
    set<unsigned int> edges=_startAdjEdgeMap[eid];
    assert(!edges.empty());
     
    unsigned int nexte=0;
    if( edges.size() == 1 )  nexte=*(edges.begin());
    else if( edges.size() > 1 )
    {
	unsigned int nexte_ccw(0), nexte_cw(0);
	double max=-2.0,min=2.0;
	
    	
	set<unsigned int>::iterator it=edges.begin();
	for(; it!=edges.end(); it++)
	{
		if(*it==edge->id()) continue;     
		double A[2], B[2], C[2];
		A[0]=edge->endPoint(0)->x;        A[1]=edge->endPoint(0)->y;
		B[0]=edge->endPoint(1)->x;        B[1]=edge->endPoint(1)->y;
		
		if(edge->endPoint(1)!=_edges[*it]->endPoint(0)) _edges[*it]->reverse();
		C[0]=_edges[*it]->endPoint(1)->x; C[1]=_edges[*it]->endPoint(1)->y;
		
		double area=orient2d(A, B, C);
		double cosb=angleCosb(A, B, C);

		if( area > 0 && max < cosb ) { nexte_ccw=*it; max=cosb; }
		else if( min > cosb ) { nexte_cw=*it; min=cosb; }
	}

	nexte = (nexte_ccw!=0) ? nexte_ccw : nexte_cw;
    } 

   return nexte; 
}

//----------------------------------------------------------------------------
//searching all monotone pieces;
//----------------------------------------------------------------------------
bool Polygon::searchMonotones()
{
	int loop=0;
	LineMap edges=_edges;
	PointbaseMap points=_points; //##MK::to assure that the per se infinite for loop inside while at some point stops for sure

	while( edges.size() > _diagonals.size() ) {
		loop++;
		Monopoly poly;
		LineMap::iterator it=edges.begin();
		Pointbase* startp=it->second->endPoint(0); //##MK::startp=startp=it->second->endPoint(0);
		Pointbase* endp=0;
		Linebase*  next=it->second;

		poly.push_back(startp->id);

		/*if(_debug)
		{
		_logfile<<"Searching for loops:"<<loop<<'\n';
		_logfile<<"vertex index:"<<startp->id<<" ";
		}*/

		long trips = 0;
		long maxtrips = 100 * _points.size(); //heuristically no polygon requires partitioning into more than two orders of magnitude more triangles than points along the 2d contour
		bool success = false;
		for( ; trips < maxtrips ; ) {
			trips++;
			endp=next->endPoint(1);
			if(next->type()!=INSERT) {
				edges.erase(next->id());
				_startAdjEdgeMap[next->endPoint(0)->id].erase(next->id());
			}
			if(endp==startp) { //only entry for successfully partitioned polygon
				success = true;
				break; 
			}
			poly.push_back(endp->id);

			//if(_debug) _logfile<<endp->id<<" ";

			unsigned int nexte=selectNextEdge(next);
			if(nexte==0) {
				cout<<"SearchMonotones::Please check your input polygon:\n";
				cout<<"1)orientations?\n2)with duplicated points?\n3)is a simple one?\n";
				cout<<"poly2tri stopped.\n";
				return false; //exit(1);
			}
			//assert( nexte > 0);
			next=edges[nexte];
			if(next->endPoint(0) !=endp ) 
				next->reverse(); 
		}

		if ( success == false ) {
			cerr << "WARNING::PolyTri suffered problems in searchMonotone!\n";
			return false;
		}

		//if(_debug) _logfile<<"\nloop closed!\n\n";

		_mpolys.push_back(poly);
	}

	return true;
}


//----------------------------------------------------------------------------
//triangulate a monotone polygon;
//----------------------------------------------------------------------------
bool Polygon::triangulateMonotone(Monopoly& mpoly)
{
	PQueue qvertex;
	Monopoly::iterator it=mpoly.begin(), itnext;
	for(; itnext=it, it!=mpoly.end(); it++) {
		itnext++;  
		if(itnext==mpoly.end()) itnext=mpoly.begin();
		Pointbase point=*_points[*it], pointnext=*_points[*itnext];
		point.left=(point > pointnext)? true:false;
		qvertex.push(point);
	}
 
	stack<Pointbase> spoint;
	for(int i=0; i<2; i++) { 
		spoint.push(qvertex.top()); qvertex.pop();
	}

	while ( qvertex.size() > 1 ) {
		Pointbase topQueuePoint=qvertex.top();
		Pointbase topStackPoint=spoint.top();
		if(topQueuePoint.left!=topStackPoint.left) {
			while ( spoint.size()  > 1 ) {
				Pointbase p1=spoint.top();
				spoint.pop(); 
				Pointbase p2=spoint.top();
				Triangle v(3);
				v[0]=topQueuePoint.id;
				v[1]=p1.id;
				v[2]=p2.id;
				_triangles.push_back(v);
#ifdef POLYTRI_DEBUG
				cout << "Add triangle " << v[0] << " " << v[1] << " " << v[2]<< endl;
#endif
			}
			spoint.pop();
			spoint.push(topStackPoint);
			spoint.push(topQueuePoint);
		}
		else {
			while( spoint.size() > 1 ) {
				Pointbase stack1Point=spoint.top();
				spoint.pop(); 
				Pointbase stack2Point=spoint.top();
				spoint.push(stack1Point);
				double pa[2], pb[2], pc[2];
				pa[0]=topQueuePoint.x; pa[1]=topQueuePoint.y;
				pb[0]=stack2Point.x;   pb[1]=stack2Point.y;
				pc[0]=stack1Point.x;   pc[1]=stack1Point.y;
			
				/*if(_debug)
				{
						_logfile<<"current top queue vertex index="<<topQueuePoint.id<<'\n';
						_logfile<<"Current top stack vertex index="<<stack1Point.id<<'\n';
						_logfile<<"Second stack vertex index="<<stack2Point.id<<'\n';
				}*/

				double area=orient2d(pa,pb,pc);
				bool left=stack1Point.left;
				if( (area > 0 && left) || (area < 0 && !left ) ) {
					Triangle v(3);
					v[0]=topQueuePoint.id;
					v[1]=stack2Point.id;
					v[2]=stack1Point.id;
					_triangles.push_back(v);
#ifdef POLYTRI_DEBUG
					cout << "Add triangle " << v[0] << " " << v[1] << " " << v[2]<< endl;
#endif
					spoint.pop();
				} else break;
			}
			spoint.push(topQueuePoint);
		}
		qvertex.pop();
	}

	Pointbase lastQueuePoint=qvertex.top();
	while( spoint.size() != 1 ) {
		Pointbase topPoint=spoint.top();
		spoint.pop();
		Pointbase top2Point=spoint.top();

		Triangle v(3);
		v[0]=lastQueuePoint.id;
		v[1]=topPoint.id;
		v[2]=top2Point.id;
		_triangles.push_back(v);
#ifdef POLYTRI_DEBUG
		cout << "Add triangle " << v[0] << " " << v[1] << " " << v[2]<< endl;
#endif
	}

	return true;
}


bool Polygon::triangulation()
{
	if ( partition2Monotone() == false ) {
		verboseTriangles();
		return false;
	}
	if ( searchMonotones() == false ) {
		verboseTriangles();
		return false;
	}
	for( Monopolys::iterator it=_mpolys.begin(); it!=_mpolys.end(); it++ ) {
		if ( triangulateMonotone(*it) == false ) {
			verboseTriangles();
			return false;
		}
	}

#ifdef POLYTRI_DEBUG
	cout << "Total number of triangles:" << _triangles.size() << '\n';
#endif
	return true;
}


void Polygon::result()
{
	for ( Triangles::iterator it = _triangles.begin(); it != _triangles.end(); it++ ) {
		cout << _points[(*it).at(0)]->x << ";" << _points[(*it).at(0)]->y << ";" << "\t";
		cout << _points[(*it).at(1)]->x << ";" << _points[(*it).at(1)]->y << "\t";
		cout << _points[(*it).at(2)]->x << ";" << _points[(*it).at(2)]->y << endl;
	}
}


vector<triangle>* Polygon::giveTriangles()
{
	vector<triangle>* these = NULL;
	these = new vector<triangle>; 
	//##MK::attempt a resize you get full page any way?

	if ( _triangles.size() > 0 ) {
		for ( Triangles::iterator it = _triangles.begin(); it != _triangles.end(); it++ ) {
			struct triangle t;
			//cout << _points[(*it).at(0)]->x << ";" << _points[(*it).at(0)]->y << ";" << "\t";
			//cout << _points[(*it).at(1)]->x << ";" << _points[(*it).at(1)]->y << "\t";
			//cout << _points[(*it).at(2)]->x << ";" << _points[(*it).at(2)]->y << endl;
			t.x1 = _points[(*it).at(0)]->x;		t.y1 =  _points[(*it).at(0)]->y;
			t.x2 = _points[(*it).at(1)]->x;		t.y2 =  _points[(*it).at(1)]->y;
			t.x3 = _points[(*it).at(2)]->x;		t.y3 =  _points[(*it).at(2)]->y;
			these->push_back( t );
		}
	}
	return these;
}

void Polygon::verboseTriangles()
{
	//vector<triangle>* these = NULL;
	if ( _triangles.size() > 0 ) {
		for ( Triangles::iterator it = _triangles.begin(); it != _triangles.end(); it++ ) {
			cout << setprecision(12) << _points[(*it).at(0)]->x << ";" << setprecision(12) << _points[(*it).at(0)]->y << ";" << "\t\t";
			cout << setprecision(12) << _points[(*it).at(1)]->x << ";" << setprecision(12) << _points[(*it).at(1)]->y << "\t\t";
			cout << setprecision(12) << _points[(*it).at(2)]->x << ";" << setprecision(12) << _points[(*it).at(2)]->y << endl;
		}
	}
}
