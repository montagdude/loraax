#include "util.h"
#include "node.h"
#include "edge.h"
#include "face.h"
#include "triface.h"
#include "quadface.h"
#include "body.h"
#include <iostream>

/******************************************************************************/
//
// Body class. Stores groups of nodes, edges, and faces making up solid bodies.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Body::Body ()
{
  _nnodes = 0;
  _nedges = 0;
  _nfaces = 0;
  _nthin = 0;
  _nte = 0;
  _ntefaces = 0;
}

/******************************************************************************/
//
// Returns node array index from label
//
/******************************************************************************/
unsigned int Body::nodeIndexFromLabel ( unsigned int lbl ) const
{
  unsigned int counter, idx;
  bool found;

  counter = 0;
  found = false;
  while (not found)
  {

    // Check for requested node
    
    if (_nodes[counter].label() == lbl)
    {
      found = true;
      idx = counter;
    }
    counter += 1;
 
    // Error if node is not found in the array

    if (counter == _nnodes+1)
    {
      conditional_stop(1, "Body::nodeIndexFromLabel", 
                       "Invalid node label specified.");
    }
  }

  return idx;
}

/******************************************************************************/
//
// Adds a node to the _nodes array given label and coordinates
//
/******************************************************************************/
void Body::addNode ( unsigned int lbl, const double & x, const double & y,
                     const double & z )
{
  Node node(lbl, x, y, z);
  _nodes.push_back(node);
  _nnodes += 1;
}

/******************************************************************************/
//
// Returns a node by array index
//
/******************************************************************************/
Node & Body::node ( unsigned int nidx ) { return _nodes[nidx]; }

/******************************************************************************/
//
// Adds a triangular face to the _tris and _faces arrays given label and
// list of node labels
//
/******************************************************************************/
void Body::addFace ( unsigned int lbl, unsigned int nlbl0, unsigned int nlbl1,
                     unsigned int nlbl2, bool thin )
{
  unsigned int nidx0, nidx1, nidx2;
  TriFace tface;

  // Get node array indices

  nidx0 = nodeIndexFromLabel(nlbl0);
  nidx1 = nodeIndexFromLabel(nlbl1);
  nidx2 = nodeIndexFromLabel(nlbl2);

  // Set up face 

  tface.setNode(0, & _nodes[nidx0]);
  tface.setNode(1, & _nodes[nidx1]);
  tface.setNode(2, & _nodes[nidx2]);
  tface.setLabel(lbl);
  if (thin) 
  { 
    tface.setThin(); 
    _nthin += 1;
  }

  // Add to _tris array and increment face counter

  _tris.push_back(tface);
  _nfaces += 1;
}

/******************************************************************************/
//
// Adds a triangular face to the _quads and _faces arrays given label and
// list of node labels
//
/******************************************************************************/
void Body::addFace ( unsigned int lbl, unsigned int nlbl0, unsigned int nlbl1,
                     unsigned int nlbl2, unsigned int nlbl3, bool thin )
{
  unsigned int nidx0, nidx1, nidx2, nidx3;
  QuadFace qface;

  // Get node array indices

  nidx0 = nodeIndexFromLabel(nlbl0);
  nidx1 = nodeIndexFromLabel(nlbl1);
  nidx2 = nodeIndexFromLabel(nlbl2);
  nidx3 = nodeIndexFromLabel(nlbl3);

  // Set up face 

  qface.setNode(0, & _nodes[nidx0]);
  qface.setNode(1, & _nodes[nidx1]);
  qface.setNode(2, & _nodes[nidx2]);
  qface.setNode(3, & _nodes[nidx3]);
  qface.setLabel(lbl);
  if (thin) 
  { 
    qface.setThin(); 
    _nthin += 1;
  }

  // Add to _quads array and increment face counter

  _quads.push_back(qface);
  _nfaces += 1;
}

/******************************************************************************/
//
// Sets pointers for generic (base class) faces to trifaces and quadfaces. Also
// sets face references for nodes.
//
/******************************************************************************/
void Body::setFacePointers ()
{
  unsigned int i, j, numtris, numquads;

  std::cout << "Setting face pointers ..." << std::endl;

  numtris = _tris.size();
  numquads = _quads.size();

  // Reserve space for _faces array

  _faces.resize(numtris + numquads);

  // Add tris first, then quads

  for ( i = 0; i < numtris; i++ ) 
  { 
    _faces[i] = & _tris[i]; 

    // Allocate space for face references
    
    for ( j = 0; j < 3; j++ ) { _tris[i].node(j).increaseNumFaces(); }
  }
  for ( i = 0; i < numquads; i++ ) 
  { 
    _faces[i+numtris] = & _quads[i]; 

    // Allocate space for face references
    
    for ( j = 0; j < 4; j++ ) { _quads[i].node(j).increaseNumFaces(); }
  }

  // Set tri and quad first/last indices

  _firsttri = 0;
  _lasttri = numtris - 1;
  _firstquad = numtris;
  _lastquad = numtris + numquads - 1;

  // Add face references for nodes (tris first, then quads)

  for ( i = 0; i < numtris; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      _tris[i].node(j).setNextFace(& _tris[i]);
    }
  }
  for ( i = 0; i < numquads; i++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      _quads[i].node(j).setNextFace(& _quads[i]);
    }
  }
}

/******************************************************************************/
//
// Returns a face by array index
//
/******************************************************************************/
Face & Body::face ( unsigned int fidx ) { return *_faces[fidx]; }

/******************************************************************************/
//
// Creates edges from face and node information
//
/******************************************************************************/
void Body::constructEdges ()
{
  unsigned int i, j, numtris, numquads;
  int edgeid;
  Node * nodelist[5];

  std::cout << "Constructing edges ..." << std::endl;

  numtris = _tris.size();
  numquads = _quads.size();

  // Loop through tri faces, constructing edges

  for ( i = 0; i < numtris; i++ )
  {
    // Get node references

    nodelist[0] = & _tris[i].node(0); 
    nodelist[1] = & _tris[i].node(1); 
    nodelist[2] = & _tris[i].node(2); 
    nodelist[3] = & _tris[i].node(0); 

    // Create edge first if it doesn't already exist

    for ( j = 0; j < 3; j++ )
    {
      edgeid = edgeFromNodes(nodelist[j]->label(), nodelist[j+1]->label());
      if (edgeid == -1)
      {
        Edge edge;
        edge.setIndex(_nedges);
        edge.setNode(0, nodelist[j]);
        edge.setNode(1, nodelist[j+1]);
        edge.increaseNumFaces(); 
        _edges.push_back(edge);
        _nedges += 1;
      }
  
      // If it does exist, just increase number of faces
  
      else { _edges[edgeid].increaseNumFaces(); }
    }
  }

  // Loop through quad faces, constructing edges

  for ( i = 0; i < numquads; i++ )
  {
    // Get node references

    nodelist[0] = & _quads[i].node(0); 
    nodelist[1] = & _quads[i].node(1); 
    nodelist[2] = & _quads[i].node(2); 
    nodelist[3] = & _quads[i].node(3); 
    nodelist[4] = & _quads[i].node(0); 

    // Create first edge if it doesn't already exist

    for ( j = 0; j < 4; j++ )
    {
      edgeid = edgeFromNodes(nodelist[j]->label(), nodelist[j+1]->label());
      if (edgeid == -1)
      {
        Edge edge;
        edge.setIndex(_nedges);
        edge.setNode(0, nodelist[j]);
        edge.setNode(1, nodelist[j+1]); 
        edge.increaseNumFaces();
        _edges.push_back(edge);
        _nedges += 1;
      }
  
      // If it does exist, just increase number of faces
  
      else { _edges[edgeid].increaseNumFaces(); }
    }
  }

  // Set face-edge connections on tris
  
  for ( i = 0; i < numtris; i++ )
  {
    // Get node references

    nodelist[0] = & _tris[i].node(0); 
    nodelist[1] = & _tris[i].node(1); 
    nodelist[2] = & _tris[i].node(2); 
    nodelist[3] = & _tris[i].node(0); 

    // Connect edges to face and vice-versa

    for ( j = 0; j < 3; j++ )
    {
      edgeid = edgeFromNodes(nodelist[j]->label(), nodelist[j+1]->label());
      if (edgeid == -1)
      {
        conditional_stop(1, "Body::constructEdges",
                         "Edge should have already been constructed!");
      }
      _tris[i].setNextEdge(& _edges[edgeid]);
      _edges[edgeid].setNextFace(& _tris[i]);
    }
  }

  // Set face-edge connections on quads
  
  for ( i = 0; i < numquads; i++ )
  {
    // Get node references

    nodelist[0] = & _quads[i].node(0); 
    nodelist[1] = & _quads[i].node(1); 
    nodelist[2] = & _quads[i].node(2); 
    nodelist[3] = & _quads[i].node(3); 
    nodelist[4] = & _quads[i].node(0); 

    // Connect edges to face and vice-versa

    for ( j = 0; j < 4; j++ )
    {
      edgeid = edgeFromNodes(nodelist[j]->label(), nodelist[j+1]->label());
      if (edgeid == -1)
      {
        conditional_stop(1, "Body::constructEdges",
                         "Edge should have already been constructed!");
      }
      _quads[i].setNextEdge(& _edges[edgeid]);
      _edges[edgeid].setNextFace(& _quads[i]);
    }
  }
}

/******************************************************************************/
//
// Returns an edge by array index
//
/******************************************************************************/
Edge & Body::edge ( unsigned int eidx ) { return _edges[eidx]; }

/******************************************************************************/
//
// Returns edge index given two node labels
//
/******************************************************************************/
int Body::edgeFromNodes ( unsigned int nlbl0, unsigned int nlbl1 ) const
{
  unsigned int i;
  int edgeid;

  edgeid = -1;   // Default if not found

  // Loop through edges, looking for the edge with the correct endpoints

  for ( i = 0; i < _nedges; i++ )
  {
    if (  ((_edges[i].node(0).label() == nlbl0) and
           (_edges[i].node(1).label() == nlbl1))
       or ((_edges[i].node(1).label() == nlbl0) and
           (_edges[i].node(0).label() == nlbl1)) )
    {
      edgeid = int(_edges[i].index());
      break;
    }
  }
    
  return edgeid;
} 

/******************************************************************************/
//
// Sets an edge as a trailing edge from node labels.  Also sets connected faces
// as upper or lower trailing edge faces.
//
/******************************************************************************/
void Body::setTrailingEdgeFromNodes ( unsigned int nlbl0, unsigned int nlbl1 )
{
  unsigned int edgeid, nfaces;
  double nz0, nz1;
  Face * face0, * face1;

  // Determine edge index in array

  edgeid = edgeFromNodes(nlbl0, nlbl1);

  // Set edge as trailing edge

  _edges[edgeid].setTrailingEdge();

  // Increment number of trailing edges

  _nte += 1;

  // Set connected faces as upper TE or lower TE

  nfaces = _edges[edgeid].numFaces();
  if (nfaces == 1)                        // Thin surfaces
  {
    _ntefaces += 1;
    face0 = & _edges[edgeid].face(0);
    face0->setUpperTrailingEdge();
  }
  else if (nfaces == 2)                   // Thick surfaces
  {
    _ntefaces += 2;
    face0 = & _edges[edgeid].face(0);
    face1 = & _edges[edgeid].face(1);
  
    // Determine top and bottom by z-components of normals

    nz0 = face0->normal()(2);
    nz1 = face1->normal()(2);

    if (nz0 > nz1)
    {
      face0->setUpperTrailingEdge();
      face1->setLowerTrailingEdge();
    }
    else if (nz0 < nz1)
    {
      face0->setLowerTrailingEdge();
      face1->setUpperTrailingEdge();
    }
    else
    {
      conditional_stop(1, "Body::setTrailingEdgeFromNodes",
                     "Trailing edge faces have identical z-normal components.");
    }
  }
  else                                    // Error (>2 faces connected to TE)
  {
    conditional_stop(1, "Body::setTrailingEdgeFromNodes",
                    "Trailing edges cannot be connected to more than 2 faces.");
  }
}

/******************************************************************************/
//
// Returns number of nodes
//
/******************************************************************************/
unsigned int Body::numNodes () const { return _nnodes; }

/******************************************************************************/
//
// Returns number of faces
//
/******************************************************************************/
unsigned int Body::numFaces () const { return _nfaces; }

/******************************************************************************/
//
// Returns number of tris
//
/******************************************************************************/
unsigned int Body::numTris () const { return _lasttri + 1; }

/******************************************************************************/
//
// Returns number of quads
//
/******************************************************************************/
unsigned int Body::numQuads () const { return _lastquad - _lasttri; }

/******************************************************************************/
//
// Returns number of faces belonging to thin (zero thickness) surfaces
//
/******************************************************************************/
unsigned int Body::numThin () const { return _nthin; }

/******************************************************************************/
//
// Returns number of trailing edge faces
//
/******************************************************************************/
unsigned int Body::numTEFaces () const { return _ntefaces; }

/******************************************************************************/
//
// Returns number of edges
//
/******************************************************************************/
unsigned int Body::numEdges () const { return _nedges; }

/******************************************************************************/
//
// Returns number of trailing edges
//
/******************************************************************************/
unsigned int Body::numTrailingEdges () const { return _nte; }

/******************************************************************************/
//
// Array index of first tri
//
/******************************************************************************/
unsigned int Body::firstTri () const { return _firsttri; }

/******************************************************************************/
//
// Array index of last tri
//
/******************************************************************************/
unsigned int Body::lastTri () const { return _lasttri; }

/******************************************************************************/
//
// Array index of first quad
//
/******************************************************************************/
unsigned int Body::firstQuad() const { return _firstquad; }

/******************************************************************************/
//
// Array index of last quad
//
/******************************************************************************/
unsigned int Body::lastQuad () const { return _lastquad; }
