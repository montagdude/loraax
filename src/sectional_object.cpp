#include "sectional_object.h"

/******************************************************************************/
//
// SectionalObject class. Super class for objects stored at discrete spanwise
// locations.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
SectionalObject::SectionalObject ()
{
  _y = 0.;
}

/******************************************************************************/
//
// Set/access data
//
/******************************************************************************/
void SectionalObject::setY ( const double & y ) { _y = y; }
const double & SectionalObject::y () const { return _y; }
