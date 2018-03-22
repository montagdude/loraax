// Header for sectional object class

#ifndef SECTIONALOBJ_H
#define SECTIONALOBJ_H

/******************************************************************************/
//
// SectionalObject class. Super class for objects stored at discrete spanwise
// locations.
//
/******************************************************************************/
class SectionalObject {

  protected:

    double _y;

  public:

    // Constructor

    SectionalObject (); 

    // Set or access data
    
    virtual void setY ( const double & y );
    virtual const double & y () const;
};

#endif
