// Header for face class

#ifndef FACE_H
#define FACE_H

#include <vector>
#include <string>
#include <Eigen/Core>

// Forward declarations

class Node;
class Edge;

/******************************************************************************/
//
// Face class.  Stores references to connected edges and nodes, computes
// source/doublet influence coefficients at a point, etc.
//
/******************************************************************************/
class Face {

  protected:

    unsigned int _lbl;             // Face label
    bool _thinflag;                // Flags for thin surface (0 thickness)
    bool _TEuflag, _TElflag;       // Flags for upper/lower lifting TE faces
    double _sigma, _mu;            // Source and doublet strength
    double _length;                // A characteristic length
    double _area;                  // Face area
    Eigen::Vector3d _norm;         // Outward-facing normal
    Eigen::Vector3d _cen;          // Face centroid
    Eigen::Matrix3d _trans, _invtrans;
                                   // Transform from inertial frame to panel
                                   // frame and vice versa
    std::vector<double> _xtrans, _ytrans;
                                   // Panel node coordinates in panel frame

    std::vector<Node*> _noderef;   // Vector of connected node pointers
    std::vector<Edge*> _edgeref;   // Vector of connected edge pointers

    unsigned int _currnodes, _curredges;
                                   // # of nodes and edges that have been set

  public:

    // Constructor

    Face ();

    // Setting and accessing face label

    void setLabel ( unsigned int );
    unsigned int label () const;

    // Access to nodes

    unsigned int numNodes () const;
    virtual void setNode ( unsigned int, Node * );
    Node & node ( unsigned int ) const;

    // Access to edges

    unsigned int numEdges () const;
    void setNextEdge ( Edge * );
    Edge & edge ( unsigned int ) const;

    // Setting and accessing source and doublet strength

    void setSourceStrength ( const double & );
    const double & sourceStrength () const;

    void setDoubletStrength ( const double & );
    const double & doubletStrength () const;

    // Returns geometric quantities

    const double & area () const;
    const Eigen::Vector3d & centroid () const;
    const Eigen::Vector3d & normal () const;

    // Generic querying functions

    const double & getScalar ( const std::string & ) const;
    const Eigen::Vector3d & getVector ( const std::string & ) const;

    // Setting and querying thin face flag (zero thickness surfaces)

    void setThin ();
    bool isThin () const;

    // Setting and querying trailing edge face properties

    void setUpperTrailingEdge ();
    void setLowerTrailingEdge ();
    bool isTrailingEdge () const;
    std::string trailingEdgeType () const;

    // Virtual functions implemented in derived classes
    
    virtual double sourcePhiCoeff ( const double &, const double &,
                                    const double &, const bool,
                                    const std::string & ) const = 0;
    virtual Eigen::Vector3d sourceVCoeff ( const double &, const double &,
                                           const double &, const bool,
                                           const std::string & ) const = 0;
    virtual double doubletPhiCoeff ( const double &, const double &,
                                     const double &, const bool,
                                     const std::string & ) const = 0;
    virtual Eigen::Vector3d doubletVCoeff ( const double &, const double &,
                                            const double &, const bool,
                                            const std::string & ) const = 0;
    virtual double inducedPotential ( const double &, const double &,
                                      const double &, const bool,
                                      const std::string & ) const = 0;
    virtual Eigen::Vector3d inducedVelocity ( const double &, const double &,
                                              const double &, const bool,
                                              const std::string & ) const = 0;

};

#endif
