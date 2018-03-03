// Contains routines to initialize and relax wake

#include <iostream>
#include "body.h"
#include "wake.h"
#include "settings.h"
#include "edge.h"
#include "node.h"

/******************************************************************************/
//
// Initializes the wake
//
/******************************************************************************/
void initialize_wake ( Body & aircraft, Wake & wake )
{
  unsigned int i, j, nedges, nwakelines;
  bool foundnode;
  Edge * edge;
  Node * node1, * node2;

  std::cout << "Initializing wake ... " << std::endl;

  // Loop through trailing edges to create wake lines

  nedges = aircraft.numEdges();
  for ( i = 0; i < nedges; i++ )
  {
    // Determine if it is a trailing edge

    edge = & aircraft.edge(i);
    if (edge->isTrailingEdge())
    {
      // Check if node 1 already has an associated wake line

      node1 = & edge->node(0);
      nwakelines = wake.numWakeLines();
      foundnode = false;
      for ( j = 0; j < nwakelines; j++ )
      {
        if (wake.wakeLine(j).TENode().label() == node1->label())
        {
          foundnode = true;
          break;
        }
      } 
      
      // Add wake line if needed

      if (not foundnode) 
      { 
        wake.addWakeLine(deforming_wake_length, num_filaments, uinf,
                         & aircraft.edge(i).node(0));
      }

      // Check if node 2 already has an associated wake line

      node2 = & edge->node(1);
      nwakelines = wake.numWakeLines();
      foundnode = false;
      for ( j = 0; j < nwakelines; j++ )
      {
        if (wake.wakeLine(j).TENode().label() == node2->label())
        {
          foundnode = true;
          break;
        }
      } 
      
      // Add wake line if needed

      if (not foundnode) 
      { 
        wake.addWakeLine(deforming_wake_length, num_filaments, uinf,
                         & aircraft.edge(i).node(1));
      }
    } /** (End check for trailing edge) **/
  } /** (End loop over edges) **/

  // Determine core radius

  // Set up horseshoe vortices

  // Print wake information

  std::cout << "Successfully initialized wake." << std::endl;
  std::cout << "Wake information: " << std::endl;
  std::cout << "Wake lines: " << wake.numWakeLines() << std::endl;
  std::cout << std::endl;
}
