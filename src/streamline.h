/************************************************************************/
/*                                                                      */
/*  written by  U.Lehnert                               08/2012         */
/*                                                                      */
/************************************************************************/

#ifndef STREAMLINE_H
#define STREAMLINE_H

#include <math.h>
#include <stdio.h>	// for debugging only

#include "vector.h"

using namespace std;

// a streamline starting at a given point and extending to xInfinity
// is is represented as a sequence of straight line segments

class Streamline
{

  public:

    // create a streamline starting at A with a predefined number of segments
    Streamline(Vector A, double xInfinity, int nSegment);
    // create a streamline starting at A with a predefined direction
    // and length of the first segment : P(1) = P(0) + startDirection
    Streamline(Vector A, double xInfinity, int nSegment, Vector startDirection);
    // create a streamline with given direction of the first panel
    // all subsequent panels are aligned with the given free-stream flow direction
    Streamline(Vector A, double xInfinity, int nSegment,
	       Vector startDirection, Vector freeStreamDirection);
    ~Streamline();

    // return some information about the subdivision into lines
    Vector point(int n);	// n is ranging from 0 to (including) segmentN
    Vector dir(int n);          // n is ranging from 0 to (excluding) segmentN
    int segmentN();

    // recompute the point position of a streamline
    // from a given set of velocities (directions) at the points
    void align(Vector *vPoint);

    // compute the velocity induced at a given point in space
    // if the streamline has a vorticity associated (i.e it is a vortex filament)
    Vector inducedVortexVelocity(Vector p);

  private:

    int NSeg;                   // number of corner points
    Vector *P;            	// point coordinates
    Vector *v;                // normalized flow direction

};

// the induced velocity of a straight vortex line does not
// belong to this class, but just fits here
Vector vortexInducedVelocity (Vector P1, Vector P2, Vector P);

#endif