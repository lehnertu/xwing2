/************************************************************************/
/*                                                                      */
/*  written by  U.Lehnert                               08/2012         */
/*                                                                      */
/************************************************************************/

#ifndef WAKESTRIPE_H
#define WAKESTRIPE_H

#include <math.h>

#include "vector.h"
#include "streamline.h"
#include "flatpanel.h"

using namespace std;

// This class represents a stripe of a wake sheet
// which orgiginates from a straight baseline given by two starting points
// like a trailing edge segmant of a wing.
// Its outline is formed by two streamlines of the flow.
// It is modeled by a number of doublet panels.

class WakeStripe
{

  public:

    WakeStripe(Streamline *A, Streamline *B);
    ~WakeStripe();

    // reference the wake stripe to the wing
    void setWingIndex(int index);
    int getWingIndex();
    // the span position of the wing, where the wake starts
    void setWingSpanPos(double s);
    double getWingSpanPos();
    // the wing chord at the wake position
    void setWingChord(double c);
    double getWingChord();

    // report the start/end point of the wake streamlines
    // parameter is 1 or 2 indicating if A or B is requested
    Vector wakeStart(int i);
    Vector wakeInfinity(int i);

    // Return a point near the center of the line AB as a wake control point
    // and normal direction at this point
    Vector wakeCP();
    Vector wakeNormal();

    // compute the velocity induced at a given point in space
    // using the panel model
    Vector inducedVelocity(Vector P);
    double inducedPotential(Vector P);

  private:

    Streamline *filamentA, *filamentB;
    Vector startA, startB;
    Vector CP;		// control point 0.01 mm behind the center of the line AB
    Vector Normal;      // normal flow direction at CP
    int wingIndex;	// the wing by which the wake is spawn
    double wingSpan;	// the span position
    double wingChord;	// the wing chord
    int NP;		// number of panels
    FlatPanel* *panels;	// array of panels span between the streamlines

};

#endif