/************************************************************************/
/*                                                                      */
/*  written by  U.Lehnert                               08/2012         */
/*                                                                      */
/************************************************************************/

#ifndef FLATPANEL_H
#define FLATPANEL_H

#include <math.h>

#include "vector.h"

using namespace std;

// *** FlatPanel ***
// flat panel given by 3-4 corner points (code will work for more,
// but flatness, convexity, area computation etc. my be a problem)

class FlatPanel
{

  public:

    FlatPanel(FlatPanel *p);	// copy constructor
    FlatPanel(int n, Vector* pp);
    FlatPanel(Vector A, Vector B, Vector C);
    FlatPanel(Vector A, Vector B, Vector C, Vector D);
    ~FlatPanel();

    // return some information about the panel
    int panelNP();
    Vector panelPoint(int i);
    Vector panelCenter();
    Vector panelNormal();
    Vector panelXi();
    Vector panelEta();
    int panelIsPlanar();
    double panelArea();
    Vector panelNue(int i);

    // This function is called once for a given point P.
    // It computes all panel parameters related to the
    // panel influence coefficients (PIC) in particular
    // the H(M,N,K) and F(M,N,K) integrals.
    // It must be called before any panel-induced
    // velocities or potentials can be computed
    void ComputePIC(Vector P);

    // some report methods provided for testing
    int panelPointIsInterior();
    double panelPointh();
    double panelPointa(int i);
    double panelPointg(int i);

    double panelH111();
    double panelH113(); // presently returns abs(h)*H113
    double panelF111(int i);

    // These functions compute the induced velocity and potential
    // assuming a constant dipole strength on the panel.
    // It applies the vortex ring equivalent for velocity computation.
    // The velocity is singular right on the panel boundaries
    // in case the function is called for points within 10^-6 of
    // the panel boundaries zero is returned.
    // The potential jump occurs right at the panel plane
    // so, self-induction control points should be moved clearly to the outside
    double ConstDoubletPotential();
    Vector ConstDoubletVelocity();

    // This function computes the induced velocity
    // assuming a constant source strength on the panel.
    // The normal velocity jumps across the panel.
    // For positive source strength the potential is negative
    // and the flow is directed away from the panel.
    // The velocity jump occurs right at the panel plane
    // so, self-induction control points should be moved clearly to the outside
    double ConstSourcePotential();
    Vector ConstSourceVelocity();

  private:

    int NP;                  // number of corner points
    Vector* CP;              // corner point coordinates
    Vector Center;           // panel center
    double Area;             // panel area
    int isPlanar;            // if panel really is co-planar

    // normal vector pointing to the outside of the body
    // seen from the outside the corner points run in
    // counter-clockwise direction
    Vector Normal;

    // unit Vectors defining the panel coordinate system
    // Xi has zero y component
    // or Eta has zero x component
    // (Xi, Eta, Normal) is a right-handed orthogonal coordinate system
    Vector Xi, Eta;

    Vector* ds;              // unit vectors along the panel edges
    Vector* nue;             // unit vectors perpendicular to the edges
                             // pointig out of te panel in the panel plane

    // all the initalization that has to be done by the constructors
    // are compiled into this methods
    void init();

    // the following variables are related to a given
    // influence point P and are computed by ComputePIC(P)
    Vector* R;		     // distance vector
    double* r;		     // absolut distance
    Vector* vv;              // R x ds
    double h, ha;	     // elevation (abs) from panel plane
    double* a;		     // projected distance from edge
    double* aa;		     // absolute value of a
    double* g;               // disctance from edge
    double* l1;              // projections of P along a panel edge
    double* l2;
    int isInterior;          // if projection of P onto the panel plane lies inside the panel

    // values of the H(M,N,K) and F(M,N,K) integrals
    // computed by ComputePIC(P) for a given influence point P
    double hH113;
    double* F111;
    double H111;

    // const. strength ssingularity distribution over the panel area
    // induced velicities and potentials
    Vector vConstDoublet;
    double fConstDoublet;
    Vector vConstSource;
    double fConstSource;

};

#endif