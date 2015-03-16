/************************************************************************/
/*                                                                      */
/*         XX      XX W           W  III  NN    N    GGGG               */
/*           XX  XX    W         W    I   N N   N   G                   */
/*             XX       W   W   W     I   N  N  N   G  GGG              */
/*           XX  XX      W W W W      I   N   N N   G    G              */
/*         XX      XX     W   W      III  N    NN    GGGG               */
/*                                                                      */
/*  Version 2.0                                    U.Lehnert  6/2011    */
/*                                                                      */
/************************************************************************/


/*

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include "streamline.h"

Streamline::Streamline(Vector A, double xInfinity, int nSegment)
{
  P = new Vector[nSegment+1];
  v = new Vector[nSegment+1];
  NSeg = nSegment;
  Vector dx = Vector((xInfinity - A.x)/NSeg, 0.0, 0.0);
  for (int i=0; i<=NSeg; i++)
  {
    P[i] = A + dx*(double)i;
    v[i] = Vector(1.0, 0.0, 0.0);
  };
  // debugging output
  // for (int i=0; i<=NSeg; i++)
  //   printf("(%g, %g, %g)\n",P[i].x,P[i].y,P[i].z);
}

Streamline::Streamline(Vector A, double xInfinity, int nSegment, Vector startDirection)
{
  P = new Vector[nSegment+1];
  v = new Vector[nSegment+1];
  NSeg = nSegment;
  P[0] = A;
  v[0] = startDirection/startDirection.norm();
  P[1] = A + startDirection;
  // projected length of the first segment
  double dx = dot(startDirection, Vector(1,0,0));
  // we have to find an enlargement factor such that the last segment ends at xInfinity
  // dx * N^f == xInf
  double xinf = xInfinity - A.x;
  double f = log(xinf/dx)/log(NSeg);
  for (int i=1; i<=NSeg; i++)
  {
    P[i] = P[0] + startDirection*exp(f*log(i));
    v[i] = v[0];
  };
  // debugging output
  // printf("Streamline :\n");
  // for (int i=0; i<=NSeg; i++)
  //   printf("(%g, %g, %g)\n",P[i].x,P[i].y,P[i].z);
  // printf("\n");
}

Streamline::Streamline(Vector A, double xInfinity, int nSegment,
	    Vector startDirection, Vector freeStreamDirection)
{
  P = new Vector[nSegment+1];
  v = new Vector[nSegment+1];
  NSeg = nSegment;
  P[0] = A;
  v[0] = startDirection/startDirection.norm();
  P[1] = A + startDirection;
  v[1] = freeStreamDirection/freeStreamDirection.norm();
  // projected length of the first segment
  double dx = dot(startDirection, Vector(1,0,0));
  double xinf = xInfinity - A.x;
  // normalize the freeStreamDirection to have the same x component as the startDirection
  Vector freeStream = freeStreamDirection * (startDirection.x / freeStreamDirection.x);
  // we have to find an enlargement factor such that the last segment ends at xInfinity
  // dx * N^f == xInf
  double f = log(xinf/dx)/log(NSeg);
  for (int i=2; i<=NSeg; i++)
  {
    P[i] = P[1] + freeStream*(exp(f*log(i))-1);
    v[i] = v[1];
  };
  // debugging output
  // printf("Streamline :\n");
  // for (int i=0; i<=NSeg; i++)
  //   printf("(%g, %g, %g)\n",P[i].x,P[i].y,P[i].z);
  // printf("\n");
}

Streamline::~Streamline()
{
  delete[] P;
  delete[] v;
}

Vector Streamline::point(int n)
{
  if ((n>=0) && (n<=NSeg))
    return P[n];
  else
    return P[0];
}

Vector Streamline::dir(int n)
{
  if ((n>=0) && (n<=NSeg))
    return v[n];
  else
    return v[0];
}

int Streamline::segmentN()
{
  return NSeg;
}

void Streamline::align(Vector *vPoint)
{
  v[0] = vPoint[0];
  for (int i=0; i<NSeg; i++)
  {
    Vector p1 = P[i];
    Vector p2 = P[i+1];
    Vector dir = vPoint[i] + vPoint[i+1];
    P[i+1] = P[i] + dir*(p2.x-p1.x)/dir.x;
    v[i+1] = vPoint[i+1];
  };
}

Vector Streamline::inducedVortexVelocity(Vector p)
{
  Vector vind = Vector(0.0, 0.0, 0.0);
  for (int i=1; i<=NSeg; i++)
    vind += vortexInducedVelocity(P[i-1],P[i],p);
  return vind;
}

Vector vortexInducedVelocity (Vector P1, Vector P2, Vector P)
{
  // Einheitsvektor der Rotationsrichtung
  Vector dI = P2 - P1;
  dI.normalize();
  // Vektor vom Start des Wirbels zum Punkt P
  Vector r = P - P1;
  // Winkel zwischen r und dI
  double cos_t1=dot(r,dI)/r.norm();
  r = P - P2;
  // Winkel zwischen r und dI
  double cos_t2=dot(r,dI)/r.norm();
  // Richtung der induzierte Geschwindigkeit nach Biot-Savart
  Vector v = cross(r,dI);
  // der Betrag des hier berechnetet Vektors v ist
  // der senkrechte Abstand des Punktes P von der Wirbellinie
  double h=v.norm();
  if (h>0.0001)
    v *= (cos_t2-cos_t1)/(4*M_PI*h*h);
  else
  {
    v = Vector(0,0,0);
    // issue a warning
    if (((cos_t2>0.0) && (cos_t1<0.0)) |
        ((cos_t2<0.0) && (cos_t1>0.0))) {
      // fatal_error(" singularity in BOUND_VORTEX()");
    };
  };
  return(v);
};
