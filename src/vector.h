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

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include "global.h"

# define over4Pi 0.25/M_PI

using namespace std;

class Vector
{

  public:

    Vector();

    // create a verctor from cartesian coordinates
    Vector(double x0, double y0, double z0);

    // create a vector from a length 3 array of double values
    Vector(double *p0);

    double x;
    double y;
    double z;

    // sum and difference of two vectors
    Vector operator+ (Vector v2);
    Vector& operator+= (Vector v2);
    Vector operator- (Vector v2);
    Vector& operator-= (Vector v2);

    // negative vector
    Vector operator- ();

    // multiplication of a vector with a factor
    Vector operator* (double faktor);
    Vector& operator*= (double faktor);
    Vector operator/ (double faktor);
    Vector& operator/= (double faktor);

    // absolute value of a vector
    double norm();

    // square of the norm
    double sqnorm();

    // make a vector unit length
    void normalize();

  private:

};

// dot product of two vectors
double dot(Vector a, Vector b);

// cross product of two vectors
Vector cross(Vector a, Vector b);

// rotate Vector v about an axis given by Vector n by an angle alpha [radian]
Vector rotate(Vector v, Vector n, double alpha);

#endif