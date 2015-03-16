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


#ifndef SPLINE_H
#define SPLINE_H

class Spline
{

public:

  Spline(int NumberOfPoints,
         double *y);

  ~Spline();

  int length();

  // this function returns the value of the spline function at parameter t
  double f(double t);

  // this function returns the value of the
  // first derivative of the spline function with respect to t
  double dfdt(double t);

  // this function computes the t values at which extrema of the function occur
  // it finds all zeros of the first derivative
  // return value is the number of extrema found
  // the search is stopped when either maxnum or all existing extrema are found
  int findExtrema(double *te, int maxnum);

  // this function computes the t values at which
  // the spline function returns y0
  // return value is the number of positions found
  // the search is stopped when either maxnum or all existing positions are found
  int findInverse(double y0, double *ti, int maxnum);

private:

  int np;               // number of points
  double *yp;           // point coordinates
  double *y2p;          // second derivatives

  // computes the second derivatives
  // so that the first derivatives are continuous across the interval boundaries
  // "natural" spline - i.e. second derivative is zero at the ends
  void computeDerivatives();

};

class Spline2D
{

  public:

    Spline2D(Spline *spx, Spline* spy);
    int length();
    double fx(double t);
    double fy(double t);

  private:

    int np;               // number of points
    Spline *sx;
    Spline *sy;

};

#endif
