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

#include "spline.h"
#include <math.h>

Spline::Spline(int NumberOfPoints,
          double *y)
{
  np = NumberOfPoints;
  yp = new double[np];
  y2p = new double[np];
  for (int i=0; i<np; i++) {
    yp[i]=y[i];
    y2p[i]=0.0;
  };
  computeDerivatives();
}

Spline::~Spline()
{
  delete yp;
  delete y2p;
}

int Spline::length()
{
  return(np);
}

void Spline::computeDerivatives()
{
  double *gam = new double[np];
  double bet, cim1;
  int i;

  y2p[0] = 0.0;
  gam[0] = 0.0;
  bet=1.0;
  cim1=0.0;
  for (i=1; i<np-1; i++) {
    gam[i] = cim1/bet;
    bet = (4.0-gam[i])/6.0;
    y2p[i] = (yp[i-1] - 2.0*yp[i] + yp[i+1] - y2p[i-1]/6.0)/bet;
    cim1=1.0/6.0;
  }
  gam[np-1] = bet/6.0;
  y2p[np-1] = 0.0;
  for (i=np-2; i>=0; --i) {
    y2p[i]=y2p[i] - gam[i+1]*y2p[i+1];
  }
  delete gam;
}

double Spline::f(double t)
{
  int it;
  double dt;
  double A, B, C, D;
  double result;

  it = floor(t);
  dt = t - (double)it;

  result=0.0;
  if ((it>=0) && (it<np-1))
  {
    A = 1.0-dt;
    B = dt;
    C = (A*A*A-A)/6.0;
    D = (B*B*B-B)/6.0;
    result = A*yp[it]+B*yp[it+1]+C*y2p[it]+D*y2p[it+1];
  };
  return(result);
}

double Spline::dfdt(double t)
{
  int it;
  double dt;
  double A, B, C, D;
  double result;

  it = floor(t);
  dt = t - (double)it;

  result=0.0;
  if ((it>=0) && (it<np-1))
  {
    A = 1.0-dt;
    B = dt;
    C = (1.0-3.0*A*A)/6.0;
    D = (-1.0+3.0*B*B)/6.0;
    result = -yp[it] + yp[it+1] + C*y2p[it] + D*y2p[it+1];
  };
  return(result);
}

int Spline::findExtrema(double *te, int maxnum)
{
  int nfound = 0;
  double ya,yb,y2a,y2b;
  double d, t;

  for (int i=0; i<np-1; i++)
  {		// loop over all segments
    ya=yp[i];
    y2a=y2p[i];
    yb=yp[i+1];
    y2b=y2p[i+1];
    d= y2a*y2a + y2a*y2b + y2b*y2b - 6.0*((y2a-y2b)*(ya-yb));
    // positive d indicates root exists
    if (d>=0.0)  {
      t=(3.0*y2a-sqrt(3.0*d))/(3.0*(y2a-y2b));
      // check first root for validity
      if ((t>=0.0) && (t<=1.0) && (nfound<maxnum-1)) {
        te[nfound++]=(double)i+t;
      }
      t=(3.0*y2a+sqrt(3.0*d))/(3.0*(y2a-y2b));
      // check second root for validity
      if ((t>=0.0) && (t<=1.0) && (nfound<maxnum-1)) {
        te[nfound++]=(double)i+t;
      }
    }
  };

  return(nfound);
}

int Spline::findInverse(double y0, double *ti, int maxnum)
{
  int nfound = 0;
  int i;
  double ya,yb;
  double t, dt;

  ya=f(0.0);
  for (i=1; i<20*(np-1); i++)
  {
    t=0.05*(double)i;
    yb=f(t);
    if (((ya>=y0)&&(yb<=y0))||((ya<=y0)&&(yb>=y0)))
    {
      // found root
      // first guess - secant
      t-=(y0-yb)/(ya-yb)*0.05;
      // polish up
      do {
        dt=(y0-f(t))/dfdt(t);
        t+=dt;
      } while (fabs(dt)>=1e-6);
      if (nfound<maxnum-1) {
        ti[nfound++]=t;
      }
    }
    ya=yb;
  }
  return(nfound);
}

Spline2D::Spline2D(Spline *spx, Spline* spy)
{
  // check for deviating number of points -> error
  sx=spx;
  sy=spy;
  np=sx->length();
}

int Spline2D::length()
{
  return(np);
}

double Spline2D::fx(double t)
{
  return(sx->f(t));
}

double Spline2D::fy(double t)
{
  return(sy->f(t));
}
