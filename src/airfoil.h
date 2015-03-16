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

#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <QtGui>
#include <QtXml/QDomDocument>

#include "spline.h"
#include "global.h"
#include "vector.h"
#include <math.h>

class Airfoil;
class Vector;

class Camberline
{

public:

  Camberline(int NumberOfPanels);
  Camberline(int NumberOfPanels,
             // specified paneling for given flap hinge
             int NumberFlapPanels,
             double flapPosition);

  // Interpolation of two camberlines
  // Both must have equal number of panels.
  // Both should have similar paneling parameters.
  // If the camberlines are mapped to Vector coordinates, also the Vectors are interpolated
  // If the interpolation fails, a camberline with a single straight panel is created
  Camberline(Camberline *cl1, double f1,
             Camberline *cl2, double f2);

  ~Camberline();

  // return number of points/panels
  int numberPoints();
  int numberPanels();

  double pointX(int i);
  double pointY(int i);

  // pointVec() requires that mapToVectors() has been executed
  Vector pointVec(int i);

  // the x-position of the camber line panels are computed during initialization
  // the y-positions have to be set using one of these methode
  void setY(int i, double y);
  void computePoints(Airfoil *af);

  // compute the aerodynamic coefficients of a camberline
  // with the given paneling
  void computeCoeff(
    double *cl0,        // lift at zero angla of attack
    double *alfa0,      // angle of attack for zero lift in radian
    double *cm0,        // momentum
    double *asf,        // angle of attack for suction-free entrance
    double *clsf        // lift coefficient for suction-free entrance
  );

  // The camberline points are mapped into real space.
  // The points (0,0) and (1,0) are mapped to the vectors nose and end
  // up is a vector indicating the y direction of the camberline.
  // up is expected to be perpendicular to the nose-end line.
  // The up direction is additionally scaled by a factor stretch.
  void mapToVectors(Vector nose, Vector end, Vector up, double stretch);

  // the direction of an emanating wake, i.e. the extension of the camber line
  // the length of the vector equals the length of the last camberline segment
  Vector wakeDirection();

private:

  int npan;
  double *xs;
  double *ys;
  double *aoa;  // panel angle of attack
  Vector *p;
};

class Outline
{

public:

  // NumberOfPanels is the number of line segments forming the outline.
  // The number of points is larger by 1.
  Outline(int NumberOfPanels);

  // Interpolation of two outlines
  // Both must have equal number of panels.
  // Both should have similar paneling parameters.
  // If the outlines are mapped to Vector coordinates, also the Vectors are interpolated
  // If the interpolation fails, an outline with two straight panels is created
  Outline(Outline *ol1, double f1,
	  Outline *ol2, double f2);

  ~Outline();

  // return number of points/panels
  int numberPoints();
  int numberPanels();

  // index i is running from 0 to NumberOfPanels
  double pointX(int i);
  double pointY(int i);

  // the x-position of the outline panels are computed during initialization
  // the y-positions have to be set using one of these methode
  void setY(int i, double y);
  void computePoints(Airfoil *af);

  // The outline points are mapped into real space.
  // The points (0,0) and (1,0) are mapped to the vectors nose and end
  // up is a vector indicating the y direction of the section.
  // up is expected to be perpendicular to the nose-end line.
  // The up direction is additionally scaled by a factor stretch.
  void mapToVectors(Vector nose, Vector end, Vector up, double stretch);

  // the following methods require that mapToVectors() has been executed
  // coordinates of an outline point
  Vector pointVec(int i);

  // the direction of an emanating wake, i.e. the extension of the camber line
  // bisecting the trailing edge angle.
  // the length of the vector equals the length of the trailing neighbouring outline segments
  Vector wakeDirection();

private:

  int npan;
  double *xs;
  double *ys;
  Vector *p;

};

class Airfoil
{

public:

  Airfoil(int NumberOfPoints,
         double *x,
         double *y);

  // load airfoil coordinates from a standard *.dat file
  Airfoil(const QString &fileName);

  // alternative constructor reading the airfoil data
  // from an XML document
  Airfoil(QDomElement af);

  ~Airfoil();

  // these 3 splines only exist if valid==TRUE
  Spline2D *outline;
  Spline *splx;
  Spline *sply;

  // return number of points
  int numPoints();

  // return airfoil surface coordinate points
  double pointX(int i);
  double pointY(int i);

  QString getName();
  bool isValid();

  // write the airfoil geometry into an XML (DOM) document
  void writeGeometryXML(
      QDomDocument doc,
      QDomElement parent);

  // sanity check for the airfoil points
  // leading edge and trailing edge should be where expected
  // circumference and area should be reasonable
  bool doSanityCheck();

  // keep track of index in the database
  void setIndex(int idx);
  int getIndex();

private:

  QString *name;
  bool valid;
  int np;
  double *xp;
  double *yp;
  double circumference;
  double area;
  bool counterclockwise;

  // compute circumference and area of the polygon
  void computeOutlineProperties();

  // normalize airfoil
  // point with lowest x (on spline) -> (0, 0)
  // first and last points both having x=1
  // first and last points symmetrical to (1, 0)
  void normalize();

  // position in database
  // used for checking consistency during data reading
  int databaseIndex;

};

#endif
