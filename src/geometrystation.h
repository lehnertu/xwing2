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

#ifndef GEOMETRYSTATION_H
#define GEOMETRYSTATION_H

#include <QtGui>
#include <QtXml/QDomDocument>

#include <math.h>
#include "global.h"
#include "vector.h"
#include "geometrysegment.h"

class GeometrySegment;

class GeometryStation
{
  
  public:

    GeometryStation();
    GeometryStation(GeometryStation* st);
    GeometryStation(double xi,
                    double yi,
                    double zi,
                    double ci,
                    double ai);
    
    // alternative constructor reading the station data
    // from an XML document
    GeometryStation(QDomElement parent);

    ~GeometryStation();

    double getX();
    double getY();
    double getZ();
    double getS();
    double getChord();
    double getAlfa();
    int getAirfoilIndex();

    void setX(double value);
    void setY(double value);
    void setZ(double value);
    void setS(double value);
    void setChord(double value);
    void setAlfa(double value);
    void setAirfoilIndex(int value);
    
    Vector getNosePoint();
    Vector getCh4Point();
    Vector getEndPoint();

    // The 2 unit vectors ChordUnitVec and UpUnitVec define the station plane.
    // If the segments right and left of the station have different dihedral angles
    // the angle of the station is defined as the median of the two.
    // In this case the airfoil is elongated vertically by a factor stretchZ
    // This factor has also to be included in the computation of the angle of
    // attack (chord direction vector) because the cross section airfoil
    // has to be rotated by a larger amount, so that the adjacent segments
    // obtain the angel of attack prescribed by the user
    double getDihedral();
    double getStretchZ();
    Vector getChordVector();
    Vector getUpVector();
    Vector getSpanVector();
    
    // the stations are connections of two segments
    // watch out - the edge stations have Null pointers
    GeometrySegment *getLeftSegment();
    GeometrySegment *getRightSegment();
    // when setting the segments also the properties like
    // dihedral angle and stretching factors are computed
    void setSegments(GeometrySegment *left,
                     GeometrySegment *right);

    // write the station geometry into an XML (DOM) document
    void writeGeometryXML(
        QDomDocument doc,
        QDomElement parent);

  private:

    // x, y, z are the coordinates of the nose point
    double x;
    double y;
    double z;
    double spanpos;	// span position in the wing
    double chord;	// in mm
    double alfa;	// in degree
    
    // The airfoil indexing starts with zero (default initialization):
    // Zero indicates the lack of a defined airfoil -> strake
    // Numbers from 1 up refer to the airfoil database present in the model class.
    // Watch out to decrease the index by one when indexing the airfoil database!
    int airfoilIndex;

    GeometrySegment *leftSegment;
    GeometrySegment *rightSegment;

};

#endif
