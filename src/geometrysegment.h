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

#ifndef GEOMETRYSEGMENT_H
#define GEOMETRYSEGMENT_H

#include <QtGui>
#include <QtXml/QDomDocument>

#include <math.h>
#include "global.h"
#include "vector.h"
#include "geometrystation.h"

class GeometryStation;

class GeometrySegment
{

  public:

    GeometrySegment(
        GeometryStation *left,
        GeometryStation *right);
    ~GeometrySegment();

    int getSpanN();
    int getChordN();
    int getFlapN();
    int getFlapGroup();
    double getFlapLeft();
    double getFlapRight();

    GeometryStation *getLeftStation();
    GeometryStation *getRightStation();

    // the segments within a wing are bi-directionally chained
    // watch out - the edge segments have Null pointers
    GeometrySegment *getLeftNeighbour();
    GeometrySegment *getRightNeighbour();
    void setLeftNeighbour(GeometrySegment *left);
    void setRightNeighbour(GeometrySegment *right);

    void setSpanN(int n);
    void setChordN(int n);
    void setFlapN(int n);
    void setFlapGroup(int n);
    void setFlapLeft(double ch);
    void setFlapRight(double ch);

    // reports about segment properties
    // length of the projection of the quarter-chord line onto the y-z plane in mm
    double segmentSpan();
    // area of the segment in dm²
    double segmentArea();
    // geometrical neutral point x coordinate
    double segmentXC();
    // rotational angle about the x axis in radian
    double segmentDihedral();

    // write the segment geometry into an XML (DOM) document
    void writeGeometryXML(
        QDomDocument doc,
        QDomElement parent);

    // alternative constructor reading the station data
    // from an XML document
    GeometrySegment(
        GeometryStation *left,
        GeometryStation *right,
        QDomElement parent);

  private:

    int nspan;
    int nchord;
    int nflap;
    int flapgroup;
    double flapleft;
    double flapright;
    GeometryStation *leftStation;
    GeometryStation *rightStation;
    GeometrySegment *leftNeighbour;
    GeometrySegment *rightNeighbour;

};

#endif
