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

#include "geometrystation.h"

// ******************************************************
// ************* Class Geometry Station *****************
// ******************************************************

GeometryStation::GeometryStation()
{
  x=0.0;
  y=0.0;
  z=0.0;
  spanpos=0.0;
  chord=100.0;
  alfa=0.0;
  airfoilIndex=0;
}

GeometryStation::GeometryStation(GeometryStation* st)
{
  x=st->getX();
  y=st->getY();
  z=st->getZ();
  spanpos=0.0;
  chord=st->getChord();
  alfa=st->getAlfa();
  airfoilIndex=st->getAirfoilIndex();
}

GeometryStation::~GeometryStation()
{
}

GeometryStation::GeometryStation(double xi,
                                 double yi,
                                 double zi,
                                 double ci,
                                 double ai)
{
  x=xi;
  y=yi;
  z=zi;
  spanpos=0.0;
  chord=ci;
  alfa=ai;
  airfoilIndex=0;
}

double GeometryStation::getX()
{
  return(x);
}

double GeometryStation::getY()
{
  // qDebug() << " y= " << y << " at " << this ;
  return(y);
}

double GeometryStation::getZ()
{
  return(z);
}

double GeometryStation::getS()
{
  return(spanpos);
}

double GeometryStation::getChord()
{
  return(chord);
}

double GeometryStation::getAlfa()
{
  return(alfa);
}

int GeometryStation::getAirfoilIndex()
{
  return(airfoilIndex);
}

void GeometryStation::setX(double value)
{
  x = value;
}

void GeometryStation::setY(double value)
{
  y = value;
}

void GeometryStation::setZ(double value)
{
  z = value;
}

void GeometryStation::setS(double span)
{
  spanpos = span;
}

void GeometryStation::setChord(double value)
{
  chord = value;
}

void GeometryStation::setAlfa(double value)
{
  alfa = value;
}

void GeometryStation::setAirfoilIndex(int value)
{
  airfoilIndex = value;
}

Vector GeometryStation::getNosePoint()
{
  Vector temp;
  temp.x = x;
  temp.y = y;
  temp.z = z;
  return (temp);
}

Vector GeometryStation::getCh4Point()
{
  Vector p1 = getNosePoint();
  Vector ec = getChordVector();
  return(p1 + ec*0.25*chord);
}

Vector GeometryStation::getEndPoint()
{
  Vector p1 = getNosePoint();
  Vector ec = getChordVector();
  return(p1 + ec*chord);
}

double GeometryStation::getDihedral()
{
  double quotient = 0.0;
  double dihedral = 0.0;
  if (leftSegment != NULL)
  {
    dihedral += leftSegment->segmentDihedral();
    quotient += 1.0;
  };
  if (rightSegment != NULL)
  {
    dihedral += rightSegment->segmentDihedral();
    quotient += 1.0;
  };
  if (quotient != 0.0) dihedral /= quotient;
  return(dihedral);
}

double GeometryStation::getStretchZ()
{
  // bei zunehmendem Anstellwinkel muß der Dehnungsfaktor wieder gegen 1 gehen
  double diff = 0.0;
  if ((leftSegment != NULL) && (rightSegment != NULL))
    diff = fabs(leftSegment->segmentDihedral()-rightSegment->segmentDihedral());
  return(1.0/cos(0.5*diff*cos(alfa*M_PI/180.0)));
}

Vector GeometryStation::getChordVector()
{
  Vector ex = Vector(1.0, 0.0, 0.0);
  Vector es = getSpanVector();
  Vector chordUnitVec = rotate(ex, es, alfa*getStretchZ()*M_PI/180.0);
  return(chordUnitVec);
}

Vector GeometryStation::getUpVector()
{
  Vector es = getSpanVector();
  Vector ec = getChordVector();
  return(cross(ec, es));
}

Vector GeometryStation::getSpanVector()
{
  Vector ex = Vector(1.0, 0.0, 0.0);
  Vector ey = Vector(0.0, 1.0, 0.0);
  Vector spanUnitVec = rotate(ey, ex, getDihedral());
  return(spanUnitVec);
}

void GeometryStation::writeGeometryXML(
        QDomDocument doc,
        QDomElement parent)
{
  QDomElement stationElement = doc.createElement("station");
  stationElement.setAttribute("x", QString("%1").arg(x,10,'f',3));
  stationElement.setAttribute("y", QString("%1").arg(y,10,'f',3));
  stationElement.setAttribute("z", QString("%1").arg(z,10,'f',3));
  stationElement.setAttribute("chord", QString("%1").arg(chord,10,'f',3));
  stationElement.setAttribute("alfa", QString("%1").arg(alfa,8,'f',3));
  stationElement.setAttribute("airfoilindex", airfoilIndex);
  parent.appendChild(stationElement);
}

GeometrySegment* GeometryStation::getLeftSegment()
{
  return(leftSegment);
}

GeometrySegment* GeometryStation::getRightSegment()
{
  return(rightSegment);
}

void GeometryStation::setSegments(GeometrySegment *left,
                                  GeometrySegment *right)
{
  leftSegment = left;
  rightSegment = right;
}

GeometryStation::GeometryStation(QDomElement elem)
{
  QString s = elem.attribute("x","0.0");
  x=s.toDouble();
  s = elem.attribute("y","0.0");
  y=s.toDouble();
  s = elem.attribute("z","0.0");
  z=s.toDouble();
  s = elem.attribute("chord","100.0");
  chord=s.toDouble();
  s = elem.attribute("alfa","0.0");
  alfa=s.toDouble();
  s = elem.attribute("airfoilindex","0");
  airfoilIndex=s.toInt();
}

