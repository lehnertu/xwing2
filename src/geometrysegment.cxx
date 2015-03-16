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

#include "geometrysegment.h"

// ******************************************************
// ************* Class Geometry Segment *****************
// ******************************************************

GeometrySegment::GeometrySegment(
        GeometryStation *left,
        GeometryStation *right)
{
  nspan = 1;
  nchord = 4;
  nflap = 1;
  flapgroup = 0;
  flapleft = 0.75;
  flapright = 0.75;
  leftStation = left;
  rightStation = right;
}

GeometrySegment::~GeometrySegment()
{
}

int GeometrySegment::getSpanN()
{
  return(nspan);
}

int GeometrySegment::getChordN()
{
  return(nchord);
}

int GeometrySegment::getFlapN()
{
  return(nflap);
}

int GeometrySegment::getFlapGroup()
{
  return(flapgroup);
}

double GeometrySegment::getFlapLeft()
{
  return(flapleft);
}

double GeometrySegment::getFlapRight()
{
  return(flapright);
}

GeometryStation* GeometrySegment::getLeftStation()
{
  return(leftStation);
}

GeometryStation* GeometrySegment::getRightStation()
{
  return(rightStation);
}

GeometrySegment* GeometrySegment::getLeftNeighbour()
{
  return(leftNeighbour);
}

GeometrySegment* GeometrySegment::getRightNeighbour()
{
  return(rightNeighbour);
}

void GeometrySegment::setLeftNeighbour(GeometrySegment *left)
{
  leftNeighbour = left;
}

void GeometrySegment::setRightNeighbour(GeometrySegment *right)
{
  rightNeighbour = right;
}

void GeometrySegment::setSpanN(int n)
{
  nspan = n;
}

void GeometrySegment::setChordN(int n)
{
  nchord = n;
}

void GeometrySegment::setFlapN(int n)
{
  nflap = n;
  if (nflap > nchord) nflap = nchord;
}

void GeometrySegment::setFlapGroup(int n)
{
  flapgroup = n;
}

void GeometrySegment::setFlapLeft(double ch)
{
  flapleft = ch;
  if (flapleft<0.0) flapleft = 0.0;
  if (flapleft>1.0) flapleft = 1.0;
}

void GeometrySegment::setFlapRight(double ch)
{
  flapright = ch;
  if (flapright<0.0) flapright = 0.0;
  if (flapright>1.0) flapright = 1.0;
}

double GeometrySegment::segmentSpan()
{
  Vector p1 = leftStation->getNosePoint();
  Vector p2 = rightStation->getNosePoint();
  Vector sp = p2-p1;
  sp.x=0.0;
  return(sp.norm());
}

double GeometrySegment::segmentArea()
{
  return(
    (leftStation->getChord() + rightStation->getChord())*
    0.00005*segmentSpan()
  );
}

// The dihedral is defined as the dihedral of the nose line.
// The quarter-chord point can only be computed after the
// dihedral direction of the station is known, which depends on the
// dihedral angle of the segments attached.
double GeometrySegment::segmentDihedral()
{
  Vector p1 = leftStation->getNosePoint();
  Vector p2 = rightStation->getNosePoint();
  Vector sp = p2-p1;
  sp.x=0.0;
  return(atan2(sp.z, sp.y));
}

double GeometrySegment::segmentXC()
{
  Vector p1 = leftStation->getCh4Point();
  Vector p2 = rightStation->getCh4Point();
  double ch1 = leftStation->getChord();
  double ch2 = rightStation->getChord();
  return( (ch1*(2.0*p1.x+p2.x) + ch2*(2.0*p2.x+p1.x)) / 60000.0 *
    segmentSpan() / segmentArea() );
}

void GeometrySegment::writeGeometryXML(
  QDomDocument doc,
  QDomElement parent)
{
  QDomElement segmentElement = doc.createElement("segment");
  segmentElement.setAttribute("nspan", nspan);
  segmentElement.setAttribute("nchord", nchord);
  segmentElement.setAttribute("nflap", nflap);
  segmentElement.setAttribute("flapgroup", flapgroup);
  segmentElement.setAttribute("flapleft", QString("%1").arg(flapleft,8,'f',4));
  segmentElement.setAttribute("flapright", QString("%1").arg(flapright,8,'f',4));
  parent.appendChild(segmentElement);
}

GeometrySegment::GeometrySegment(
  GeometryStation *left,
  GeometryStation *right,
  QDomElement elem)
{
  QString s;
  s = elem.attribute("nspan","1");
  nspan=s.toInt();
  s = elem.attribute("nchord","4");
  nchord=s.toInt();
  s = elem.attribute("nflap","1");
  nflap=s.toInt();
  s = elem.attribute("flapgroup","0");
  flapgroup=s.toInt();
  s = elem.attribute("flapleft","0.75");
  flapleft=s.toDouble();
  s = elem.attribute("flapright","0.75");
  flapright=s.toDouble();
  leftStation = left;
  rightStation = right;
}

