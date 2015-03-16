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

#include "geometrywing.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>

// ******************************************************
// ************* Class Geometry Wing  *******************
// ******************************************************

GeometryWing::GeometryWing()
{
  initDefault(2);
  update();
}

GeometryWing::GeometryWing(int ns)
{
  initDefault(ns);
  update();
}

void GeometryWing::initDefault(int ns)
{
  name = QString("no name");
  NumberOfStations=ns;
  if (NumberOfStations < 2) NumberOfStations=2;
  stations = new GeometryStation*[NumberOfStations];
  for (int i=0; i<NumberOfStations; i++)
  {
    stations[i] = new GeometryStation(0.0, -500.0+1000.0*(double)i/(double)(ns-1),0.0,100.0,0.0);
  };
  NumberOfSegments=NumberOfStations-1;
  segments = new GeometrySegment*[NumberOfSegments];
  for (int i=0; i<NumberOfSegments; i++)
  {
    segments[i] = new GeometrySegment(stations[i], stations[i+1]);
  };
  leftClosure = TRUE;
  rightClosure = TRUE;
}

void GeometryWing::update()
{
  // create pointer for every segment to its left and right neighbours
  for (int i=0; i<NumberOfSegments; i++)
  {
    if (i==0)
    {
      segments[i]->setLeftNeighbour(NULL);
    } else {
      segments[i]->setLeftNeighbour(segments[i-1]);
    };
    if (i==NumberOfSegments-1)
    {
      segments[i]->setRightNeighbour(NULL);
    } else {
      segments[i]->setRightNeighbour(segments[i+1]);
    };
  };
  // create pointers for each station to the left and right adjacent segments
  stations[0]->setSegments(NULL, segments[0]);
  stations[NumberOfStations-1]->setSegments(segments[NumberOfStations-2], NULL);
  for (int i=1; i<NumberOfStations-1; i++)
    stations[i]->setSegments(segments[i-1], segments[i]);
  // set span positions for all stations
  double span = - wingSpan()/2.0;
  stations[0]->setS(span);
  for (int i=0; i<NumberOfSegments; i++)
  {
    span += segments[i]->segmentSpan();
    segments[i]->getRightStation()->setS(span);
  };
}

GeometryWing::~GeometryWing()
{
  for (int i=1; i<NumberOfStations; i++)
  {
    delete stations[i-1];
  };
  delete stations;
  for (int i=1; i<NumberOfSegments; i++)
  {
    delete segments[i-1];
  };
  delete segments;
}

QString GeometryWing::getName()
{
  return(name);
}

void GeometryWing::setName(QString s)
{
  name = s;
}

int GeometryWing::numberOfStations()
{
  return(NumberOfStations);
}

int GeometryWing::numberOfSegments()
{
  return(NumberOfSegments);
}

GeometryStation* GeometryWing::getStation(int n)
{
  if (n<1) return(stations[0]);
  if (n>NumberOfStations) return(stations[NumberOfStations-1]);
  return(stations[n-1]);
}

GeometrySegment* GeometryWing::getSegment(int n)
{
  if (n<1) return(segments[0]);
  if (n>NumberOfSegments) return(segments[NumberOfSegments-1]);
  return(segments[n-1]);
}

void GeometryWing::addStation(int n)
{
  // create a new station
  // which will have the logical number n
  // n must be larger than 1 (there would be no station left from it, otherwise)
  // n can at maximum be NumberOfStations (a new station between the last two)
  GeometryStation *newstation;
  if ((n > 1) && (n <= NumberOfStations))
  {
    newstation = new GeometryStation(
      0.5*(stations[n-2]->getX()+stations[n-1]->getX()),
      0.5*(stations[n-2]->getY()+stations[n-1]->getY()),
      0.5*(stations[n-2]->getZ()+stations[n-1]->getZ()),
      0.5*(stations[n-2]->getChord()+stations[n-1]->getChord()),
      0.5*(stations[n-2]->getAlfa()+stations[n-1]->getAlfa()));
    // recreate array of stations
    GeometryStation **newarray = new GeometryStation*[NumberOfStations+1];
    NumberOfStations+=1;
    for (int i=0; i<n-1; i++)
      newarray[i] = stations[i];
    newarray[n-1] = newstation;
    for (int i=n; i<NumberOfStations; i++)
      newarray[i] = stations[i-1];
    // replace the old with the new array
    delete stations;
    stations = newarray;
    // create new segments adjacent to the new station
    GeometrySegment *oldsegment = segments[n-2];
    GeometrySegment *leftsegment = new GeometrySegment(stations[n-2], stations[n-1]);
    leftsegment->setSpanN(oldsegment->getSpanN());
    leftsegment->setChordN(oldsegment->getChordN());
    leftsegment->setFlapN(oldsegment->getFlapN());
    leftsegment->setFlapGroup(oldsegment->getFlapGroup());
    leftsegment->setFlapLeft(oldsegment->getFlapLeft());
    leftsegment->setFlapRight(0.5*(oldsegment->getFlapLeft() + oldsegment->getFlapRight()));
    GeometrySegment *rightsegment = new GeometrySegment(stations[n-1], stations[n]);
    rightsegment->setSpanN(oldsegment->getSpanN());
    rightsegment->setChordN(oldsegment->getChordN());
    rightsegment->setFlapN(oldsegment->getFlapN());
    rightsegment->setFlapGroup(oldsegment->getFlapGroup());
    rightsegment->setFlapLeft(0.5*(oldsegment->getFlapLeft() + oldsegment->getFlapRight()));
    rightsegment->setFlapRight(oldsegment->getFlapRight());
    delete oldsegment;
    // recreate array of segments
    NumberOfSegments = NumberOfStations-1;
    GeometrySegment **newsegarray = new GeometrySegment*[NumberOfSegments];
    for (int i=0; i<n-2; i++)
      newsegarray[i] = segments[i];
    newsegarray[n-2] = leftsegment;
    newsegarray[n-1] = rightsegment;
    for (int i=n; i<NumberOfSegments; i++)
      newsegarray[i] = segments[i-1];
    // replace the old with the new array
    delete segments;
    segments = newsegarray;
    update();
  }
}

void GeometryWing::deleteStation(int n)
{
  if ((NumberOfStations>2) && (n>=2) && (n<NumberOfStations))
  {
    // unite the two segments surrounding the station to be deleted
    GeometrySegment *newsegment = new GeometrySegment(stations[n-2], stations[n]);
    newsegment->setSpanN(segments[n-2]->getSpanN());
    newsegment->setChordN(segments[n-2]->getChordN());
    newsegment->setFlapN(segments[n-2]->getFlapN());
    newsegment->setFlapGroup(segments[n-2]->getFlapGroup());
    newsegment->setFlapLeft(segments[n-2]->getFlapLeft());
    newsegment->setFlapRight(segments[n-1]->getFlapRight());
    delete segments[n-2];
    delete segments[n-1];
    segments[n-2] = newsegment;
    for (int i=n; i<NumberOfSegments; i++)
      segments[i-1] = segments[i];
    NumberOfSegments-=1;
    // delete the station object and reorganize the list
    delete stations[n-1];
    for (int i=n; i<NumberOfStations; i++)
      stations[i-1] = stations[i];
    NumberOfStations-=1;
    update();
  };
}

bool GeometryWing::getLeftClosure()
{
  return(leftClosure);
}

void GeometryWing::setLeftClosure(bool state)
{
  leftClosure = state;
}

bool GeometryWing::getRightClosure()
{
  return(rightClosure);
}

void GeometryWing::setRightClosure(bool state)
{
  rightClosure = state;
}

double GeometryWing::wingSpan()
{
  double span=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    span += segments[i]->segmentSpan();
  };
  return(span);
}

double GeometryWing::wingArea()
{
  double area=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    area += segments[i]->segmentArea();
  };
  return(area);
}

double GeometryWing::wingLambda()
{
  double s = wingSpan();
  double A = wingArea();
  return(0.0001*s*s/A);
}

double GeometryWing::wingLiftingArea()
{
  double area=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    double cosdih = cos(segments[i]->segmentDihedral());
    area += segments[i]->segmentArea() * cosdih*cosdih;
  };
  return(area);
}

double GeometryWing::wingXC()
{
  double area=0.0;
  double sum=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    double cosdih = cos(segments[i]->segmentDihedral());
    area += segments[i]->segmentArea() * cosdih*cosdih;
    sum += segments[i]->segmentArea() * cosdih*cosdih * segments[i]->segmentXC();
  };
  return(sum/area);
}

double GeometryWing::wingLmu()
{
  double area=0.0;
  double sum=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    double ch1 = segments[i]->getLeftStation()->getChord();
    double ch2 = segments[i]->getRightStation()->getChord();
    double span = segments[i]->segmentSpan();
    double cosdih = cos(segments[i]->segmentDihedral());
    area += segments[i]->segmentArea() * cosdih*cosdih;
    sum += 0.5*(ch1*ch1+ch2*ch2)*span * cosdih*cosdih;
  };
  return(0.0001*sum/area);
}

double GeometryWing::wingSideArea()
{
  double area=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    double sindih = sin(segments[i]->segmentDihedral());
    area += segments[i]->segmentArea() * sindih*sindih;
  };
  return(area);
}

double GeometryWing::wingXS()
{
  double area=0.0;
  double sum=0.0;
  for (int i=0; i<NumberOfSegments; i++)
  {
    double sindih = sin(segments[i]->segmentDihedral());
    area += segments[i]->segmentArea() * sindih*sindih;
    sum += segments[i]->segmentArea() * sindih*sindih * segments[i]->segmentXC();
  };
  return(sum/area);
}

void GeometryWing::sourceVTK(
      vtkPolyData *polyData,
      int highlightSegment,
      int highlightStation)
{
  Vector P;
  int lineIndex=0; // for keeping track of the lines already created
  // first set all point coordinates
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(2*NumberOfStations);
  for (int st=0; st<NumberOfStations; st++)
    {
      P = stations[st]->getNosePoint();
      pts->InsertPoint(2*st, P.x, P.y, P.z);
      P = stations[st]->getEndPoint();
      pts->InsertPoint(2*st+1, P.x, P.y, P.z);
    };
  polyData->SetPoints(pts);
  pts->Delete();
  // all lines that connect the given points
  vtkCellArray *lines = vtkCellArray::New();
  vtkCellArray *polys = vtkCellArray::New();
  // define an array of colors for the lines and faces
  vtkSmartPointer<vtkUnsignedCharArray> lineColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  lineColors->SetName("Colors");
  lineColors->SetNumberOfComponents(3);
  lineColors->SetNumberOfTuples(2*NumberOfSegments+2*NumberOfStations);
  // define the colors
  // unsigned char outlineColor[3] = {0,255,0};
  // unsigned char highlightColor[3] = {255,0,0};
  // double color[3] = {0, 0, 0};
  // draw leading edge and trailing edge lines
  for (int seg=0; seg<NumberOfSegments; seg++)
    {
      bool isHighlighted = FALSE;
      if ((highlighting == SegmentHighlight) && (seg+1 == highlightSegment))
        isHighlighted = TRUE;
      if (highlighting == WingHighlight) isHighlighted = TRUE;
      // nose line
      lines->InsertNextCell(2);
      lines->InsertCellPoint(lineIndex);
      lines->InsertCellPoint(lineIndex+2);
      if (isHighlighted) {
        lineColors->InsertTuple3(lineIndex,255,0,255);
      } else {
        lineColors->InsertTuple3(lineIndex,200,200,200);
      };
      lineIndex++;
      // trailing edge line
      lines->InsertNextCell(2);
      lines->InsertCellPoint(lineIndex);
      lines->InsertCellPoint(lineIndex+2);
      if (isHighlighted) {
        lineColors->InsertTuple3(lineIndex,255,0,255);
      } else {
        lineColors->InsertTuple3(lineIndex,200,200,200);
      };
      lineIndex++;
    };
   // draw station lines
   for (int st=0; st<NumberOfStations; st++)
    {
      bool isHighlighted = FALSE;
      if ((highlighting == StationHighlight) && (st+1 == highlightStation))
        isHighlighted = TRUE;
      if ((highlighting == SegmentHighlight) && (st+1 == highlightSegment))
        isHighlighted = TRUE;
      if ((highlighting == SegmentHighlight) && (st+1 == highlightSegment+1))
        isHighlighted = TRUE;
      if (highlighting == WingHighlight) isHighlighted = TRUE;
      lines->InsertNextCell(2);
      lines->InsertCellPoint(2*st);
      lines->InsertCellPoint(2*st+1);
      if (isHighlighted) {
        lineColors->InsertTuple3(lineIndex,255,0,255);
      } else {
        lineColors->InsertTuple3(lineIndex,200,200,200);
      };
      lineIndex++;
    };
   // draw segment chord-plane polygons
   for (int seg=0; seg<NumberOfSegments; seg++)
    {
      bool isHighlighted = FALSE;
      if ((highlighting == SegmentHighlight) && (seg+1 == highlightSegment))
        isHighlighted = TRUE;
      if (highlighting == WingHighlight) isHighlighted = TRUE;
      vtkQuad *quad = vtkQuad::New();
      quad->GetPointIds()->SetId(0,2*seg);
      quad->GetPointIds()->SetId(1,2*seg+2);
      quad->GetPointIds()->SetId(2,2*seg+3);
      quad->GetPointIds()->SetId(3,2*seg+1);
      // polys->InsertNextCell(quad->GetCellType(),quad->GetPointIds());
      polys->InsertNextCell(quad);
      quad->Delete();
      if (isHighlighted) {
        lineColors->InsertTuple3(lineIndex,255,0,255);
      } else {
        lineColors->InsertTuple3(lineIndex,200,200,200);
      };
      lineIndex++;
    };
  polyData->SetLines(lines);
  polyData->SetPolys(polys);
  // attach the line color dataset
  polyData->GetCellData()->AddArray(lineColors);
  lines->Delete();
  polys->Delete();
  // lineColors->Delete();
}

void GeometryWing::setHighlighting (geometryHighlightMode mode)
{
  highlighting = mode;
}

void GeometryWing::writeGeometryXML(
        QDomDocument doc,
        QDomElement parent)
{
  QDomElement wingElement = doc.createElement("wing");
  wingElement.setAttribute("name", this->getName());
  if (leftClosure)
    wingElement.setAttribute("leftClosure", "True");
  else
    wingElement.setAttribute("leftClosure", "False");
  if (rightClosure)
    wingElement.setAttribute("rightClosure", "True");
  else
    wingElement.setAttribute("rightClosure", "False");
  for (int i=0; i<NumberOfStations; i++)
  {
    stations[i]->writeGeometryXML(doc, wingElement);
  };
  for (int i=0; i<NumberOfSegments; i++)
  {
    segments[i]->writeGeometryXML(doc, wingElement);
  };
  parent.appendChild(wingElement);
}

GeometryWing::GeometryWing(QDomElement parent)
{
  QDomNode node;
  QDomElement elem;
  QString s;
  GeometryStation *newstation;
  GeometrySegment *newsegment;
  this->setName(parent.attribute("name","no name"));
  s = parent.attribute("leftClosure","True");
  this->setLeftClosure(s=="True" || s=="TRUE" || s=="true");
  s = parent.attribute("rightClosure","True");
  this->setRightClosure(s=="True" || s=="TRUE" || s=="true");
  // browse nodes for station entries and count them
  node=parent.firstChild();
  NumberOfStations=0;
  while( !node.isNull() )
  {
    QDomElement elem=node.toElement();
    if( !elem.isNull() )
    {
      if( elem.tagName() == "station" ) NumberOfStations++;
    }
    node=node.nextSibling();
  };
  if (NumberOfStations<2)
  {
    Globals::MainTextDisplay->append(
        QString("\n")+
        QString("error reading wing data\n")+
        QString("less than 2 valid stations found - using default\n\n"));
    initDefault(2);
    return;
  };
  // construct the stations
  stations = new GeometryStation*[NumberOfStations];
  // browse once more and read stations
  node=parent.firstChild();
  int ns=0;
  while( !node.isNull() && (ns<NumberOfStations) )
  {
    elem=node.toElement();
    if( !elem.isNull() )
    {
      if( elem.tagName() == "station" )
      {
        newstation = new GeometryStation(elem);
        stations[ns++] = newstation;
      }
    }
    node=node.nextSibling();
  };
  // all stations not properly read from file are initialized as default
  while( ns<NumberOfStations )
  {
    Globals::MainTextDisplay->append(
        QString("warning : initializing station from default\n"));
    newstation = new GeometryStation();
    stations[ns++] = newstation;
  };
  // construct the segments
  NumberOfSegments=NumberOfStations-1;
  segments = new GeometrySegment*[NumberOfSegments];
  // browse once more and read segments
  node=parent.firstChild();
  ns=0;
  while( !node.isNull() && (ns<NumberOfSegments) )
  {
    elem=node.toElement();
    if( !elem.isNull() )
    {
      if( elem.tagName() == "segment" )
      {
        newsegment = new GeometrySegment(stations[ns],stations[ns+1],elem);
        segments[ns++] = newsegment;
      }
    }
    node=node.nextSibling();
  };
  // all segments not properly read from file are initialized as default
  while( ns<NumberOfSegments )
  {
    Globals::MainTextDisplay->append(
        QString("warning : initializing segment from default\n"));
    newsegment= new GeometrySegment(stations[ns],stations[ns+1]);
    segments[ns++] = newsegment;
  };
  update();
}

