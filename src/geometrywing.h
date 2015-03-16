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

#ifndef GEOMETRYWING_H
#define GEOMETRYWING_H

#include <QtGui>
#include <QtXml/QDomDocument>

// these defines are necessary for VTK only when building with qmake
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include <vtkPolyData.h>

#include <math.h>
#include "global.h"
#include "vector.h"
#include "geometrystation.h"
#include "geometrysegment.h"

class GeometryStation;
class GeometrySegment;

class GeometryWing
{

  public:

    GeometryWing();
    GeometryWing(int ns);       // constructor for a wing with a given number of stations
    ~GeometryWing();

  private:

    // supplementary routine for default initialization
    // used in the constructors
    void initDefault(int ns);

  public:

    // the individual segments and stations of a wing are linked
    // this is done here
    // update should be called whenever something has changed
    void update();

    QString getName();
    void setName(QString s);

    int numberOfStations();
    int numberOfSegments();

    // n is the logical number (1...n-1) of the new station
    void addStation(int n);
    void deleteStation(int n);

    bool getLeftClosure();
    void setLeftClosure(bool state);
    bool getRightClosure();
    void setRightClosure(bool state);

    // reports about wing properties

    // length of the projection of the quarter-chord line onto the y-z plane in mm
    double wingSpan();

    // total wing area in dm²
    double wingArea();

    // aspect ratio
    double wingLambda();

    // effective lifting wing area in dm²
    double wingLiftingArea();

    // reference chord lµ in mm for lifting area
    double wingLmu();

    // geometrical neutral point
    double wingXC();

    // effective side wing area in dm²
    double wingSideArea();

    // geometrical neutral point of side area
    double wingXS();

    // how to highlight parts of the model when rendering
    void setHighlighting (geometryHighlightMode mode);

    // generate source data (vtkPolyData)
    // for visualization using VTK
    void sourceVTK(
      vtkPolyData *polyData,
      int highlightSegment = 0,   // logic segment number for highlighting
      int highlightStation = 0);   // logic station number for highlighting

    // return a pointer to the station object
    // n is the logical number running from 1 to NumberOfStations
    GeometryStation* getStation(int n);
    // n is the logical number running from 1 to NumberOfSegments
    GeometrySegment* getSegment(int n);

    // write the wing geometry into an XML (DOM) document
    void writeGeometryXML(
                QDomDocument doc,
                QDomElement parent);

    // alternative constructor reading the wing
    // from an XML document
    GeometryWing(QDomElement parent);

  private:

    QString name;

    int NumberOfStations;
    GeometryStation **stations; // array of pointers to stations

    int NumberOfSegments;
    GeometrySegment **segments; // array of pointers to segments

    geometryHighlightMode highlighting;

    bool leftClosure;
    bool rightClosure;

};

#endif
