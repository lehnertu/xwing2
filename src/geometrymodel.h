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

#ifndef GEOMETRYMODEL_H
#define GEOMETRYMODEL_H

#include <QtGui>
#include <QtXml/QDomDocument>

// these defines are necessary for VTK only when building with qmake
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include <vtkAppendPolyData.h>

#include <math.h>
#include "global.h"
#include "vector.h"
#include "airfoil.h"
#include "geometrystation.h"
#include "geometrysegment.h"
#include "geometrywing.h"

class Airfoil;
class GeometryStation;
class GeometrySegment;
class GeometryWing;

class Model
{

  public:

    Model();
    ~Model();

    // alternative constructor reading the model
    // from an XML document
    Model(QDomElement parent);

    // update should be called whenever something has changed
    // that could have side-effects throughout the data structure
    void update();

    // handle global properties
    QString getName();
    void setName(QString s);
    QString getAuthor();
    void setAuthor(QString s);
    double getMass();
    void setMass(double m);
    double getRefSpan();
    void setRefSpan(double span);
    double getRefArea();
    void setRefArea(double area);
    double getRefChord();
    void setRefChord(double chord);

    double getXinfinity();
    void setXinfinity(double xinf);
    int getWakePanelNumber();
    void setWakePanelNumber(int n);
    double getAOA();
    void setAOA(double aoa);

    // manage some wings belonging to the model
    // n is the logical number running from 1 to NumberOfWings
    int numberOfWings();
    void addWing();
    void deleteWing(int n);
    // return a pointer to the wing object
    GeometryWing* getWing(int n);

    // manage airfoil database belonging to the model
    // index running from 1 to NumberOfAirfoils
    int numberOfAirfoils();
    void addAirfoil(Airfoil *af);
    Airfoil* airfoilDatabaseAt(int index);
    bool airfoilIsUsed(int index);
    void deleteAirfoil(int index);

    // generate source data (vtkPolyData)
    // for visualization using VTK
    void sourceVTK(
      vtkAppendPolyData *polyData,
      int highlightWing = 0,
      int highlightSegment = 0,
      int highlightStation = 0);

    // how to highlight parts of the model when rendering
    void setHighlighting (geometryHighlightMode mode);

    // write the model geometry into an XML (DOM) document
    void writeGeometryXML(
                QDomDocument doc,
                QDomElement parent);

  private:

    QString ModelName;
    QString ModelAuthor;

    double ModelMass;		// all-up weight of the model in g
    double ReferenceSpan;	// reference span for normalization in mm
    double ReferenceArea;	// reference planform area for normalization in dm²
    double ReferenceChord;      // reference chord for stability margin etc.

    double xInfinity;		// the extension of the trailing vortex sheet
    int WakePanelNumber;	// the number of panels to model a wake stripe
				// i.e. number of segments of the wake filaments
    double flowAOA;		// angle of attack for flow modeling

    int NumberOfWings;
    GeometryWing** wings;       // list of pointers to the wing objects

    // The airfoil indexing starts with zero (default initialization):
    // Zero indicates the lack of a defined airfoil -> strake
    // Numbers from 1 up refer to the airfoil database present in the model class.
    // Watch out to decrease the index by one when indexing the airfoil database!
    int NumberOfAirfoils;
    QList<Airfoil *> *airfoilDatabase;

    geometryHighlightMode highlighting;

    // check if all the references into the airfoil database
    // point to valid airfoils
    // faulty references are reset to "Strak" and reported
    void checkAirfoilReferences();

};

#endif
