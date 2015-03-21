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

#ifndef SOURCEDOUBLETMODEL_H
#define SOURCEDOUBLETMODEL_H

#include <QtGui>

// these defines are necessary for VTK only when building with qmake
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include <vtkPolyData.h>
#include "vtkChartXY.h"

#include <math.h>
#include "global.h"
#include "vector.h"
#include "flatpanel.h"
#include "geometrymodel.h"
#include "geometrystation.h"
#include "geometrysegment.h"
#include "geometrywing.h"
#include "streamline.h"
#include "wakestripe.h"

/*

This class implements an aerodynamic model which entirely
consists of independant (albeit neighboughring) panels.
The flow around the body is modeled by constant-strength doublet
distributions on the panels (essentially vortex rings).
Control points at which a zero normal flow condition is
applied are placed at the panel centers.

Wakes are simulated as strips of doublet panels
with constant strength equal to the doublet strength
at the trailng edge of the wing.
The wake will adapt to the flow direction
but at present is is implemented as a single panel
reaching to infinity.

*/

class SourceDoubletModel
{

  public:

    SourceDoubletModel(Model *geometrymodel);
    ~SourceDoubletModel();

    // return whether the model is created properly
    bool isValid();
    // return whether the solution for the singularity strengths has bee found
    bool isSolved();

    // generate source data (vtkPolyData)
    // for visualization using VTK
    void sourcePanelsVTK(vtkPolyData *polyData,
			 bool showPanelLines=TRUE,
			 bool showNormals=FALSE,
			 bool showFlowDirection=FALSE,
			 flowRenderMode renderMode=RenderCP
			);
    void sourceWakeVTK(vtkPolyData *polyData);

    // generate export data of all panels in STL format
    void exportSTL(QTextStream *ts);

    // setup calculation for given angle of attack
    // compute doublet strengths for flow model
    void runModelAOA(double aoa);

    // return the total flow velocity at a number of points
    // computed from the solved flow model
    void flowField(int np, Vector *x, Vector *v);
    Vector flowPoint(Vector x);
    
    // compute relaxed wake position aligned
    // with the computed flow field
    void relaxWake();
    
    // compute lift and drag from the local velocities and singularity strengths
    void analyzeCirculation();

    // compute aerodynamic coefficients by surface pressure integration
    void analyzePressure();

    // compute lift and drag from the far wake (Trefftz plane)
    void analyzeWake();

    // report methods for solution parameters
    void getSigmaRange (double *min, double *max);
    void getMuRange (double *min, double *max);
    void getCpRange (double *min, double *max);

    // generate text output
    void printPaneling();
    void printSolution();
    void printCirculation();

    // generate graphics output
    void sourceGammaPlot(vtkChartXY *chart);
    void sourceClCdPlot(vtkChartXY *chart);

private:

    Model *model;
    bool valid;
    bool validSolution;

    int NumberOfPanels;
    QList<FlatPanel*> *mesh;
    // We keep a list denoting to which wing a panel belongs
    QList<int> *wingref;
    // The model containes a number of variable singularity strengths.
    // We list which type of singularity is solved for on each panel.
    QList<variableSigularityType> *varType;

    // attached wake simulation
    // From the corner points of the trailing panels free vortex lines emanate.
    // Two neighboughring vortices define a wake stripe alike a horse-shoe vortex
    int NumberOfFilaments;
    QList<Streamline*> *wakelines;
    int NumberOfWakes;
    QList<WakeStripe*> *wake;
    // every wake strip has two parent panels,
    // one with the same orientation, one with the opposite
    QList<int> *wakeParentPositive;
    QList<int> *wakeParentNegative;

    // Flow boundary conditions are enforced at a number of control points.
    // This model has one control point per surface panel located at the panel center.
    // The panel control points are displaced 0.1mm to the in/outside
    // to stay well clear off the singularity jump across the panel
    // all the following arrays only exist if valid == TRUE
    Vector *InnerControlPoint;	// coordinates of the control points
    Vector *normal;		// surface normals at all CP
    boundaryConditionType *BC;	// the type of boundary condition to be enforced at a given CP

    // all the following arrays only exist if validSolution == TRUE
    Vector vInfinity;		// undisturbed velocity vector (free stream)
    double *sigSolution;	// the source strength for all panels
    double sigMin, sigMax;
    double *muSolution;		// the doublet strength for all panels
    double muMin, muMax;
    double *muWake;		// the doublet strength for all wakes
    Vector *vSolution;		// flow velocity at the control points
    double vMin, vMax;
    double *cpSolution;		// the pressure coefficient
    double cpMin, cpMax;
    double *phiSolution;	// perturbation potential at the inside of the panels
    double phiMin, phiMax;

    // these arrays store the influence matrices of the panels and wakes
    // for solid panels there are sources and doublets
    // for wakes there are the doublets, only
    Vector *sourceInducedVelocity;
    Vector *doubletInducedVelocity;
    Vector *wakeInducedVelocity;
    double *sourceInducedPotential;
    double *doubletInducedPotential;
    double *wakeInducedPotential;

    // wake properties in the Trefftz plane
    bool Trefftz_data_avail;	// whether these memory fields are allocated
    Vector *Trefftz_x1;		// left and right intersection points
    Vector *Trefftz_x2;		// of the vortex through the plane
    Vector *Trefftz_xc;		// center (=control) point
    double *Trefftz_gamma;	// vortex strength projected into the plane
    Vector *Trefftz_vi;		// induced velocity
    double *Trefftz_cl;
    double *Trefftz_cd;

    // Add a number of panels to the list
    // describing the (cambered) median surface of the segment.
    // The number of added panels is returned.
    // Zero or negative panel number return values indicate a triangulation error.
    int createSegmentModel(int wingref, GeometrySegment *segment);

};

#endif