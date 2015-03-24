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

#include <stdio.h>
#include "sourcedoubletmodel.h"
#include "airfoil.h"

#include <vtkCellArray.h>
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include "vtkTable.h"
#include "vtkPlot.h"
#include "vtkAxis.h"

// include <gsl/gsl_math.h>
// include <gsl/gsl_linalg.h>
// include <gsl/gsl_statistics.h>

// OpenBLAS
#include "cblas.h"
#include "lapacke.h"

// for parallel threads
#include <omp.h>

// one should keep (and check) an error flag
// indicating possible errors occured during mesh generation
SourceDoubletModel::SourceDoubletModel(Model *geometrymodel)
{
  GeometryWing *wing;		// current wing
  GeometryStation *st;		// station
  int afi;			// airfoil index
  int i1, i2;
  Vector A, B, C, D;
  Globals::MainTextDisplay->append(
      QString("\n******************************"));
  Globals::MainTextDisplay->append(
      QString("creating Surface Panel Model"));
  Globals::MainTextDisplay->append(
      QString("******************************"));
  model = geometrymodel;
  mesh = new QList<FlatPanel*>();
  wingref = new QList<int>();
  varType = new QList<variableSigularityType>();
  wakeref = new QList<int>();
  wakelines = new QList<Streamline*>();
  wake = new QList<WakeStripe*>();
  // we need the free stream vector to create the wake mesh
  double aoa = model->getAOA();
  vInfinity = Vector(cos(M_PI/180.0*aoa), 0.0, sin(M_PI/180.0*aoa));
  valid = TRUE;
  validSolution = FALSE;
  Trefftz_data_avail = FALSE;
  NumberOfPanels = 0;
  NumberOfWakes = 0;
  NumberOfFilaments = 0;
  // panel all wings
  int nwings = model->numberOfWings();
  for (int w=1; w<=nwings; w++)
  {
    Globals::MainTextDisplay->append(
      QString("\npaneling of wing No. %1").arg(w));
    wing = model->getWing(w);
    // we compute the difference of numbers before and after creating the wing model
    int wingnpan = -NumberOfPanels;
    int wingnwake = -NumberOfWakes;
    int wingnfil = -NumberOfFilaments;
    // the ends of a wing must have airfoils assigned - check!
    // this is necessary for the interpolation - we cannot extrapolate
    st = wing->getStation(1);
    afi = st->getAirfoilIndex();
    if (afi==0)
    {
      Globals::MainTextDisplay->append(
	QString("wing No. %1 has no airfoil assigned at the left tip - cannot panel").arg(w));
      valid = FALSE;
      NumberOfPanels = 0;
      NumberOfWakes = 0;
      NumberOfFilaments = 0;
      mesh->clear();
      wingref->clear();
      varType->clear();
      wakeref->clear();
      wakelines->clear();
      wake->clear();
      return;
    }
    else
    {
      // comput the left wing closure
      if (wing->getLeftClosure())
      {
	int nchord = wing->getSegment(1)->getChordN();
	Outline *ol = new Outline(2*nchord);
	ol->computePoints(model->airfoilDatabaseAt(afi));
	ol->mapToVectors(st->getNosePoint(), st->getEndPoint(), st->getUpVector(), st->getStretchZ());
	// now generate some panels
	i1=0;				// index of rear upper-side point
	i2=ol->numberPoints()-1;	// index of rear lower-side point
	A = ol->pointVec(i1);
	B = ol->pointVec(i1+1);
	C = ol->pointVec(i2-1);
	D = ol->pointVec(i2);
	if ((A-D).norm()>1e-6)				// avoid degenerate panels
	{
	  mesh->append(new FlatPanel(A,B,C,D));
	  wingref->append(w);
	  varType->append(VariableSource);
          // to be modified !!!
          wakeref->append(0);
	}
	else
	{
	  mesh->append(new FlatPanel(A,B,C));
	  wingref->append(w);
	  varType->append(VariableSource);
          // to be modified !!!
          wakeref->append(0);
	};
	NumberOfPanels++;
	i1+=1;
	i2-=1;
	while (i2-i1>=2)
	{
	  if (i2-i1==2) {
	    A = ol->pointVec(i1);
	    B = ol->pointVec(i1+1);
	    D = ol->pointVec(i2);
	    mesh->append(new FlatPanel(A,B,D));
	    wingref->append(w);
	    varType->append(VariableSource);
            // to be modified !!!
            wakeref->append(0);
	  } else {
	    A = ol->pointVec(i1);
	    B = ol->pointVec(i1+1);
	    C = ol->pointVec(i2-1);
	    D = ol->pointVec(i2);
	    mesh->append(new FlatPanel(A,B,C,D));
	    wingref->append(w);
	    varType->append(VariableSource);
            // to be modified !!!
            wakeref->append(0);
	  };
	  NumberOfPanels++;
	  i1+=1;
	  i2-=1;
	};
      };
    };
    st = wing->getStation(wing->numberOfStations());
    afi = st->getAirfoilIndex();
    if (afi==0)
    {
      Globals::MainTextDisplay->append(
	QString("wing No. %1 has no airfoil assigned at the right tip - cannot panel").arg(w));
      valid = FALSE;
      NumberOfPanels = 0;
      NumberOfWakes = 0;
      NumberOfFilaments = 0;
      mesh->clear();
      wingref->clear();
      varType->clear();
      wakeref->clear();
      wakelines->clear();
      wake->clear();
      return;
    }
    else
    {
      // comput the right wing closure
      if (wing->getRightClosure())
      {
	int nchord = wing->getSegment(wing->numberOfSegments())->getChordN();
	Outline *ol = new Outline(2*nchord);
	ol->computePoints(model->airfoilDatabaseAt(afi));
	ol->mapToVectors(st->getNosePoint(), st->getEndPoint(), st->getUpVector(), st->getStretchZ());
	// now generate some panels
	i1=0;			// index of rear upper-side pint
	i2=ol->numberPoints()-1;	// index of rear lower-side point
	A = ol->pointVec(i1);
	B = ol->pointVec(i1+1);
	C = ol->pointVec(i2-1);
	D = ol->pointVec(i2);
	if ((A-D).norm()>1e-6)				// avoid degenerate panels
	{
	  mesh->append(new FlatPanel(D,C,B,A));
	  wingref->append(w);
	  varType->append(VariableSource);
          // to be modified !!!
          wakeref->append(0);
        }
	else
	{
	  mesh->append(new FlatPanel(A,C,B));
	  wingref->append(w);
	  varType->append(VariableSource);
          // to be modified !!!
          wakeref->append(0);
	};
	NumberOfPanels++;
	i1+=1;
	i2-=1;
	while (i2-i1>=2)
	{
	  if (i2-i1==2) {
	    A = ol->pointVec(i1);
	    B = ol->pointVec(i1+1);
	    D = ol->pointVec(i2);
	    mesh->append(new FlatPanel(A,D,B));
	    wingref->append(w);
	    varType->append(VariableSource);
            // to be modified !!!
            wakeref->append(0);
	  } else {
	    A = ol->pointVec(i1);
	    B = ol->pointVec(i1+1);
	    C = ol->pointVec(i2-1);
	    D = ol->pointVec(i2);
	    mesh->append(new FlatPanel(A,D,C,B));
	    wingref->append(w);
	    varType->append(VariableSource);
            // to be modified !!!
            wakeref->append(0);
	  };
	  NumberOfPanels++;
	  i1+=1;
	  i2-=1;
	};
      };
    };
    // panel all segments of this wing
    int nsegm = wing->numberOfSegments();
    // generate panels segment by segment
    for (int seg=1; seg<=nsegm; seg++)
    {
      GeometrySegment *segment = wing->getSegment(seg);
      int npan = createSegmentModel(w, segment);
      if (npan<1)
      {
        Globals::MainTextDisplay->append(
          QString("error paneling segment No. %1, no panels generated").arg(seg));
	valid = FALSE;
	NumberOfPanels = 0;
	NumberOfWakes = 0;
	NumberOfFilaments = 0;
	mesh->clear();
	wingref->clear();
	varType->clear();
        wakeref->clear();
	wakelines->clear();
	wake->clear();
        return;
      };
    };
    wingnpan += NumberOfPanels;
    wingnwake += NumberOfWakes;
    wingnfil += NumberOfFilaments;
    // count non-flat panels
    int nonflat=0;
    for (int i=NumberOfPanels-wingnpan; i<NumberOfPanels; i++)
      if (!(mesh->at(i)->panelIsPlanar())) nonflat++;
    Globals::MainTextDisplay->append(QString("  %1 panels generated").arg(wingnpan));
    if (nonflat>0)
      Globals::MainTextDisplay->append(QString("  ! warning : %1 non-flat panels generated").arg(nonflat));
    Globals::MainTextDisplay->append(QString("  %1 wake filaments generated").arg(wingnfil));
    Globals::MainTextDisplay->append(QString("  %1 wake stripes generated").arg(wingnwake));
  };
  // check if numbers are correct
  if (NumberOfPanels != mesh->count())
  {
    Globals::MainTextDisplay->append(
      QString("\nThis should never happen - program error"));
    Globals::MainTextDisplay->append(
      QString("The length of the list of panels differes from the counted number of generated panels."));
    valid = FALSE;
  };
  if (NumberOfPanels != varType->count())
  {
    Globals::MainTextDisplay->append(
      QString("\nThis should never happen - program error"));
    Globals::MainTextDisplay->append(
      QString("The length of the list of panel variable types differes from the number of generated panels."));
    valid = FALSE;
  };
  if (NumberOfWakes != wake->count())
  {
    Globals::MainTextDisplay->append(
      QString("\nThis should never happen - program error"));
    Globals::MainTextDisplay->append(
      QString("The length of the list of wakes differes from the counted number of generated stripes."));
    valid = FALSE;
  };
  if (NumberOfFilaments != wakelines->count())
  {
    Globals::MainTextDisplay->append(
      QString("\nThis should never happen - program error"));
    Globals::MainTextDisplay->append(
      QString("The length of the list of wake filaments differes from the counted number of generated filaments."));
    valid = FALSE;
  };
  if (valid)
  {
    NumberCP = NumberOfPanels + NumberOfWakes;
    ControlPoint = new Vector[NumberCP];
    normal = new Vector[NumberCP];
    BC = new boundaryConditionType[NumberCP];
    // geometry setup of the boundary conditions
    // collect control points and normal directions
    int indexCP=0;
    for (int icp=0; icp<NumberOfPanels; icp++)
    {
      FlatPanel *p = mesh->at(icp);
      // we use control points displaced to the inside of the panel
      ControlPoint[indexCP] = p->panelCenter() + p->panelNormal()*0.0001;
      normal[indexCP] = p->panelNormal();
      BC[indexCP] = PerturbationPotentialBC;
      indexCP++;
    };
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      WakeStripe* wk = wake->at(iw);
      ControlPoint[indexCP] = wk->wakeCP();
      normal[indexCP] = wk->wakeNormal();
      BC[indexCP] = NormalVelocityBC;
      indexCP++;
    };
    Globals::MainTextDisplay->append(QString("\nPaneling done.\n"));
    Globals::MainTextDisplay->append(QString("in total %1 free parameters.\n").arg(NumberCP));
  } else {
    NumberOfPanels = 0;
    NumberOfWakes = 0;
    NumberOfFilaments = 0;
    mesh->clear();
    wingref->clear();
    varType->clear();
    wakeref->clear();
    wakelines->clear();
    wake->clear();
    Globals::MainTextDisplay->append(QString("\nPaneling failed.\n"));
  };
}

SourceDoubletModel::~SourceDoubletModel()
{
  for (int i=0; i<NumberOfPanels; i++)
    delete mesh->at(i);
  delete mesh;
  delete varType;
  // printf("delete panels ok.\n");
  for (int i=0; i<NumberOfWakes; i++)
    delete wake->at(i);
  delete wake;
  // printf("delete wake ok.\n");
  for (int i=0; i<NumberOfFilaments; i++)
    delete wakelines->at(i);
  delete wakelines;
  // printf("delete filaments ok.\n");
  // free all matrices and vectors belonging to the solution
  if (valid==TRUE)
  {
    delete[] ControlPoint;
    delete[] normal;
    delete[] BC;
  };
  if (validSolution==TRUE)
  {
    delete[] sigSolution;
    delete[] muSolution;
    delete[] vSolution;
    delete[] cpSolution;
    delete[] phiSolution;
    delete[] sourceInducedVelocity;
    delete[] doubletInducedVelocity;
    delete[] sourceInducedPotential;
    delete[] doubletInducedPotential;
  };
  if (Trefftz_data_avail==TRUE)
  {
    delete[] Trefftz_x1;
    delete[] Trefftz_x2;
    delete[] Trefftz_xc;
    delete[] Trefftz_gamma;
    delete[] Trefftz_vi;
    delete[] Trefftz_cl;
    delete[] Trefftz_cd;
  };
}

int SourceDoubletModel::createSegmentModel(int wingno, GeometrySegment *segment)
{
  // corner points of generated panels
  Vector A,B, C, D;
  // streamlines emanating from the trailing edge
  Streamline *leftFilament, *rightFilament;
  // leftspan and rightspan are the distances from the left/right end
  // of the current segment to the next station with a known airfoil
  double leftspan = 0.0;
  double centerspan = segment->segmentSpan();
  double rightspan = 0.0;
  GeometrySegment *test = segment;
  GeometryStation *lst = test->getLeftStation();
  int lafi = lst->getAirfoilIndex();
  while (lafi==0)
  {
    // if the left station of the segment has no airfoil assigned
    // move on to the next and add the segment span to the accumulated distance
    test = test->getLeftNeighbour();
    if (test==NULL)
    {
      Globals::MainTextDisplay->append(
        QString("\nThis should never happen - program error!"));
      Globals::MainTextDisplay->append(
        QString("didn't find any airfoil on the left side - a check should have failed earlier on"));
      return -1;
    }
    lst = test->getLeftStation();
    lafi = lst->getAirfoilIndex();
    leftspan += test->segmentSpan();
  }
  test = segment;
  GeometryStation *rst = test->getRightStation();
  int rafi = rst->getAirfoilIndex();
  while (rafi==0)
  {
    // if the right station of the segment has no airfoil assigned
    // move on to the next and add the segment span to the accumulated distance
    test = test->getRightNeighbour();
    if (test==NULL)
    {
      Globals::MainTextDisplay->append(
        QString("\nThis should never happen - program error!"));
      Globals::MainTextDisplay->append(
        QString("didn't find any airfoil on the right side - a check should have failed earlier on"));
      return -1;
    }
    rst = test->getRightStation();
    rafi = rst->getAirfoilIndex();
    rightspan += test->segmentSpan();
  }
  // compute the section outlines between which the interpolation is to be done
  int nchord = segment->getChordN();
  Outline *leftSection = new Outline(2*nchord);
  leftSection->computePoints(model->airfoilDatabaseAt(lafi));
  Outline *rightSection = new Outline(2*nchord);
  rightSection->computePoints(model->airfoilDatabaseAt(rafi));
  // generate interpolated outlines for the stations
  // normalization of the interpolation factors is done
  // by the outline interpolating constructor
  lst = segment->getLeftStation();
  Outline *leftStationOL = new Outline(
    leftSection, rightspan+centerspan, rightSection, leftspan);
  leftStationOL->mapToVectors(
    lst->getNosePoint(), lst->getEndPoint(), lst->getUpVector(), lst->getStretchZ());
  double leftStationSpan = lst->getS();
  double leftStationChord = lst->getChord();
  rst = segment->getRightStation();
  Outline *rightStationOL = new Outline(
    leftSection, rightspan, rightSection, leftspan+centerspan);
  rightStationOL->mapToVectors(
    rst->getNosePoint(), rst->getEndPoint(), rst->getUpVector(), rst->getStretchZ());
  double rightStationSpan = rst->getS();
  double rightStationChord = rst->getChord();
  // the segment is subdivided in stripes of panels
  int nstripes = segment->getSpanN();
  int npan = 0;
  for (int stripe=0; stripe<nstripes; stripe++)
  {
    // interpolate the camberlines
    // Now that the interpolation for the segment endstation airfoils has been done
    // the leftSection/rightSection can be reused for the stripe left and right sections
    delete leftSection;
    leftSection = new Outline(
      leftStationOL, nstripes-stripe, rightStationOL, stripe);
    delete rightSection;
    rightSection = new Outline(
      leftStationOL, nstripes-stripe-1, rightStationOL, stripe+1);
    // spanwise position of the stripe (wing referenced)
    double stripespan =
      ((double)(nstripes-stripe)-0.5)/(double)nstripes*leftStationSpan +
      ((double)stripe+0.5)/(double)nstripes*rightStationSpan;
    // wing chord at the position of the stripe start
    double stripechord =
      ((double)(nstripes-stripe)-0.5)/(double)nstripes*leftStationChord +
      ((double)stripe+0.5)/(double)nstripes*rightStationChord;
    // now generate some panels
    for (int i=0; i<2*nchord; i++)
    {
      A = leftSection->pointVec(i);
      B = leftSection->pointVec(i+1);
      C = rightSection->pointVec(i);
      D = rightSection->pointVec(i+1);
      mesh->append(new FlatPanel(A,C,D,B));
      wingref->append(wingno);
      varType->append(VariableDoublet);
      // to be modified !!!
      wakeref->append(0);
      NumberOfPanels++;
      npan++;
    };
    // create trailing wake
    Vector start = leftSection->pointVec(0);
    Vector direction = leftSection->wakeDirection();
    bool exists=FALSE;
    for (int iw=0; iw<NumberOfFilaments; iw++)
    {
      if ((wakelines->at(iw)->point(0)-start).norm()<1.0e-3)
      {
	exists = TRUE;
	leftFilament = wakelines->at(iw);
      }
    };
    // if the filament already exists (usually true for the left
    // filament of a stripe) there is no need to create a new one
    if (!exists)
    {
      leftFilament = new Streamline(start,1000.0*(model->getXinfinity()),model->getWakePanelNumber(),
		direction,vInfinity);
      wakelines->append(leftFilament);
      NumberOfFilaments++;
    };
    start = rightSection->pointVec(0);
    direction = leftSection->wakeDirection();
    exists=FALSE;
    for (int iw=0; iw<NumberOfFilaments; iw++)
    {
      if ((wakelines->at(iw)->point(0)-start).norm()<1.0e-3)
      {
	exists = TRUE;
	rightFilament = wakelines->at(iw);
      }
    };
    if (!exists)
    {
      rightFilament = new Streamline(start,1000.0*(model->getXinfinity()),model->getWakePanelNumber(),
		direction,vInfinity);
      wakelines->append(rightFilament);
      NumberOfFilaments++;
    };
    WakeStripe *wk = new WakeStripe(leftFilament,rightFilament);
    wk->setWingIndex(wingno);
    wk->setWingSpanPos(stripespan);
    wk->setWingChord(stripechord);
    wake->append(wk);
    NumberOfWakes++;
  };
  // clean up
  delete leftSection;
  delete rightSection;
  delete leftStationOL;
  delete rightStationOL;
  return npan;
}

bool SourceDoubletModel::isValid()
{
  return(valid);
}

bool SourceDoubletModel::isSolved()
{
  return(validSolution);
}

int SourceDoubletModel::numberWakes()
{
  return(NumberOfWakes);
}

void SourceDoubletModel::sourcePanelsVTK(vtkPolyData *polyData,
					 bool showPanelLines,
					 bool showNormals,
					 bool showFlowDirection,
					 flowRenderMode renderMode)
{
  FlatPanel *pan;	// a panel
  int NP;		// number of points for this panel
  int index;
  Vector p;
  double color;
  // restore initial state and release memory
  polyData->Initialize();
  // all point coordinates
  vtkPoints *pts = vtkPoints::New();
  // all lines that connect the given points
  vtkCellArray *lines = vtkCellArray::New();
  vtkFloatArray *scalar = vtkFloatArray::New();
  // precompute the number of points needed for rendering
  int nop=0;
  for (int i=0; i<NumberOfPanels; i++)
  {
    nop+=mesh->at(i)->panelNP();
  };
  if (showNormals) nop+=2*NumberOfPanels;
  if (showFlowDirection && validSolution) nop+=2*NumberOfPanels;
  pts->SetNumberOfPoints(nop);
  scalar->SetNumberOfValues(nop);
  // set points of panel corners
  index=0;
  for (int i=0; i<NumberOfPanels; i++)
  {
    pan=mesh->at(i);
    NP=pan->panelNP();
    for (int np=0; np<NP; np++)
    {
      p=pan->panelPoint(np);
      if (validSolution)
      {
	switch (renderMode)
	{
	  // any scaling is done in the UI
	  case RenderCP : 	color = cpSolution[i]; break;
	  case RenderSource : 	color = sigSolution[i]; break;
	  case RenderDoublet : 	color = muSolution[i]; break;
	  default : 		color = 0.0;
	};
	scalar->SetValue(index, color);
      }
      else scalar->SetValue(index, 0.0);
      pts->InsertPoint(index++, p.x, p.y, p.z);
    };
  };
  // set points for panel normals
  if (showNormals)
  {
    for (int i=0; i<NumberOfPanels; i++)
    {
      p=mesh->at(i)->panelCenter();
      scalar->SetValue(index, -1000.0);
      // index keeps running
      pts->InsertPoint(index++, p.x, p.y, p.z);
      p+=mesh->at(i)->panelNormal()*20.0;
      scalar->SetValue(index, -1000.0);
      // index keeps running
      pts->InsertPoint(index++, p.x, p.y, p.z);
    };
  };
  // set points for flow direction
  if (showFlowDirection && validSolution)
  {
    for (int i=0; i<NumberOfPanels; i++)
    {
      p=mesh->at(i)->panelCenter();
      scalar->SetValue(index, 1000.0);
      // index keeps running
      pts->InsertPoint(index++, p.x, p.y, p.z);
      p+=vSolution[i]*20.0;
      scalar->SetValue(index, 1000.0);
      // index keeps running
      pts->InsertPoint(index++, p.x, p.y, p.z);
    };
  };
  polyData->SetPoints(pts);
  pts->Delete();
  polyData->GetPointData()->SetScalars(scalar);
  scalar->Delete();

  // connect the given points with lines
  index=0;
  for (int i=0; i<NumberOfPanels; i++)
  {
    pan=mesh->at(i);
    NP=pan->panelNP();
    for (int np=0; np<NP-1; np++)
    {
      if (showPanelLines)
      {
	lines->InsertNextCell(2);
	lines->InsertCellPoint(index);
	lines->InsertCellPoint(index+1);
      };
      index++;
    };
    if (showPanelLines)
    {
      lines->InsertNextCell(2);
      lines->InsertCellPoint(index);
      lines->InsertCellPoint(index-NP+1);
    };
    index++;
  };
  // index keeps running
  if (showNormals)
  {
    for (int i=0; i<NumberOfPanels; i++)
    {
      lines->InsertNextCell(2);
      lines->InsertCellPoint(index++);
      lines->InsertCellPoint(index++);
    };
  };
  // index keeps running
  if (showFlowDirection && validSolution)
  {
    for (int i=0; i<NumberOfPanels; i++)
    {
      lines->InsertNextCell(2);
      lines->InsertCellPoint(index++);
      lines->InsertCellPoint(index++);
    };
  };
  polyData->SetLines(lines);
  lines->Delete();

  if (renderMode != WireFrame)
  {
    // polygons
    vtkCellArray *polys = vtkCellArray::New();
    index=0;
    for (int i=0; i<NumberOfPanels; i++)
    {
      pan=mesh->at(i);
      NP=pan->panelNP();
      if (NP==4)
      {
	vtkQuad *quad = vtkQuad::New();
	quad->GetPointIds()->SetId(0,index++);
	quad->GetPointIds()->SetId(1,index++);
	quad->GetPointIds()->SetId(2,index++);
	quad->GetPointIds()->SetId(3,index++);
	polys->InsertNextCell(quad);
	quad->Delete();
      };
      if (NP==3)
      {
	vtkTriangle *tri = vtkTriangle::New();
	tri->GetPointIds()->SetId(0,index++);
	tri->GetPointIds()->SetId(1,index++);
	tri->GetPointIds()->SetId(2,index++);
	polys->InsertNextCell(tri);
	tri->Delete();
      };
    };
    polyData->SetPolys(polys);
    polys->Delete();
  };

}

void SourceDoubletModel::sourceWakeVTK(vtkPolyData *polyData)
{
  int NumberOfPoints = 0;		// number of points
  for (int i=0; i<NumberOfFilaments; i++)
  {
    Streamline* wl = wakelines->at(i);
    NumberOfPoints += wl->segmentN()+1;
  };
  // first set all point coordinates
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(NumberOfPoints);
  int index=0;
  for (int i=0; i<NumberOfFilaments; i++)
  {
    Streamline* wl = wakelines->at(i);
    int np = wl->segmentN()+1;
    for (int ip=0; ip<np; ip++)
    {
      Vector p = wl->point(ip);
      pts->InsertPoint(index++, p.x, p.y, p.z);
    };
  };
  polyData->SetPoints(pts);
  pts->Delete();
  // all lines that connect the given points
  vtkCellArray *lines = vtkCellArray::New();
  index=0;
  for (int i=0; i<NumberOfFilaments; i++)
  {
    Streamline* wl = wakelines->at(i);
    int ns = wl->segmentN();
    for (int line=0; line<ns; line++)
    {
      lines->InsertNextCell(2);
      lines->InsertCellPoint(index++);
      lines->InsertCellPoint(index);
    };
    index++;
  };
  polyData->SetLines(lines);
  lines->Delete();
}

void SourceDoubletModel::exportSTL(QTextStream *ts)
{
  Vector p;				// a point
  if (valid)
  {
    for (int i=0; i<NumberOfPanels; i++)
    {
      FlatPanel *pan = mesh->at(i);	// a panel
      int NP = pan->panelNP();		// number of points for this panel
      if (NP==3)
      {					// output triangles directly
	Vector n = pan->panelNormal();	// outer normal direction
	*ts << QString("  facet normal %1 %2 %3")
	  .arg(n.x,10,'f',6).arg(n.y,10,'f',6).arg(n.z,10,'f',6);
	*ts << endl;
	*ts << QString("    outer loop") << endl;
	p = pan->panelPoint(0);
	*ts << QString("      vertex %1 %2 %3")
	  .arg(p.x,10,'f',3).arg(p.y,10,'f',3).arg(p.z,10,'f',3);
	*ts << endl;
	p = pan->panelPoint(1);
	*ts << QString("      vertex %1 %2 %3")
	  .arg(p.x,10,'f',3).arg(p.y,10,'f',3).arg(p.z,10,'f',3);
	*ts << endl;
	p = pan->panelPoint(2);
	*ts << QString("      vertex %1 %2 %3")
	  .arg(p.x,10,'f',3).arg(p.y,10,'f',3).arg(p.z,10,'f',3);
	*ts << endl;
	*ts << QString("    endloop") << endl;
	*ts << QString("  endfacet") << endl;
      }
      else
      {
	for (int k=0; k<NP; k++)		// resolve panel into NP triangles
	{
	  Vector n = pan->panelNormal();	// outer normal direction
	  Vector c = pan->panelCenter();	// panel center
	  *ts << QString("  facet normal %1 %2 %3")
	    .arg(n.x,10,'f',6).arg(n.y,10,'f',6).arg(n.z,10,'f',6);
	  *ts << endl;
	  *ts << QString("    outer loop") << endl;
	  Vector p1 = pan->panelPoint(k);	// first point of the triangle
	  Vector p2;				// second point of the triangle
	  if (k == NP-1)			// last triangle
            p2 = pan->panelPoint(0);
	  else
            p2 = pan->panelPoint(k+1);
	  *ts << QString("      vertex %1 %2 %3")
	    .arg(p1.x,10,'f',3).arg(p1.y,10,'f',3).arg(p1.z,10,'f',3);
	  *ts << endl;
	  *ts << QString("      vertex %1 %2 %3")
	    .arg(p2.x,10,'f',3).arg(p2.y,10,'f',3).arg(p2.z,10,'f',3);
	  *ts << endl;
	  // third point is always the panel center
	  *ts << QString("      vertex %1 %2 %3")
	    .arg(c.x,10,'f',3).arg(c.y,10,'f',3).arg(c.z,10,'f',3);
	  *ts << endl;
	  *ts << QString("    endloop") << endl;
	  *ts << QString("  endfacet") << endl;
	};
      };
    };
  };
}

void SourceDoubletModel::runModelAOA(double aoa)
{
  FlatPanel *p;
  WakeStripe *wk;
  Vector cp;		// control point
  Vector n;		// panel normal at the control point
  Vector vind;		// induced velocity
  double phi;		// induced potential
  double *M;		// influence matrix
  double *w;		// normal component of induced velocity
  double *perturbation;	// perturbation potential
  // check preconditions
  if (!valid) return;

  Globals::MainTextDisplay->append(
      QString("\n**********************************************************"));
  Globals::MainTextDisplay->append(
      QString("running Surface Panel Model for %1 deg angle of attack").arg(aoa));
  Globals::MainTextDisplay->append(
      QString("**********************************************************\n"));
  Globals::MainTextDisplay->repaint();
  Globals::MainStatusBar->showMessage(QString("computation started ..."));
  Globals::MainStatusBar->repaint();

  vInfinity = Vector(cos(M_PI/180.0*aoa), 0.0, sin(M_PI/180.0*aoa));

  if (validSolution)
  {
    delete[] sigSolution;
    delete[] muSolution;
    delete[] vSolution;
    delete[] cpSolution;
    delete[] phiSolution;
    delete[] sourceInducedVelocity;
    delete[] doubletInducedVelocity;
    delete[] sourceInducedPotential;
    delete[] doubletInducedPotential;
    validSolution = FALSE;
  };
  sigSolution = new double[NumberOfPanels];
  muSolution = new double[NumberCP];
  vSolution = new Vector[NumberCP];
  cpSolution = new double[NumberCP];
  // the potential is only computed for the surface control points
  phiSolution = new double[NumberOfPanels];
  // only the solid boundary panels have a source strength
  // velocities are computed for all panels
  // the perturbation potential is only computed for the solid boundary panels
  sourceInducedVelocity = new Vector[NumberOfPanels*NumberCP];
  sourceInducedPotential = new double[NumberOfPanels*NumberOfPanels];
  // all panels and wakes have a doublet strength
  doubletInducedVelocity = new Vector[NumberCP*NumberCP];
  doubletInducedPotential = new double[NumberCP*NumberOfPanels];
  validSolution = TRUE;

  // First we compute the aerodynamic influence coefficients (AIC)
  // for all perturbation singularities at all control points
  Globals::MainStatusBar->showMessage(QString("setting up influence matrix ..."));
  Globals::MainStatusBar->repaint();
  // show a progress bar
  QProgressDialog progress("computing influence matrix ...", "Abort",
			   0, NumberCP, Globals::MainTextDisplay);
  progress.setWindowModality(Qt::WindowModal);
  progress.setMinimumDuration(0);
  progress.show();

  // The inflence matrices are stored in a row-major order.
  // The rows contain the coefficients for all panels to a single control point.
  // The computation is done in a column-major order.
  // The panels index [ipan/iw] is the column number and the outer loop index.
  // This loop can be parallelized.
  // The control point index [icp] is the row number and the inner loop index.
  // The inner loop cannot be parallelized, because the coefficient computation
  // for a particular panel is not re-entrance safe.
  int counter = 0;
  bool abort = FALSE;
  #pragma omp parallel private(cp,p,wk,phi) shared(counter,abort)
  { // begin parallel domain
    #pragma omp single
    {
      printf("computing on %d cores\n",omp_get_num_procs());
      Globals::MainTextDisplay->append(
        QString("computing on %1 cores").arg(omp_get_num_procs()));
      printf("computing with %d threads\n",omp_get_num_threads());
      Globals::MainTextDisplay->append(
        QString("computing with %1 threads").arg(omp_get_num_threads()));
    };
    // loop over all panels (column-wise loop)
    #pragma omp for
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
    {
      // make sure abort is read from memory
      #pragma omp flush(abort)
      if (!abort)
      {
	// set a pointer to the panel being analyzed
	p = mesh->at(ipan);
        // debugging output
        // printf(" thread No. %d : panel %d\n",omp_get_thread_num(),ipan);
	// compute the influence at all control points
	for (int icp=0; icp<NumberCP; icp++)
	{
	  cp = ControlPoint[icp];
	  p->ComputePIC(cp);
	  sourceInducedVelocity[icp*NumberOfPanels+ipan] = p->ConstSourceVelocity();
	  doubletInducedVelocity[icp*NumberCP+ipan] = p->ConstDoubletVelocity();
	  // induced potentials are computed for the surface control points only
	  if (icp<NumberOfPanels)
	  {
	    sourceInducedPotential[icp*NumberOfPanels+ipan] = p->ConstSourcePotential();
	    phi = p->ConstDoubletPotential();
	    if (icp==ipan)
	      // The self induced potential is computed for the outside control point.
	      // At the inside of the panel it has the opposite sign.
	      doubletInducedPotential[icp*NumberCP+ipan] = -phi;
	    else
	      doubletInducedPotential[icp*NumberCP+ipan] = phi;
	  };
	};
      }
      else
      {
	ipan=NumberOfPanels;
      };
      // all threads increment the counter when they have finished one loop
      // only one thread exclusively is allowed in this section of code
      #pragma omp critical
      {
	#pragma omp flush(counter)
	counter++;
	#pragma omp flush(counter)
      }
      // only the first thread updates the progress bar
      if(omp_get_thread_num() == 0)
      {
	progress.setValue(counter);
	if (progress.wasCanceled()) abort=TRUE;
	#pragma omp flush(abort)
      };
    }; // end of panel loop
    // Now we append the influence of the wakes to the same control points.
    // At the end the matrix M has (NumberOfPanels+NumberOfWakes) columns.
    // loop over all wakes (column-wise loop continued)
    #pragma omp for
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      // make sure abort is read from memory
      #pragma omp flush(abort)
      if (!abort)
      {
	// set a pointer to the wake being analyzed
	wk = wake->at(iw);
	// compute the influence at all control points
	for (int icp=0; icp<NumberCP; icp++)
	{
	  cp = ControlPoint[icp];
	  doubletInducedVelocity[icp*NumberCP+NumberOfPanels+iw] = wk->inducedVelocity(cp);
	  // induced potentials are computed for the surface control points only
	  if (icp<NumberOfPanels)
	    doubletInducedPotential[icp*NumberCP+NumberOfPanels+iw] = wk->inducedPotential(cp);
	};
      }
      else
      {
	iw=NumberOfWakes;
      };
      // all threads increment the counter when they have finished one loop
      // only one thread exclusively is allowed in this section of code
      #pragma omp critical
      {
	#pragma omp flush(counter)
	counter++;
	#pragma omp flush(counter)
      }
      // only the first thread updates the progress bar
      if(omp_get_thread_num() == 0)
      {
	progress.setValue(counter);
	if (progress.wasCanceled()) abort=TRUE;
	#pragma omp flush(abort)
      };
    }; // end of wake loop
  } // end parallel domain

  if (abort)
  {
    Globals::MainTextDisplay->append(QString("computation aborted\n"));
    delete[] sigSolution;
    delete[] muSolution;
    delete[] vSolution;
    delete[] cpSolution;
    delete[] phiSolution;
    delete[] sourceInducedVelocity;
    delete[] doubletInducedVelocity;
    delete[] sourceInducedPotential;
    delete[] doubletInducedPotential;
    validSolution = FALSE;
    return;
  };

  Globals::MainStatusBar->showMessage(QString("setting up system of equations ..."));
  Globals::MainStatusBar->repaint();
  progress.reset();
  progress.setLabelText("setting up system of equations ...");
  progress.show();

  // The source strength is known from the beginning.
  // It makes the interior normal velocity equal to the free-stream normal velocity
  // while the normal velocity is zero at the outside of the panel
  // sigma = - vinf . n
  for (int ipan=0; ipan<NumberOfPanels; ipan++)
  {
    // if (varType->at(ipan) == VariableSource)
    //   sigSolution[ipan] = 0.0;
    // else
    sigSolution[ipan] = -dot(vInfinity,normal[ipan]);
    muSolution[ipan] = 0.0;
  }

  // The inflence matrices are computed in a row-major order.
  // The rows contain the coefficients for all panels to a single control point.
  // This way, the control point index is the row number and the outer loop index.
  // The panels index is the column number and the inner loop index.
  M = new double[NumberCP*NumberCP];
  w = new double[NumberCP];
  perturbation = new double[NumberCP];  // perturbation potential
  // compute the induction of all panels at the given control points
  for (int icp=0; icp<NumberCP; icp++)
  {
    progress.setValue(icp);
    n = normal[icp];
    // solid boundaries and wake CPs have different boundary conditions
    if (icp<NumberOfPanels)
    // for the solid boundary panels we impose the
    // zero internal perturbation potential boundary condition
    {
      perturbation[icp] = 0.0;
      for (int ipan=0; ipan<NumberOfPanels; ipan++)
      {
	// the sources add to the known perturbations
	perturbation[icp] += -sourceInducedPotential[icp*NumberOfPanels+ipan] * sigSolution[ipan];
	// the doublets have to be solved for
	M[icp*NumberCP+ipan] = doubletInducedPotential[icp*NumberCP+ipan];
      };
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	// wakes have no source terms
	// the doublets have to be solved for
	M[icp*NumberCP+NumberOfPanels+iw] = doubletInducedPotential[icp*NumberCP+NumberOfPanels+iw];
      };
    }
    else
    // for the wake control points we impose the
    // zero normal velocity boundary condition
    {
      // the perturbations contain the normal components of the free-stream velocity
      // plus the contribution of the source terms.
      perturbation[icp] = -dot(vInfinity,n);
      for (int ipan=0; ipan<NumberOfPanels; ipan++)
      {
	// the sources add to the known perturbations
	perturbation[icp] += -dot(sourceInducedVelocity[icp*NumberOfPanels+ipan],n) * sigSolution[ipan];
	// the doublets have to be solved for (by velocity BC)
	M[icp*NumberCP+ipan] = dot(doubletInducedVelocity[icp*NumberCP+ipan],n);
      };
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	// wakes have no source terms
	// the doublets have to be solved for (by velocity BC)
	M[icp*NumberCP+NumberOfPanels+iw] = dot(doubletInducedVelocity[icp*NumberCP+NumberOfPanels+iw],n);
      };
    };
  };

  // now we solve the system of linear equations
  Globals::MainStatusBar->showMessage(QString("solving system of equations ..."));
  Globals::MainStatusBar->repaint();

  // use OpenBLAS for solving the linear system of equations
  int typepar = openblas_get_parallel();
  if (typepar==0) Globals::MainTextDisplay->append(QString("OpenBLAS type 0: sequential"));
  if (typepar==1) Globals::MainTextDisplay->append(QString("OpenBLAS type 1: pthread"));
  if (typepar==2) Globals::MainTextDisplay->append(QString("OpenBLAS type 2: OpenMP"));
  int *ipiv = (int*)malloc(sizeof(int)*NumberCP);
  int info = LAPACKE_dgesv(
      LAPACK_ROW_MAJOR,		// storage ordering of the matrix
      NumberCP,			// LDA : the leading array dimension = number of rows
      1,			// number of right-hand side vectors
      M,			// the matrix
      NumberCP,			// LDB : the leading dimension of the RHS vectors
      ipiv,			// workspace for the pivot vector
      perturbation,		// right-hand side matrix
      1);			// number of right-hand side vectors
  free(ipiv);
  if (info!=0) Globals::MainTextDisplay->append(QString("OpenBLAS error: %1").arg(info));
  Globals::MainStatusBar->showMessage(QString("system of equations done"));
  // perturbation now contains the solution vector
  // copy the solution into the permanent storage
  for (int ipan=0; ipan<NumberOfPanels; ipan++)
    muSolution[ipan]=perturbation[ipan];
  for (int iw=0; iw<NumberOfWakes; iw++)
    muSolution[NumberOfPanels+iw]=perturbation[NumberOfPanels+iw];

  // compute the resulting flow velocities and perturbation potentials
  Vector *SIV = sourceInducedVelocity;
  Vector *DIV = doubletInducedVelocity;
  double *SIP = sourceInducedPotential;
  double *DIP = doubletInducedPotential;
  SIV = sourceInducedVelocity;
  DIV = doubletInducedVelocity;
  for (int icp=0; icp<NumberCP; icp++)
  {
    n = normal[icp];
    vind = vInfinity;
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
      vind += *SIV++ * sigSolution[ipan];
    for (int ipan=0; ipan<NumberCP; ipan++) // this also includes the wakes
      vind += *DIV++ * muSolution[ipan];
    vSolution[icp] = vind;
    // we also check the normal velocity component
    w[icp] = fabs(dot(vind,n));
    // compute the pressure coefficient
    cpSolution[icp] = 1-vind.sqnorm();
  };
  SIP = sourceInducedPotential;
  DIP = doubletInducedPotential;
  for (int icp=0; icp<NumberOfPanels; icp++)
  {
    phi = 0.0;
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
      phi += *SIP++ * sigSolution[ipan];
    for (int ipan=0; ipan<NumberCP; ipan++) // this also includes the wakes
      phi += *DIP++ * muSolution[ipan];
    // induced perturbation potential
    phiSolution[icp] = phi;
  };

  // check extrema of the solution vectors
  double wMax = w[0];
  muMin = muSolution[0];
  muMax = muSolution[0];
  sigMin = sigSolution[0];
  sigMax = sigSolution[0];
  cpMin = cpSolution[0];
  cpMax = cpSolution[0];
  phiMin = phiSolution[0];
  phiMax = phiSolution[0];
  for (int icp=1; icp<NumberOfPanels; icp++)
  {
    if(wMax < w[icp]) wMax = w[icp];
    if(muMin > muSolution[icp]) muMin = muSolution[icp];
    if(muMax < muSolution[icp]) muMax = muSolution[icp];
    if(sigMin > sigSolution[icp]) sigMin = sigSolution[icp];
    if(sigMax < sigSolution[icp]) sigMax = sigSolution[icp];
    if(cpMin > cpSolution[icp]) cpMin = cpSolution[icp];
    if(cpMax < cpSolution[icp]) cpMax = cpSolution[icp];
    if(phiMin > phiSolution[icp]) phiMin = phiSolution[icp];
    if(phiMax < phiSolution[icp]) phiMax = phiSolution[icp];
  };

  delete[] w;
  delete[] perturbation;
  delete[] M;

  Globals::MainStatusBar->showMessage(QString("solution ready."));
  Globals::MainStatusBar->repaint();
  Globals::MainTextDisplay->append(QString("solution ready.\n"));
  Globals::MainTextDisplay->append(
    QString("source strength ranges %1 ... %2").arg(sigMin).arg(sigMax));
  Globals::MainTextDisplay->append(
    QString("doublet strength ranges %1 ... %2").arg(muMin).arg(muMax));
  Globals::MainTextDisplay->append(
    QString("max. residual normal velocity = %1").arg(wMax));
  Globals::MainTextDisplay->append(
    QString("residual perturbation potential ranges %1 ... %2").arg(phiMin).arg(phiMax));
  Globals::MainTextDisplay->append(
    QString("pressure coefficient cp ranges %1 ... %2").arg(cpMin).arg(cpMax));
  Globals::MainTextDisplay->update();

}

void SourceDoubletModel::flowField(int np, Vector *x, Vector *v)
{
  if (validSolution)
  {
  }
  else
  {
  }
}

Vector SourceDoubletModel::flowPoint(Vector x)
{
  Vector v = Vector(0.0, 0.0, 0.0);
  if (validSolution)
  {
    v = vInfinity;
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
    {
      // set a pointer to the panel being analyzed
      FlatPanel *p = mesh->at(ipan);
      p->ComputePIC(x);
      v += p->ConstSourceVelocity() * sigSolution[ipan];
      v += p->ConstDoubletVelocity() * muSolution[ipan];
    };
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      // set a pointer to the wake being analyzed
      WakeStripe *wk = wake->at(iw);
      v += wk->inducedVelocity(x) * muSolution[NumberOfPanels+iw];
    };
  }
  return v;
}

void SourceDoubletModel::relaxWake()
{
  if (validSolution)
  {
    Globals::MainTextDisplay->append(QString("\n*** computing aligned wake position ***\n"));
  }
  if (validSolution)
  {
    Globals::MainTextDisplay->append(QString("\nparameters of modeled flow :\n"));
    QFont previous = Globals::MainTextDisplay->currentFont();
    QFont actual = QFont(previous);
    actual.setStyleHint(QFont::TypeWriter);
    actual.setFixedPitch(TRUE);
    actual.setKerning(FALSE);
    actual.setPointSize(8);
    Globals::MainTextDisplay->setCurrentFont(actual);
    Globals::MainTextDisplay->setTabStopWidth(40);
    Globals::MainTextDisplay->append(QString("       \t    doublet  \t              residual    \t    \tvelocity"));
    Globals::MainTextDisplay->append(QString("       \t    strength \t              normal"));
    Globals::MainTextDisplay->append(QString("       \t             \t              velocity"));
    Globals::MainTextDisplay->append(QString("----------------------------------------------------------------------------------------------------------------------------------------"));
    for (int icp=0; icp<NumberCP; icp++)
    {
      Vector n = normal[icp];
      Vector x = ControlPoint[icp];
      Vector v = flowPoint(x);
      Globals::MainTextDisplay->append(QString("%1 :\t  %2     \t  %3          \t( %4, %5, %6  )")
	.arg(icp,4)
	.arg(muSolution[icp],12,'e',3)
	.arg(dot(v,n),12,'e',3)
	.arg(v.x,8,'f',4).arg(v.y,8,'f',4).arg(v.z,8,'f',4));
    };
    Globals::MainTextDisplay->setCurrentFont(previous);
    Globals::MainTextDisplay->append(QString("\n"));
  }
}

void SourceDoubletModel::analyzePressure()
{
  Globals::MainTextDisplay->append(QString("\n*** analyzing the surface pressure ***\n"));
  // the lift direction is perpendicular to the free-flow direction
  Vector eLift = cross(vInfinity, Vector(0.0, 1.0, 0.0));
  Vector totalForce = Vector(0.0, 0.0, 0.0);
  // the total momentum about the origin (0,0,0) of the coordinate system
  Vector totalMomentum = Vector(0.0, 0.0, 0.0);
  int NumberOfWings = Globals::GeometryModel->numberOfWings();
  for (int wing=0; wing<NumberOfWings; wing++)
  {
    Vector wingForce = Vector(0.0, 0.0, 0.0);
    double surface = 0.0;
    for (int ip=0; ip<NumberOfPanels; ip++)
      if (wingref->at(ip) == wing+1)
      {
	FlatPanel *p = mesh->at(ip);
	// p-p(inf) = rho/2 V cp
	surface += p->panelArea();
	Vector force = p->panelNormal() * p->panelArea() * -cpSolution[ip];
	wingForce += force;
	Vector R = p->panelCenter();
	totalMomentum += cross(R,force);
      };
    double area = Globals::GeometryModel->getWing(wing+1)->wingArea();
    area *= 1e4; // convert dm² to mm²
    // the coefficient of lift is defined by
    // L = rho/2 * v_Inf^2 * c_L * area = F
    double cL = 2.0*dot(eLift,wingForce)/area;
    double cD = 2.0*dot(vInfinity,wingForce)/area;
    totalForce += wingForce;
    Globals::MainTextDisplay->append(
      QString("wing No. %1 : %2").arg(wing+1).arg(Globals::GeometryModel->getWing(wing+1)->getName()));
    Globals::MainTextDisplay->append(
      QString("surface area = %1").arg(surface / 1e4, 8,'f',2) + QString::fromUtf8(" dm²"));
    Globals::MainTextDisplay->append(
      QString("ref. area    = %1").arg(area / 1e4, 8,'f',2) + QString::fromUtf8(" dm²"));
      Globals::MainTextDisplay->append(
        QString("wing force = (%1,  %2,  %3)").arg(2e-4*wingForce.x,8,'f',4).arg(2e-4*wingForce.y,8,'f',4).arg(2e-4*wingForce.z,8,'f',4));
    Globals::MainTextDisplay->append(
      QString("cL = %1").arg(cL,10,'f',6));
    Globals::MainTextDisplay->append(
      QString("cD = %1\n").arg(cD,10,'f',6));
  };
  Globals::MainTextDisplay->append(QString("plane total :"));
  // total area in mm²
  double RefArea = Globals::GeometryModel->getRefArea() * 1.0e4;
  Globals::MainTextDisplay->append(
    QString("ref. area = %1").arg(1.0e-4*RefArea) + QString::fromUtf8(" dm²"));
  Globals::MainTextDisplay->append(
    QString("total force = (%1,  %2,  %3)").arg(2e-4*totalForce.x,8,'f',4).arg(2e-4*totalForce.y,8,'f',4).arg(2e-4*totalForce.z,8,'f',4));
  double cL = 2.0*dot(eLift,totalForce)/RefArea;
  Globals::MainTextDisplay->append(
    QString("cL = %1").arg(cL,10,'f',6));
  // double RefLength = Globals::GeometryModel->getRefChord();
  Globals::MainTextDisplay->append(
    QString("total momentum about (0,0,0) : cM = (%1,  %2,  %3) mm")
      .arg(totalMomentum.x/RefArea).arg(totalMomentum.y/RefArea).arg(totalMomentum.z/RefArea));
  double cMy = totalMomentum.y/RefArea;
  Globals::MainTextDisplay->append(
    QString("center of lift xL = %1 mm").arg(-2.0*cMy/cL));
}

void SourceDoubletModel::analyzeWake()
{
  if (validSolution)
  {
    Globals::MainTextDisplay->append(QString("\n*** analyzing the far wake ***\n"));
    if (Trefftz_data_avail)
    {
      delete[] Trefftz_x1;
      delete[] Trefftz_x2;
      delete[] Trefftz_xc;
      delete[] Trefftz_gamma;
      delete[] Trefftz_vi;
      delete[] Trefftz_cl;
      delete[] Trefftz_cd;
      Trefftz_data_avail = FALSE;
    };
    Trefftz_x1 = new Vector[NumberOfWakes];
    Trefftz_x2 = new Vector[NumberOfWakes];
    Trefftz_xc = new Vector[NumberOfWakes];
    Trefftz_gamma = new double[NumberOfWakes];
    Trefftz_vi = new Vector[NumberOfWakes];
    Trefftz_cl = new double[NumberOfWakes];
    Trefftz_cd = new double[NumberOfWakes];
    Trefftz_data_avail = TRUE;
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      WakeStripe *wk = wake->at(iw);
      // get the ends of the wake
      Vector x1 = wk->wakeInfinity(1);
      Vector x2 = wk->wakeInfinity(2);
      // project the wake points into a plane through the origin (0,0,0)
      // and perpendicular to the free-stream velocity.
      // This is analoguous to the Trefftz plane at infinity
      // if we ignore all but the wake vortices.
      Trefftz_x1[iw] = x1 - vInfinity * dot(x1,vInfinity);
      Trefftz_x2[iw] = x2 - vInfinity * dot(x2,vInfinity);
      Trefftz_xc[iw] = (Trefftz_x1[iw]+Trefftz_x2[iw])*0.5;
      // get the circulation strength of the associated horse-shoe vortex
      // gamma units are m/s mm
      Trefftz_gamma[iw] = -muSolution[NumberOfPanels+iw];
    };
    // compute total induced velocities at all center points
    // i.e. (half) the induction of the two-ends infinite vortex streamlines
    for (int icp=0; icp<NumberOfWakes; icp++)
    {
      Vector cp = Trefftz_xc[icp];
      Trefftz_vi[icp]=vInfinity;
      // sum up the induced velocities from all wake stripes
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	// left vortex filament
	Vector r = cp-Trefftz_x1[iw];
	Vector v = cross(vInfinity,r);
	Trefftz_vi[icp] -= v * Trefftz_gamma[iw] / (4.0*M_PI*r.sqnorm());
	// right vortex filament
	r = cp-Trefftz_x2[iw];
	v = cross(vInfinity,r);
	Trefftz_vi[icp] += v * Trefftz_gamma[iw] / (4.0*M_PI*r.sqnorm());
      };
    };
    double cL = 0.0;
    double cD = 0.0;
    Vector eLift = cross(vInfinity, Vector(0.0, 1.0, 0.0));
    double ref_area = 1.0e4*Globals::GeometryModel->getRefArea();
    double ref_span = Globals::GeometryModel->getRefSpan();
    for (int icp=0; icp<NumberOfWakes; icp++)
    {
      WakeStripe *wk = wake->at(icp);
      double chord = wk->getWingChord();
      // S units are mm
      Vector S = Trefftz_x2[icp] - Trefftz_x1[icp];
      // vi units are m/s
      // gamma units are m/s mm
      // force/rho units are (m/s)² mm²
      // F = rho * (vi x Gamma) * ds
      Vector force = cross(Trefftz_vi[icp], S) * Trefftz_gamma[icp];
      // F = rho/2 * v_Inf^2 * c_x * area
      // v_Inf^2 == 1
      double cl = 2.0*dot(eLift,force) / S.norm() / chord;
      Trefftz_cl[icp] = cl;
      cL += 2.0*dot(eLift,force) / ref_area;
      double cd = 2.0*dot(vInfinity,force) / S.norm() / chord;
      Trefftz_cd[icp] = cd;
      cD += 2.0*dot(vInfinity,force) / ref_area;
    };
    Globals::MainTextDisplay->append(
      QString("cL = %1").arg(cL,10,'f',6));
    Globals::MainTextDisplay->append(
      QString("cD = %1   (k=%2)\n")
	.arg(cD,10,'f',6)
	.arg(cD*M_PI*ref_span*ref_span/ref_area/cL/cL,10,'f',6) );
  };
}

void SourceDoubletModel::getSigmaRange (double *min, double *max)
{
  if (validSolution)
  {
    *min = sigMin;
    *max = sigMax;
  }
  else
  {
    *min = 0.0;
    *max = 0.0;
  };
}

void SourceDoubletModel::getMuRange (double *min, double *max)
{
  if (validSolution)
  {
    *min = muMin;
    *max = muMax;
  }
  else
  {
    *min = 0.0;
    *max = 0.0;
  };
}

void SourceDoubletModel::getCpRange (double *min, double *max)
{
  if (validSolution)
  {
    *min = cpMin;
    *max = cpMax;
  }
  else
  {
    *min = 0.0;
    *max = 0.0;
  };
}

void SourceDoubletModel::printPaneling()
{
  if (valid)
  {
    Globals::MainTextDisplay->append(QString("\npaneling of model surface :\n"));
    QFont previous = Globals::MainTextDisplay->currentFont();
    QFont actual = QFont(previous);
    actual.setStyleHint(QFont::TypeWriter);
    actual.setFixedPitch(TRUE);
    actual.setKerning(FALSE);
    actual.setPointSize(8);
    Globals::MainTextDisplay->setCurrentFont(actual);
    Globals::MainTextDisplay->setTabStopWidth(40);
    Globals::MainTextDisplay->append(QString("  \tpanel center                                   \tpanel normal"));
    Globals::MainTextDisplay->append(QString("-----------------------------------------------------------------------------------------------------"));
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
    {
      FlatPanel *p = mesh->at(ipan);
      Vector pc = p->panelCenter();
      Vector pn = p->panelNormal();
      Globals::MainTextDisplay->append(QString("  %1 :\t(%2   %3   %4  )  \t(%5   %6   %7  )")
	.arg(ipan,4)
	.arg(pc.x,8,'f',2).arg(pc.y,8,'f',2).arg(pc.z,8,'f',2)
	.arg(pn.x,8,'f',4).arg(pn.y,8,'f',4).arg(pn.z,8,'f',4));
    };
    Globals::MainTextDisplay->setCurrentFont(previous);
    Globals::MainTextDisplay->append(QString("\n"));
  }
}

void SourceDoubletModel::printSolution()
{
  if (validSolution)
  {
    Globals::MainTextDisplay->append(QString("\nparameters of modeled flow :\n"));
    QFont previous = Globals::MainTextDisplay->currentFont();
    QFont actual = QFont(previous);
    actual.setStyleHint(QFont::TypeWriter);
    actual.setFixedPitch(TRUE);
    actual.setKerning(FALSE);
    actual.setPointSize(8);
    Globals::MainTextDisplay->setCurrentFont(actual);
    Globals::MainTextDisplay->setTabStopWidth(40);
    Globals::MainTextDisplay->append(QString("  \t  surface   \t  source      \t  doublet     \t  residual    \t  residual    \t    surface"));
    Globals::MainTextDisplay->append(QString("  \t  pressure  \t  strength    \t  strength    \t  normal      \t  perturbation\t    velocity"));
    Globals::MainTextDisplay->append(QString("  \t  (cp)      \t              \t              \t  velocity    \t  potential   \t    "));
    Globals::MainTextDisplay->append(QString("--------------------------------------------------------------------------------------------------------------------------------"));
    for (int icp=0; icp<NumberCP; icp++)
    {
      Vector n = normal[icp];
      Vector v = vSolution[icp];
      if (icp<NumberOfPanels)
	Globals::MainTextDisplay->append(QString("%1 :\t%2\t%3\t%4\t%5\t%6\t( %7 %8 %9  )")
	  .arg(icp,4)
	  .arg(cpSolution[icp],10,'f',5)
	  .arg(sigSolution[icp],12,'e',3)
	  .arg(muSolution[icp],12,'e',3)
	  .arg(dot(v,n),12,'e',3)
	  .arg(phiSolution[icp],12,'e',3)
	  .arg(v.x,8,'f',4).arg(v.y,8,'f',4).arg(v.z,8,'f',4));
      else
	Globals::MainTextDisplay->append(QString("%1 :\t%2\t          \t%4\t%5\t        \t( %7 %8 %9  )")
	  .arg(icp,4)
	  .arg(cpSolution[icp],10,'f',5)
	  .arg(muSolution[icp],12,'e',3)
	  .arg(dot(v,n),12,'e',3)
	  .arg(v.x,8,'f',4).arg(v.y,8,'f',4).arg(v.z,8,'f',4));
    };
    Globals::MainTextDisplay->setCurrentFont(previous);
    Globals::MainTextDisplay->append(QString("\n"));
  }
}

void SourceDoubletModel::printCirculation()
{
  if (validSolution && Trefftz_data_avail)
  {
    Globals::MainTextDisplay->append(QString("\nflow properties of far wake :\n"));
    QFont previous = Globals::MainTextDisplay->currentFont();
    QFont actual = QFont(previous);
    actual.setStyleHint(QFont::TypeWriter);
    actual.setFixedPitch(TRUE);
    actual.setKerning(FALSE);
    actual.setPointSize(8);
    Globals::MainTextDisplay->setCurrentFont(actual);
    Globals::MainTextDisplay->setTabStopWidth(40);
    Globals::MainTextDisplay->append(QString("          wingref\t  span pos.   \t     gamma             \t     cl           cd                             velocity"));
    Globals::MainTextDisplay->append(QString("----------------------------------------------------------------------------------------------------------------------------------------------------------------"));
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      WakeStripe *wk = wake->at(iw);
      Globals::MainTextDisplay->append(QString("  %1 :  %2\t%3\t  %4  \t%5  %6          (%7  %8  %9 )")
	.arg(iw,4)
	.arg(wk->getWingIndex(),10)
	.arg(wk->getWingSpanPos(),12)
	.arg(muSolution[NumberOfPanels+iw],12,'e',3)
	.arg(Trefftz_cl[iw],8,'f',4)
	.arg(Trefftz_cd[iw],10,'f',6)
	.arg(Trefftz_vi[iw].x,8,'f',4).arg(Trefftz_vi[iw].y,8,'f',4).arg(Trefftz_vi[iw].z,8,'f',4)  );
    };
    Globals::MainTextDisplay->setCurrentFont(previous);
    Globals::MainTextDisplay->append(QString("\n"));
  }
}

WakeStripe* SourceDoubletModel::getWake(int iw)
{
  if (valid && (iw>=0) && (iw<NumberOfWakes))
    return(wake->at(iw));
  else
    return(0);
}

void SourceDoubletModel::sourceGammaPlot(vtkChartXY *chart)
{
  char label[40];
  // remove all previous plots
  // chart->ClearPlots();
  if (validSolution && Trefftz_data_avail)
  {
    // keep track of the extreme circulation values
    double gMin = 0.0;
    double gMax = 0.0;
    // loop through the wings
    int NumberOfWings = Globals::GeometryModel->numberOfWings();
    for (int wing=0; wing<NumberOfWings; wing++)
    {
      // count the number of wakes belonging to this wing
      int nw = 0;
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	WakeStripe *wk = wake->at(iw);
	if (wk->getWingIndex() == wing+1) nw++;
      };
      // create a table with two columns for the model data
      vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
      vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
      arrS->SetName("span");
      table->AddColumn(arrS);
      vtkSmartPointer<vtkFloatArray> arr1 = vtkSmartPointer<vtkFloatArray>::New();
      sprintf(label,"Wing No.%3d (SPM)",wing+1);
      arr1->SetName(label);
      table->AddColumn(arr1);
      // enter the postion and circulation values into the table
      table->SetNumberOfRows(nw);
      int it = 0; // index in the table
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	WakeStripe *wk = wake->at(iw);
	if (wk->getWingIndex() == wing+1)
	{
	  double s = wk->getWingSpanPos();
	  table->SetValue(it, 0, s);
	  double g = -muSolution[NumberOfPanels+iw];
	  if (g<gMin) gMin=g;
	  if (g>gMax) gMax=g;
	  table->SetValue(it, 1, g);
	  it++;
	};
      };
      // add a plot for this wing
      vtkPlot *pl = chart->AddPlot(vtkChart::POINTS);
      pl->SetInputData(table, 0, 1);
    };

    // create a table for the reference ellipse
    vtkSmartPointer<vtkTable> elltab = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> ellS = vtkSmartPointer<vtkFloatArray>::New();
    ellS->SetName("span");
    elltab->AddColumn(ellS);
    vtkSmartPointer<vtkFloatArray> ellG = vtkSmartPointer<vtkFloatArray>::New();
    ellG->SetName("elliptical");
    elltab->AddColumn(ellG);
    double ellMax=gMax;
    if (abs(gMin)>ellMax) ellMax=gMin;
    // we use 101 points to draw an ellipse with the reference span
    elltab->SetNumberOfRows(101);
    double refspan = Globals::GeometryModel->getRefSpan();
    for (int i=0; i<101; i++)
    {
      double s = (i-50.0)/100.0 * refspan;
      double g = sqrt(1.0-(2.0*s/refspan)*(2.0*s/refspan))*ellMax;
      elltab->SetValue(i, 0, s);
      elltab->SetValue(i, 1, g);
    };
    // add a plot of the reference elliptical circulation
    vtkPlot *ellpl = chart->AddPlot(vtkChart::LINE);
    ellpl->SetInputData(elltab, 0, 1);
    ellpl->SetColor(0, 0, 128, 255);

    chart->SetDrawAxesAtOrigin(false);
    chart->SetShowLegend(true);
    chart->GetAxis(vtkAxis::LEFT)->SetMinimum(gMin);
    chart->GetAxis(vtkAxis::LEFT)->SetMaximum(gMax);
    chart->GetAxis(vtkAxis::LEFT)->AutoScale();
    // chart->GetAxis(vtkAxis::LEFT)->SetNotation(3);
    // chart->GetAxis(vtkAxis::LEFT)->SetPrecision(3);
    // chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
    chart->GetAxis(vtkAxis::BOTTOM)->AutoScale();
    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("span position [mm]");
  };
}

void SourceDoubletModel::sourceClCdPlot(vtkChartXY *chart)
{
  char label[40];
  // remove all previous plots
  // chart->ClearPlots();
  // keep track of the extreme values
  if (validSolution && Trefftz_data_avail)
  {
    double clMin = 0.0;
    double clMax = 0.0;
    double cdMin = 0.0;
    double cdMax = 0.0;
    // the lift direction is perpendicular to the free-flow direction
    // Vector eLift = cross(vInfinity, Vector(0.0, 1.0, 0.0));
    // loop through the wings
    int NumberOfWings = Globals::GeometryModel->numberOfWings();
    for (int wing=0; wing<NumberOfWings; wing++)
    {
      // count the number of wakes belonging to this wing
      int nw = 0;
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	WakeStripe *wk = wake->at(iw);
	if (wk->getWingIndex() == wing+1) nw++;
      };
      // create a table with 3 columns for the model data
      vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
      vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
      arrS->SetName("span");
      table->AddColumn(arrS);
      vtkSmartPointer<vtkFloatArray> arr1 = vtkSmartPointer<vtkFloatArray>::New();
      sprintf(label,"Wing No.%3d cl (SPM)",wing+1);
      arr1->SetName(label);
      table->AddColumn(arr1);
      vtkSmartPointer<vtkFloatArray> arr2 = vtkSmartPointer<vtkFloatArray>::New();
      sprintf(label,"Wing No.%3d cd (SPM)",wing+1);
      arr2->SetName(label);
      table->AddColumn(arr2);
      // enter the postion and circulation values into the table
      table->SetNumberOfRows(nw);
      int it = 0; // index in the table
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	WakeStripe *wk = wake->at(iw);
	if (wk->getWingIndex() == wing+1)
	{
	  double s = wk->getWingSpanPos();
	  table->SetValue(it, 0, s);
	  double cl = Trefftz_cl[iw];
	  if (cl<clMin) clMin=cl;
	  if (cl>clMax) clMax=cl;
	  table->SetValue(it, 1, cl);
	  double cd = Trefftz_cd[iw];
	  if (cd<cdMin) cdMin=cd;
	  if (cd>cdMax) cdMax=cd;
	  table->SetValue(it, 2, cd);
	  it++;
	};
      };
      // add plots for this wing
      vtkPlot *pl1 = chart->AddPlot(vtkChart::LINE);
      pl1->SetInputData(table, 0, 1);
      // use bottom and left axis
      chart->SetPlotCorner( pl1, 0);
      vtkPlot *pl2 = chart->AddPlot(vtkChart::LINE);
      pl2->SetInputData(table, 0, 2);
      // use bottom and right axis
      chart->SetPlotCorner( pl2, 1);
    };
    chart->SetDrawAxesAtOrigin(false);
    chart->SetShowLegend(true);
    chart->SetBorders(20,10,20,10);
    chart->GetAxis(vtkAxis::LEFT)->SetMinimum(clMin);
    chart->GetAxis(vtkAxis::LEFT)->SetMaximum(clMax);
    chart->GetAxis(vtkAxis::LEFT)->AutoScale();
    chart->GetAxis(vtkAxis::LEFT)->SetGridVisible(true);
    chart->GetAxis(vtkAxis::LEFT)->SetTitle("lift coefficient");
    chart->GetAxis(vtkAxis::RIGHT)->SetMinimum(cdMin);
    chart->GetAxis(vtkAxis::RIGHT)->SetMaximum(cdMax);
    chart->GetAxis(vtkAxis::RIGHT)->SetGridVisible(false);
    chart->GetAxis(vtkAxis::RIGHT)->AutoScale();
    chart->GetAxis(vtkAxis::RIGHT)->SetTitle("drag coefficient");
    chart->GetAxis(vtkAxis::BOTTOM)->AutoScale();
    chart->GetAxis(vtkAxis::BOTTOM)->SetGridVisible(true);
    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("span position [mm]");
  };
}