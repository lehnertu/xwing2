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
#include "vortexlatticemodel.h"
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

// OpenBLAS
// include "cblas.h"
#include "lapacke.h"

// for parallel threads
#include <omp.h>

// one should keep (and check) an error flag
// indicating possible errors occured during mesh generation
VortexLatticeModel::VortexLatticeModel(Model *geometrymodel)
{
  GeometryWing *wing;		// current wing
  GeometryStation *st;		// station
  int afi;			// airfoil index
  Vector A, B, C, D;
  Globals::MainTextDisplay->append(
      QString("\n*******************************"));
  Globals::MainTextDisplay->append(
      QString("creating Vortex Lattice Model"));
  Globals::MainTextDisplay->append(
      QString("*******************************"));
  model = geometrymodel;
  mesh = new QList<FlatPanel*>();
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
      wakelines->clear();
      wake->clear();
      return;
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
      wakelines->clear();
      wake->clear();
      return;
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
    // geometry setup of the boundary conditions
    // collect control points and normal directions
    int indexCP=0;
    for (int icp=0; icp<NumberOfPanels; icp++)
    {
      FlatPanel *p = mesh->at(icp);
      ControlPoint[indexCP] = p->panelCenter();
      normal[indexCP] = p->panelNormal();
      indexCP++;
    };
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      WakeStripe* wk = wake->at(iw);
      ControlPoint[indexCP] = wk->wakeCP();
      normal[indexCP] = wk->wakeNormal();
      indexCP++;
    };
    Globals::MainTextDisplay->append(QString("\nPaneling done.\n"));
    Globals::MainTextDisplay->append(QString("in total %1 free parameters.\n").arg(NumberCP));
  } else {
    NumberOfPanels = 0;
    NumberOfWakes = 0;
    NumberOfFilaments = 0;
    mesh->clear();
    wakelines->clear();
    wake->clear();
    Globals::MainTextDisplay->append(QString("\nPaneling failed.\n"));
  };
}

VortexLatticeModel::~VortexLatticeModel()
{
  for (int i=0; i<NumberOfPanels; i++)
    delete mesh->at(i);
  delete mesh;
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
  };
  if (validSolution==TRUE)
  {
    delete[] muSolution;
    delete[] vSolution;
    delete[] doubletInducedVelocity;
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

int VortexLatticeModel::createSegmentModel(int wingref, GeometrySegment *segment)
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
  // compute the section camber lines between which the interpolation is to be done
  int nchord = segment->getChordN();
  Camberline *leftCL = new Camberline(nchord);
  leftCL->computePoints(model->airfoilDatabaseAt(lafi));
  Camberline *rightCL = new Camberline(nchord);
  rightCL->computePoints(model->airfoilDatabaseAt(rafi));
  // generate interpolated camber lines for the stations
  // normalization of the interpolation factors is done
  // by the camberline interpolating constructor
  lst = segment->getLeftStation();
  Camberline *leftStationCL = new Camberline(
    leftCL, rightspan+centerspan, rightCL, leftspan);
  leftStationCL->mapToVectors(
    lst->getNosePoint(), lst->getEndPoint(), lst->getUpVector(), lst->getStretchZ());
  double leftStationSpan = lst->getS();
  double leftStationChord = lst->getChord();
  rst = segment->getRightStation();
  Camberline *rightStationCL = new Camberline(
    leftCL, rightspan, rightCL, leftspan+centerspan);
  rightStationCL->mapToVectors(
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
    // the leftCL/rightCL can be reused for the stripe left and right sections
    delete leftCL;
    leftCL = new Camberline(
      leftStationCL, nstripes-stripe, rightStationCL, stripe);
    delete rightCL;
    rightCL = new Camberline(
      leftStationCL, nstripes-stripe-1, rightStationCL, stripe+1);
    // spanwise position of the stripe (wing referenced)
    double stripespan =
      ((double)(nstripes-stripe)-0.5)/(double)nstripes*leftStationSpan +
      ((double)stripe+0.5)/(double)nstripes*rightStationSpan;
    // wing chord at the position of the stripe start
    double stripechord =
      ((double)(nstripes-stripe)-0.5)/(double)nstripes*leftStationChord +
      ((double)stripe+0.5)/(double)nstripes*rightStationChord;
    // now generate some panels frome nose to trailing edge
    for (int i=0; i<nchord; i++)
    {
      A = leftCL->pointVec(i);
      B = leftCL->pointVec(i+1);
      C = rightCL->pointVec(i);
      D = rightCL->pointVec(i+1);
      mesh->append(new FlatPanel(A,B,D,C));
      NumberOfPanels++;
      npan++;
    };
    // create trailing wake
    Vector start = leftCL->pointVec(nchord);
    Vector direction = leftCL->wakeDirection();
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
    start = rightCL->pointVec(nchord);
    direction = leftCL->wakeDirection();
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
    wk->setWingIndex(wingref);
    wk->setWingSpanPos(stripespan);
    wk->setWingChord(stripechord);
    wake->append(wk);
    NumberOfWakes++;
  };
  // clean up
  delete leftCL;
  delete rightCL;
  delete leftStationCL;
  delete rightStationCL;
  return npan;
}

bool VortexLatticeModel::isValid()
{
  return(valid);
}

bool VortexLatticeModel::isSolved()
{
  return(validSolution);
}

void VortexLatticeModel::sourcePanelsVTK(vtkPolyData *polyData,
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

void VortexLatticeModel::sourceWakeVTK(vtkPolyData *polyData)
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

void VortexLatticeModel::runModelAOA(double aoa)
{
  FlatPanel *p;
  WakeStripe *wk;
  Vector cp;		// control point
  Vector n;		// panel normal at the control point
  Vector vind;		// induced velocity
  double *M;		// influence matrix
  double *w;		// normal component of induced velocity
  double *perturbation;	// perturbation potential
  // check preconditions
  if (!valid) return;

  Globals::MainTextDisplay->append(
      QString("\n**********************************************************"));
  Globals::MainTextDisplay->append(
      QString("running Vortex Lattice Model for %1 deg angle of attack").arg(aoa));
  Globals::MainTextDisplay->append(
      QString("**********************************************************\n"));
  Globals::MainTextDisplay->repaint();
  Globals::MainStatusBar->showMessage(QString("computation started ..."));
  Globals::MainStatusBar->repaint();

  vInfinity = Vector(cos(M_PI/180.0*aoa), 0.0, sin(M_PI/180.0*aoa));

  if (validSolution)
  {
    delete[] muSolution;
    delete[] vSolution;
    delete[] doubletInducedVelocity;
    validSolution = FALSE;
  };
  muSolution = new double[NumberCP];
  vSolution = new Vector[NumberCP];
  // velocities are computed for all panels
  // all panels and wakes have a doublet strength
  doubletInducedVelocity = new Vector[NumberCP*NumberCP];
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
  // The inner loop cannot be parrallelized, because the coefficient computation
  // for a particular panel is not re-entrance safe.
  int counter = 0;
  bool abort = FALSE;
  #pragma omp parallel private(cp,p,wk) shared(counter,abort)
  { // begin parallel domain
    #pragma omp single
    {
      Globals::MainTextDisplay->append(
        QString("computing on %1 cores").arg(omp_get_num_procs()));
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
	// compute the influence at all control points
	for (int icp=0; icp<NumberCP; icp++)
	{
	  cp = ControlPoint[icp];
	  p->ComputePIC(cp);
	  doubletInducedVelocity[icp*NumberCP+ipan] = p->ConstDoubletVelocity();
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
    delete[] muSolution;
    delete[] vSolution;
    delete[] doubletInducedVelocity;
    validSolution = FALSE;
    return;
  };

  Globals::MainStatusBar->showMessage(QString("setting up system of equations ..."));
  Globals::MainStatusBar->repaint();
  progress.reset();
  progress.setLabelText("setting up system of equations ...");
  progress.show();

  for (int ipan=0; ipan<NumberOfPanels; ipan++)
  {
    muSolution[ipan] = 0.0;
  }

  // The inflence matrices are computed in a row-major order.
  // The rows contain the coefficients for all panels to a single control point.
  // Here we assemble a linear system of equations M.w=perturbation , with
  // each row containing the boundary condition at one particular control point.
  M = new double[NumberCP*NumberCP];
  w = new double[NumberCP];
  perturbation = new double[NumberCP];  // perturbation potential
  // compute the induction of all panels at the given control points
  // icp is the index of the control point
  for (int icp=0; icp<NumberCP; icp++)
  {
    progress.setValue(icp);
    n = normal[icp];
    // we impose a zero normal velocity boundary condition
    // at all control points, body and wake included
    // the perturbations contain the normal components of the free-stream velocity
    perturbation[icp] = -dot(vInfinity,n);
    // we assemble the influence matrix
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
    {
      // the doublets have to be solved for (by velocity BC)
      M[icp*NumberCP+ipan] = dot(doubletInducedVelocity[icp*NumberCP+ipan],n);
    };
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      // the doublets have to be solved for (by velocity BC)
      M[icp*NumberCP+NumberOfPanels+iw] = dot(doubletInducedVelocity[icp*NumberCP+NumberOfPanels+iw],n);
    };
  };

  // now we solve the system of linear equations
  Globals::MainStatusBar->showMessage(QString("solving system of equations ..."));
  Globals::MainStatusBar->repaint();

  // use OpenBLAS for solving the linear system of equations
  /*
  int typepar = openblas_get_parallel();
  if (typepar==0) Globals::MainTextDisplay->append(QString("OpenBLAS type 0: sequential"));
  if (typepar==1) Globals::MainTextDisplay->append(QString("OpenBLAS type 1: pthread"));
  if (typepar==2) Globals::MainTextDisplay->append(QString("OpenBLAS type 2: OpenMP"));
  */
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
  Globals::MainStatusBar->showMessage(QString("system of equations done."));
  // perturbation now contains the solution vector
  // copy the solution into the permanent storage
  for (int ipan=0; ipan<NumberOfPanels; ipan++)
    muSolution[ipan]=perturbation[ipan];
  for (int iw=0; iw<NumberOfWakes; iw++)
    muSolution[NumberOfPanels+iw]=perturbation[NumberOfPanels+iw];

  // compute the resulting flow velocities
  Vector *DIV = doubletInducedVelocity;
  for (int icp=0; icp<NumberCP; icp++)
  {
    n = normal[icp];
    vind = vInfinity;
    for (int ipan=0; ipan<NumberCP; ipan++) // this also includes the wakes
      vind += *DIV++ * muSolution[ipan];
    vSolution[icp] = vind;
    // we also check the normal velocity component
    w[icp] = fabs(dot(vind,n));
  };

  // check extrema of the solution vectors
  double wMax = w[0];
  muMin = muSolution[0];
  muMax = muSolution[0];
  for (int icp=1; icp<NumberOfPanels; icp++)
  {
    if(wMax < w[icp]) wMax = w[icp];
    if(muMin > muSolution[icp]) muMin = muSolution[icp];
    if(muMax < muSolution[icp]) muMax = muSolution[icp];
  };

  delete[] w;
  delete[] perturbation;
  delete[] M;

  Globals::MainStatusBar->showMessage(QString("solution ready."));
  Globals::MainStatusBar->repaint();
  Globals::MainTextDisplay->append(QString("solution ready.\n"));
  Globals::MainTextDisplay->append(
    QString("doublet strength ranges %1 ... %2").arg(muMin).arg(muMax));
  Globals::MainTextDisplay->append(
    QString("max. residual normal velocity = %1").arg(wMax));
  Globals::MainTextDisplay->update();

}

void VortexLatticeModel::flowField(int np, Vector *x, Vector *v)
{
  if (validSolution)
  {
    for (int ip=0; ip<np; ip++)
      v[ip] = vInfinity;
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
    {
      // set a pointer to the panel being analyzed
      FlatPanel *p = mesh->at(ipan);
      for (int ip=0; ip<np; ip++)
      {
	p->ComputePIC(x[ip]);
	v[ip] += p->ConstDoubletVelocity() * muSolution[ipan];
      };
    };
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      // set a pointer to the wake being analyzed
      WakeStripe *wk = wake->at(iw);
      for (int ip=0; ip<np; ip++)
	v[ip] += wk->inducedVelocity(x[ip]) * muSolution[NumberOfPanels+iw];
    };
  }
  else
  {
    for (int ip=0; ip<np; ip++)
      v[ip] = Vector(0.0, 0.0, 0.0);
  }
}

void VortexLatticeModel::relaxWake()
{
  if (validSolution)
  {
    // !!! there exist several wakes - one per wing

    Globals::MainTextDisplay->append(QString("\n*** computing aligned wake position ***\n"));
    // The length of the wake as defined in the model.
    // this is the number of wake panels used in a WakeStripe.
    int NWakePanels=model->getWakePanelNumber();
    // The number of points which are used for the definition of the wake sheet
    // (i.e. the definition points of the filaments)
    int NFilPt=NumberOfFilaments*(NWakePanels+1);
    // The number of points at which the flow velocity is computed.
    // This includes all wake control points, all wake panel center
    // and some additional points outside the limiting wake filaments

    /*
    // the position of the wake filament
    Vector *x = new Vector[NP];
    Vector *ptx = x;            // a running pointer into the array of points
    for (int ifil=0; ifil<NumberOfFilaments; ifil++)
    {
      Streamline *filament=wakelines->at(ifil);
      for (int isp=0; isp<NPW; isp++)
        *ptx++ = filament->point(isp);
    };
    // compute the total flow velocity at those points
    Vector *v = new Vector[NP];
    Vector *ptv = v;            // a running pointer into the array of velocities
    flowField(NP, x, v);
    // for every filament we have to subtract the velocity contribution
    // induced by this particular filament

    // ... TBD ...
    // das geht so nicht. Man sollte sich auf das Geschwindigkeitsfeld,
    // das an den Kontrollpunkten berechnet wurde beschränken

    // now we have to compute the relaxed wake position
    Vector *relaxed = new Vector[NP];
    Vector *ptr = relaxed;            // a running pointer into the array of velocities
    for (int ifil=0; ifil<NumberOfFilaments; ifil++)
    {
      Streamline *filament = wakelines->at(ifil);
      // set the pointer to the velocities of this streamline
      ptv = v + ifil*NPW;
      // we do not change the direction at the start point
      // because this had been defined to have zero normal velocity
      // with respect to the trailing edge normal
      *ptv = filament->dir(0);
      // align the stram line to the given set of velocities
      filament->align(ptv);
      for (int isp=0; isp<NPW; isp++)
        *ptr++ = filament->point(isp);
    };
    Globals::MainTextDisplay->append(QString("\nflow velocities in the wake :\n"));
    QFont previous = Globals::MainTextDisplay->currentFont();
    QFont actual = QFont(previous);
    actual.setStyleHint(QFont::TypeWriter);
    actual.setFixedPitch(TRUE);
    actual.setKerning(FALSE);
    actual.setPointSize(8);
    Globals::MainTextDisplay->setCurrentFont(actual);
    Globals::MainTextDisplay->setTabStopWidth(80);
    Globals::MainTextDisplay->append(QString("  point coordinates : velocity - vInf = ( %1  %2  %3 )")
      .arg(vInfinity.x,8,'f',4).arg(vInfinity.y,8,'f',4).arg(vInfinity.z,8,'f',4));
    Globals::MainTextDisplay->append(QString("----------------------------------------------------------------------------------------------------------------------------------------"));
    for (int icp=0; icp<NP; icp++)
    {
      Globals::MainTextDisplay->append(QString("  ( %1  %2  %3 )  \t:  ( %4  %5  %6 )  \t:  ( %7  %8  %9 )")
        .arg(x[icp].x,9,'f',2).arg(x[icp].y,9,'f',2).arg(x[icp].z,9,'f',2)
	.arg(v[icp].x,8,'f',4).arg(v[icp].y,8,'f',4).arg(v[icp].z,8,'f',4)
        .arg(relaxed[icp].x,9,'f',2).arg(relaxed[icp].y,9,'f',2).arg(relaxed[icp].z,9,'f',2));
    };
    Globals::MainTextDisplay->append(QString("\n"));
    Globals::MainTextDisplay->setCurrentFont(previous);
    Globals::MainTextDisplay->append(QString("\n"));
    // after printing, we actually move the wake to the new position

    // ... TBD ...
    // the filaments can be moved, but
    // the wake stripes have to be deleted and recreated

    delete x;
    delete v;
    delete relaxed;
    */
  }
}

void VortexLatticeModel::analyzeCirculation()
{
  if (validSolution)
  {
    Globals::MainTextDisplay->append(QString("\n*** analyzing the local circulation ***\n"));
    // the lift direction is perpendicular to the free-flow direction
    Vector eLift = cross(vInfinity, Vector(0.0, 1.0, 0.0));
    Vector totalForce = Vector(0.0, 0.0, 0.0);
    int NumberOfWings = Globals::GeometryModel->numberOfWings();
    for (int wing=0; wing<NumberOfWings; wing++)
    {
      // A const-doublet panel (i.e a vortex ring) does not have a net lift.
      // We only have to take into account the vortex lines which spawn the wake panels.
      Vector wingForce = Vector(0.0, 0.0, 0.0);
      for (int iw=0; iw<NumberOfWakes; iw++)
      {
	WakeStripe *wk = wake->at(iw);
	// We sum up the lift for each wing separately.
	if (wk->getWingIndex() == wing+1)
	{
	  // get the start points of the wake stripe baseline
	  Vector x1 = wk->wakeStart(1);
	  Vector x2 = wk->wakeStart(2);
	  // span units are mm
	  Vector S = x2-x1;
	  // get the flow velocity vector at the wake control point
	  // v units are m/s
	  Vector v = vSolution[NumberOfPanels+iw];
	  // get the circulation strength of the associated horse-shoe vortex
	  // gamma units are m/s mm
	  double gamma = -muSolution[NumberOfPanels+iw];
          // The force generated by the vortex is perpendicular to the flow
	  // and the axis direction of the vortex
          // F = rho * v_F x Gamma * ds
	  // force units are (m/s)² mm²
	  Vector force = cross(v, S) * gamma;
	  wingForce += force;
	};
      };
      totalForce += wingForce;
      // total wing area in dm²
      double area = Globals::GeometryModel->getWing(wing+1)->wingArea();
      area *= 1e4; // convert dm² to mm²
      // the coefficient of lift is defined by
      // L = rho/2 * v_Inf^2 * c_L * area = F
      double cL = 2.0*dot(eLift,wingForce)/area;
      double cD = 2.0*dot(vInfinity,wingForce)/area;
      Globals::MainTextDisplay->append(
	QString("wing No. %1 : %2").arg(wing+1).arg(Globals::GeometryModel->getWing(wing+1)->getName()));
      Globals::MainTextDisplay->append(
	QString("wing area = %1").arg(area / 1e4, 8,'f',2) + QString::fromUtf8(" dm²"));
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
    // the total momentum about the origin (0,0,0) of the coordinate system
    Vector totalMomentum = Vector(0.0, 0.0, 0.0);
    // A const-doublet panel (i.e a vortex ring) does not have a net lift.
    // It has a momentum, though.
    for (int ipan=0; ipan<NumberOfPanels; ipan++)
    {
      // get the flow velocity vector at the panel control point
      // v units are m/s
      Vector v = vSolution[ipan];
      // the panel normal
      Vector n = normal[ipan];
      // the vortex strength
      double gamma = -muSolution[ipan];
      // the panel area
      double Area = mesh->at(ipan)->panelArea();
      totalMomentum += cross (v, n) * gamma * Area;
    };
    // Add the contribution of the true lift
    for (int iw=0; iw<NumberOfWakes; iw++)
    {
      WakeStripe *wk = wake->at(iw);
      // get the start points of the wake stripe baseline
      Vector x1 = wk->wakeStart(1);
      Vector x2 = wk->wakeStart(2);
      // span units are mm
      Vector S = x2-x1;
      // Vector from the origin tho the wake line center
      Vector R = (x1+x2)*0.5;
      // get the flow velocity vector at the wake control point
      // v units are m/s
      Vector v = vSolution[NumberOfPanels+iw];
      // get the circulation strength of the associated horse-shoe vortex
      // gamma units are m/s mm
      double gamma = -muSolution[NumberOfPanels+iw];
      // The force generated by the vortex is perpendicular to the flow
      // and the axis direction of the vortex
      // F = rho * v_F x Gamma * ds
      // force units are (m/s)² mm²
      Vector force = cross(v, S) * gamma;
      totalMomentum += cross(R,force);
    };
    // double RefLength = Globals::GeometryModel->getRefChord();
    Globals::MainTextDisplay->append(
      QString("total momentum about (0,0,0) : cM = (%1,  %2,  %3) mm")
	.arg(totalMomentum.x/RefArea).arg(totalMomentum.y/RefArea).arg(totalMomentum.z/RefArea));
    double cMy = totalMomentum.y/RefArea;
    Globals::MainTextDisplay->append(
      QString("center of lift xL = %1 mm").arg(-2.0*cMy/cL));
  };
}

void VortexLatticeModel::analyzeWake()
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

void VortexLatticeModel::getMuRange (double *min, double *max)
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

void VortexLatticeModel::printPaneling()
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
	.arg(ipan,6)
	.arg(pc.x,8,'f',2).arg(pc.y,8,'f',2).arg(pc.z,8,'f',2)
	.arg(pn.x,8,'f',4).arg(pn.y,8,'f',4).arg(pn.z,8,'f',4));
    };
    Globals::MainTextDisplay->setCurrentFont(previous);
    Globals::MainTextDisplay->append(QString("\n"));
  }
}

void VortexLatticeModel::printSolution()
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
    Globals::MainTextDisplay->append(QString("       \t    doublet  \t              residual    \t    \tvelocity"));
    Globals::MainTextDisplay->append(QString("       \t    strength \t              normal"));
    Globals::MainTextDisplay->append(QString("       \t             \t              velocity"));
    Globals::MainTextDisplay->append(QString("----------------------------------------------------------------------------------------------------------------------------------------"));
    for (int icp=0; icp<NumberCP; icp++)
    {
      Vector n = normal[icp];
      Vector v = vSolution[icp];
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

void VortexLatticeModel::printCirculation()
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

void VortexLatticeModel::sourceGammaPlot(vtkChartXY *chart)
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
      sprintf(label,"Wing No.%3d (VLM)",wing+1);
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

void VortexLatticeModel::sourceClCdPlot(vtkChartXY *chart)
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
      sprintf(label,"Wing No.%3d cl (VLM)",wing+1);
      arr1->SetName(label);
      table->AddColumn(arr1);
      vtkSmartPointer<vtkFloatArray> arr2 = vtkSmartPointer<vtkFloatArray>::New();
      sprintf(label,"Wing No.%3d cd (VLM)",wing+1);
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
