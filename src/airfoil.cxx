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

#include "airfoil.h"

// OpenBLAS
// include "cblas.h"
#include "lapacke.h"

/******************************************************/
/****** Class Camberline*******************************/
/******************************************************/

Camberline::Camberline(int NumberOfPanels)
{
  npan = NumberOfPanels;
  xs = new double[npan+1];
  ys = new double[npan+1];
  aoa = new double[npan];
  p = new Vector[npan+1];
  // Theta is running from Pi down to offTheta.
  // the panel density increases at Theta=Pi where the derivative
  // of the cosine vanishes qudratically
  // at the traling edge (start and end of the loop)
  // the panel density is about 3 times the minimum density
  double offTheta = 0.05*M_PI;
  double dTheta = (M_PI - offTheta) / NumberOfPanels;
  // While Theta is running from Pi down to offTheta
  // the cosine is running from -1.0 to maxCos.
  // This is mapped to x=0 to x=1
  double maxCos=cos(offTheta);
  double Theta=offTheta;
  for (int i=0; i<=npan; i++)
  {
    xs[i] = (1.0 + cos(M_PI - (double)i*dTheta)) / (1.0+maxCos);
    ys[i] = 0.0;
    if (i<npan) aoa[i] = 0.0;
    p[i] = Vector(xs[i], 0.0, 0.0);
    Theta += dTheta;
  }
}

Camberline::Camberline(int NumberOfPanels,
             int NumberFlapPanels,
             double flapPosition)
{
  double theta, dTheta;
  npan = NumberOfPanels;
  xs = new double[npan+1];
  ys = new double[npan+1];
  aoa = new double[npan];
  p = new Vector[npan+1];
  // point index of the flap hinge
  int flapBreak = NumberOfPanels - NumberFlapPanels;
  double flapTheta = acos(1.0-2.0*flapPosition);
  // panels ahead of the hinge
  dTheta = flapTheta / (NumberOfPanels-NumberFlapPanels);
  theta=0.0;
  for (int i=0; i<flapBreak; i++)
  {
    xs[i] = 0.5*(1.0-cos(theta));
    ys[i] = 0.0;
    if (i<npan) aoa[i] = 0.0;
    p[i] = Vector(xs[i], 0.0, 0.0);
    theta+=dTheta;
  }
  // panels after the hinge
  dTheta = (M_PI-flapTheta) / NumberFlapPanels;
  for (int i=flapBreak; i<=npan; i++)
  {
    xs[i] = 0.5*(1.0-cos(theta));
    ys[i] = 0.0;
    if (i<npan) aoa[i] = 0.0;
    p[i] = Vector(xs[i], 0.0, 0.0);
    theta+=dTheta;
  }
}

Camberline::Camberline(Camberline *cl1, double f1,
                       Camberline *cl2, double f2)
{
  npan = cl1->numberPanels();
  if (npan == cl2->numberPanels())
  {
    xs = new double[npan+1];
    ys = new double[npan+1];
    aoa = new double[npan];
    p = new Vector[npan+1];
    double c1 = f1/(f1+f2);
    double c2 = f2/(f1+f2);
    for (int i=0; i<=npan; i++)
    {
      xs[i] = c1*cl1->pointX(i) + c2*cl2->pointX(i);
      ys[i] = c1*cl1->pointY(i) + c2*cl2->pointY(i);
      p[i] = cl1->pointVec(i)*c1 + cl2->pointVec(i)*c2;
    };
    for (int i=0; i<npan; i++) {
      aoa[i] = atan(-(ys[i+1]-ys[i])/(xs[i+1]-xs[i]));
    };
  } else {
    npan = 1;
    xs = new double[npan+1];
    ys = new double[npan+1];
    aoa = new double[npan];
    p = new Vector[npan+1];
    xs[0] = 0.0;
    xs[1] = 1.0;
    ys[0] = 0.0;
    ys[1] = 0.0;
    aoa[0] = 0.0;
    p[0] = Vector(0.0, 0.0, 0.0);
    p[1] = Vector(1.0, 0.0, 0.0);
  };
}

Camberline::~Camberline()
{
  delete[] xs;
  delete[] ys;
  delete[] aoa;
  delete[] p;
}

int Camberline::numberPoints()
{
  return(npan+1);
}

int Camberline::numberPanels()
{
  return(npan);
}

double Camberline::pointX(int i)
{
  double r;
  r = 0.0;
  if ((i>=0) && (i<=npan)) r=xs[i];
  return(r);
}

double Camberline::pointY(int i)
{
  double r;
  r = 0.0;
  if ((i>=0) && (i<=npan)) r=ys[i];
  return(r);
}

Vector Camberline::pointVec(int i)
{
  Vector X = Vector(0.0, 0.0, 0.0);
  if ((i>=0) && (i<=npan)) X=p[i];
  return(X);
}

void Camberline::setY(int i, double y)
{
  if ((i>=0) && (i<=npan))
    ys[i] = y;
}

void Camberline::computePoints(Airfoil *af)
{
  int i, n;
  double *yp = new double[10];
  double ycl;
  ys[0]=0.0;
  ys[npan]=0.0;
  for (i=1; i<npan; i++) {
    ycl=0.0;
    n=af->splx->findInverse(xs[i], yp, 10);
    if (n==2) ycl=0.5*(af->sply->f(yp[0])+af->sply->f(yp[1]));
    ys[i]=ycl;
  }
}

void Camberline::computeCoeff(double *cl0, double *alfa0,
                              double* cm0, double* asf, double *clsf)
{
  int i, k;
  double *xvortex = new double[npan];
  double *xcontrol = new double[npan];
  double gamma0, dclda;

  for (i=0; i<npan; i++) {
    xvortex[i] = 0.75*xs[i] + 0.25*xs[i+1];
    xcontrol[i] = 0.25*xs[i] + 0.75*xs[i+1];
    aoa[i] = atan(-(ys[i+1]-ys[i])/(xs[i+1]-xs[i]));
    // printf(" dx=%9.6f  a=%9.6f deg\n",xs[i+1]-xs[i],aoa[i]*180.0/M_PI);
  }
  // The inflence matrix M is computed in a row-major order.
  // The rows contain the coefficients for all panels to a single control point.
  // This way, the control point index is the row number and the outer loop index.
  // The panels index is the column number and the inner loop index.
  // The induced velocity of any vortex is proportional
  // to the inverse distance to the respective control point.
  double *indu = new double[npan*npan];
  for (i=0; i<npan; i++)          // lines (control points)
    for (k=0; k<npan; k++)        // columns (vortices)
      indu[npan*i+k] = 1.0/(xcontrol[i]-xvortex[k]);

  // Now we solve the system of linear equations
  // Indu * gamma = aoa
  // Indu * gamma = 1    (for the slope of the lift coeffiecient)
  // for this we create a right-hand side matrix with 2 vectors
  double *b = new double[2*npan];
  for (i=0; i<npan; i++)          // lines (control points)
  {
    b[2*i] = aoa[i];
    b[2*i+1] = 1.0;
  }

  // use OpenBLAS LAPACKE to solve the system
  int *ipiv = new int[npan];
  int info = LAPACKE_dgesv(
      LAPACK_ROW_MAJOR,		// storage ordering of the matrix
      npan,			// LDA : the leading array dimension = number of rows
      2,			// number of right-hand side vectors
      indu,			// the matrix
      npan,			// LDB : the leading dimension of the RHS vectors
      ipiv,			// workspace for the pivot vector
      b,			// right-hand side matrix
      2);			// number of right-hand side vectors

  if (info==0)
  {
    *cl0=0.0;
    *cm0=0.0;
    dclda=0.0;
    for (i=0; i<npan; i++) {
      *cl0+=2.0*b[2*i]*2.0*M_PI;
      *cm0+=2.0*(0.25-xvortex[i])*b[2*i]*2.0*M_PI;;
      dclda+=2.0*b[2*i+1]*2.0*M_PI;
    };
    gamma0=b[0];
    *alfa0 = -*cl0 / dclda;
    // suction-free entrance
    // zero vorticity at leading edge
    *asf = -gamma0 / b[1];
    *clsf = *cl0 + dclda * *asf;
  }
  else
  {
    Globals::MainTextDisplay->append(
	QString("error (linear algebra) computing thin-airfoil theory parameters"));
    *cl0=0.0;
    *cm0=0.0;
    *alfa0 = 0.0;
    *asf = 0.0;
    *clsf = 0.0;
  };

  delete[] xvortex;
  delete[] xcontrol;
  delete[] indu;
  delete[] b;
  delete[] ipiv;
}

void Camberline::mapToVectors(Vector nose, Vector end, Vector up, double stretch)
{
  Vector xvec = end-nose;
  Vector yvec = up;
  yvec.normalize();
  yvec*=xvec.norm()*stretch;
  for (int i=0; i<=npan; i++)
    p[i] = nose + xvec*xs[i] + yvec*ys[i];
}

Vector Camberline::wakeDirection()
{
  return( p[npan] - p[npan-1] );
}

/******************************************************/
/****** Class Outline*******************************/
/******************************************************/

Outline::Outline(int NumberOfPanels)
{
  npan = NumberOfPanels;
  xs = new double[npan+1];
  ys = new double[npan+1];
  p = new Vector[npan+1];
  // Theta is running from offTheta to 2PI-offTheta
  // the panel density increases at Theta=Pi where the derivative
  // of the cosine vanishes qudratically
  // at the traling edge (start and end of the loop)
  // the panel density is about 3 times the minimum density
  double offTheta = 0.05*M_PI;
  double dTheta = (2.0*M_PI - 2*offTheta) / NumberOfPanels;
  // the cosine is running from maxCos down to -1.0
  // this is mapped to x=1 down to x=0
  double maxCos=cos(offTheta);
  double Theta=offTheta;
  for (int i=0; i<=npan; i++)
  {
    xs[i] = (cos(Theta)+1.0) / (maxCos+1.0);
    ys[i] = 0.0;
    p[i] = Vector(xs[i], 0.0, 0.0);
    Theta += dTheta;
  }
}

Outline::Outline(Outline *ol1, double f1,
		 Outline *ol2, double f2)
{
  npan = ol1->numberPanels();
  if (npan == ol2->numberPanels())
  {
    xs = new double[npan+1];
    ys = new double[npan+1];
    p = new Vector[npan+1];
    double c1 = f1/(f1+f2);
    double c2 = f2/(f1+f2);
    for (int i=0; i<=npan; i++)
    {
      xs[i] = c1*ol1->pointX(i) + c2*ol2->pointX(i);
      ys[i] = c1*ol1->pointY(i) + c2*ol2->pointY(i);
      p[i] = ol1->pointVec(i)*c1 + ol2->pointVec(i)*c2;
    };
  } else {
    npan = 2;
    xs = new double[npan+1];
    ys = new double[npan+1];
    p = new Vector[npan+1];
    xs[0] = 1.0;
    xs[1] = 0.0;
    xs[2] = 1.0;
    ys[0] = 0.0;
    ys[1] = 0.0;
    ys[2] = 0.0;
    p[0] = Vector(0.0, 0.0, 0.0);
    p[1] = Vector(1.0, 0.0, 0.0);
    p[2] = Vector(0.0, 0.0, 0.0);
  };
}

Outline::~Outline()
{
  delete[] xs;
  delete[] ys;
  delete[] p;
}

int Outline::numberPoints()
{
  return(npan+1);
}

int Outline::numberPanels()
{
  return(npan);
}

double Outline::pointX(int i)
{
  double r;
  r = 0.0;
  if ((i>=0) && (i<=npan)) r=xs[i];
  return(r);
}

double Outline::pointY(int i)
{
  double r;
  r = 0.0;
  if ((i>=0) && (i<=npan)) r=ys[i];
  return(r);
}

void Outline::setY(int i, double y)
{
  if ((i>=0) && (i<=npan))
    ys[i] = y;
}

void Outline::computePoints(Airfoil *af)
{
  int i, n;
  int np2;			// number of panels on upper side
  double yt;
  double *t = new double[10];	// spline parameters for inverse
  double yol;
  ys[0]=0.0;
  ys[npan]=0.0;
  np2= npan/2;
  for (i=1; i<=np2; i++)	// loop over upper side
  {
    n=af->splx->findInverse(xs[i], t, 10);
    if (n<1)
    {
      yol=0.0; 			// no point found
    } else {
      yol=af->sply->f(t[0]);
      for (int k=1; k<n; k++)	// take largest value for upper side
      {
	yt=af->sply->f(t[k]);
	if (yt>yol) yol=yt;
      };
    };
    ys[i]=yol;
  }
  for (i=np2+1; i<npan; i++)	// loop over lower side
  {
    n=af->splx->findInverse(xs[i], t, 10);
    if (n<1)
    {
      yol=0.0; 			// no point found
    } else {
      yol=af->sply->f(t[0]);
      for (int k=1; k<n; k++)	// take lowest value for lower side
      {
	yt=af->sply->f(t[k]);
	if (yt<yol) yol=yt;
      };
    };
    ys[i]=yol;
  }
}

void Outline::mapToVectors(Vector nose, Vector end, Vector up, double stretch)
{
  Vector xvec = end-nose;
  Vector yvec = up;
  yvec.normalize();
  yvec*=xvec.norm()*stretch;
  for (int i=0; i<=npan; i++)
    p[i] = nose + xvec*xs[i] + yvec*ys[i];
}

Vector Outline::pointVec(int i)
{
  Vector X = Vector(0.0, 0.0, 0.0);
  if ((i>=0) && (i<=npan)) X=p[i];
  return(X);
}

Vector Outline::wakeDirection()
{
  return((p[0]+p[npan]-p[1]-p[npan-1])*0.5);
}

/******************************************************/
/****** Class Airfoil *********************************/
/******************************************************/

Airfoil::Airfoil(int NumberOfPoints,
          double *x,
          double *y)
{
  name = new QString("<no name>");
  np = NumberOfPoints;
  xp = new double[np];
  yp = new double[np];
  for (int i=0; i<np; i++) {
    xp[i]=x[i];
    yp[i]=y[i];
  };
  if (doSanityCheck())
  {
    splx = new Spline(np,xp);
    sply = new Spline(np,yp);
    outline = new Spline2D(splx,sply);
    valid=TRUE;
  } else {
    Globals::MainTextDisplay->append(QString("airfoil sanity check failed"));
    valid=FALSE;
  };
  if (valid) {
    normalize();
    if (!valid)
      Globals::MainTextDisplay->append(QString("airfoil normalization failed"));
    computeOutlineProperties();
  };
}

Airfoil::Airfoil(QDomElement af)
{
  QDomNode node;
  QDomElement pt;
  // read airfoil attributes
  QString s = af.attribute("np","0");
  np=s.toInt();
  if (np<0) np=0;
  if (np>999) np=999;
  s = af.attribute("name","<noname>");
  name = new QString(s);
  s = af.attribute("index","0");
  databaseIndex=s.toInt();
  // now read points
  xp = new double[np];
  yp = new double[np];
  node = af.firstChild();
  for (int i=0; i<np; i++)
  {
    if (!node.isNull())
    {
      pt=node.toElement();
      if (!pt.isNull())
        if (pt.tagName() == "point")
        {
          s = pt.attribute("x","0.0");
          xp[i] = s.toDouble();
          s = pt.attribute("y","0.0");
          yp[i] = s.toDouble();
        };
    };
    node=node.nextSibling();
  };
  valid=TRUE;
  if (doSanityCheck())
  {
    splx = new Spline(np,xp);
    sply = new Spline(np,yp);
    outline = new Spline2D(splx,sply);
    this->normalize();
    if (!valid)
      Globals::MainTextDisplay->append(QString("airfoil normalization failed"));
    if (!doSanityCheck())
    {
      Globals::MainTextDisplay->append(QString("airfoil sanity check failed after normalization"));
      valid=FALSE;
    } else
    {
      Globals::MainStatusBar->showMessage(QString("airfoil successfully read : ")+this->getName());
      valid=TRUE;
    };
  } else {
    Globals::MainTextDisplay->append(QString("airfoil sanity check failed"));
    valid=FALSE;
  };
}

Airfoil::Airfoil(const QString &fileName)
{
  double x, y;
  double xx[999];
  double yy[999];
  double cl0, alfa0, cm0, asf, clsf;

  Globals::MainTextDisplay->append("reading airfoil : " + fileName);
  QFile file(fileName);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    Globals::MainStatusBar->showMessage(QString("cannot read : "+fileName+" - "+file.errorString()));
    return;
  }
  QTextStream in(&file);
  name = new QString(in.readLine());
  Globals::MainTextDisplay->append(this->getName());

  np=0;
  do {
    in >> x;
    in >> y;
    xx[np]=x;
    yy[np]=y;
    in.readLine();
    np++;
  } while (!in.atEnd());
  Globals::MainTextDisplay->append(QString::number(np) + QString(" points"));

  xp = new double[np];
  yp = new double[np];
  for (int i=0; i<np; i++) {
    xp[i]=xx[i];
    yp[i]=yy[i];
  };

  if (doSanityCheck())
  {
    splx = new Spline(np,xp);
    sply = new Spline(np,yp);
    outline = new Spline2D(splx,sply);
    this->normalize();
    if (!doSanityCheck())
    {
      Globals::MainTextDisplay->append(QString("airfoil sanity check failed after normalization"));
      valid=FALSE;
    } else
    {
      Camberline *cl = new Camberline(200);
      cl->computePoints(this);
      cl->computeCoeff(&cl0, &alfa0, &cm0, &asf, &clsf);
      Globals::MainTextDisplay->append(QString("cl0 = ") + QString::number(cl0,'f',3));
      Globals::MainTextDisplay->append(QString("alpha0 = ") +
	QString::number(180.0/M_PI*alfa0,'f',3) + QString(" deg"));
      Globals::MainTextDisplay->append(QString("cm0 = ") + QString::number(cm0,'f',4));
      Globals::MainTextDisplay->append(QString("alpha(sf) = ") +
	QString::number(180.0/M_PI*asf,'f',3) + QString(" deg"));
      Globals::MainTextDisplay->append(QString("cl(sf) = ") + QString::number(clsf,'f',3));
      Globals::MainStatusBar->showMessage(QString("airfoil successfully read : ")+this->getName());
      Globals::MainTextDisplay->append(QString("\n"));
      valid=TRUE;
    };
  } else {
    Globals::MainTextDisplay->append(QString("airfoil sanity check failed"));
    valid=FALSE;
  };
}

Airfoil::~Airfoil()
{
  delete[] xp;
  delete[] yp;
  delete name;
  if (valid)
  {
    delete splx;
    delete sply;
    delete outline;
  };
}

int Airfoil::numPoints()
{
  return(np);
}

double Airfoil::pointX(int i)
{
  if ((i>=0) && (i<np))
    return(xp[i]);
  else
    return(0.0);
}

double Airfoil::pointY(int i)
{
  if ((i>=0) && (i<np))
    return(yp[i]);
  else
    return(0.0);
}

bool Airfoil::doSanityCheck()
{
  bool isSane=TRUE;
  if (fabs(xp[0]-1.0)>0.02) isSane=FALSE;
  if (fabs(yp[0])>0.02) isSane=FALSE;
  if (fabs(xp[np-1]-1.0)>0.02) isSane=FALSE;
  if (fabs(yp[np-1])>0.02) isSane=FALSE;
  if (!isSane)
  {
    Globals::MainTextDisplay->append(
      QString("first or last point not near the expected trailing edge (1,0)"));
    return(isSane);
  };
  if (fabs(xp[0]-xp[np-1])>0.001)
  {
    isSane=FALSE;
    Globals::MainTextDisplay->append(
      QString("upper and lower side have different length of the trailing edge"));
    return(isSane);
  };
  double xmin=xp[0];
  double ymin=yp[0];
  for (int i=1; i<np; i++)
    if (xp[i]<xmin)
    {
      xmin=xp[i];
      ymin=yp[i];
    };
  if (fabs(xmin)>0.02) isSane=FALSE;
  if (fabs(ymin)>0.02) isSane=FALSE;
  if (!isSane)
  {
    Globals::MainTextDisplay->append(
      QString("didn't find a nose point near (0,0)"));
    return(isSane);
  };
  computeOutlineProperties();
  if ((circumference<1.95) || (circumference>2.3))
  {
    isSane=FALSE;
    Globals::MainTextDisplay->append(
      QString("circumference=%1 - outside the expected range").arg(circumference));
    return(isSane);
  };
  if ((area<0.02) || (area>0.2))
  {
    isSane=FALSE;
    Globals::MainTextDisplay->append(
      QString("airfoil should have a reasonable cross section - area=%1").arg(area));
    return(isSane);
  };
  // if (isSane) Globals::MainTextDisplay->append(QString("airfoil sanity check OK"));
  return(isSane);
}

void Airfoil::computeOutlineProperties()
{
  double x1=xp[0];
  double y1=yp[0];
  circumference=0.0;
  area=0.0;
  for (int i=1; i<np; i++)
  {
    double x2=xp[i];
    double y2=yp[i];
    circumference+=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    area+=(x1-x2)*(y1+y2)*0.5;
    x1=x2;
    y1=y2;
  };
  counterclockwise=(area<0.0);
  // area=fabs(area);
};

void Airfoil::normalize()
{
  int i, nex;
  double noset, nosex, nosey;
  double r, xte, yte;

  double *te = new double[2*np];
  nex = splx->findExtrema(te, 2*np);
  if (nex<1)
  {
    Globals::MainTextDisplay->append(
      QString("cannot normalize - didn't find airfoil nose point"));
    valid = FALSE;
  } else {
    nosex=1.0;
    noset=0.0;
    nosey=0.0;
    for (i=0; i<nex; i++) {
      if (splx->f(te[i])<nosex) {
	noset=te[i];
	nosex=splx->f(noset);
	nosey=sply->f(noset);
      };
    };

    xte=xp[0];
    if (xp[np-1]>xte) {
      xte=xp[np-1];
      r=(xte-xp[1])/(xp[0]-xp[1]);
      xp[0]=xp[1]+r*(xp[0]-xp[1]);
      yp[0]=yp[1]+r*(yp[0]-yp[1]);
    } else {
      r=(xte-xp[np-2])/(xp[np-1]-xp[np-2]);
      xp[np-1]=xp[np-2]+r*(xp[np-1]-xp[np-2]);
      yp[np-1]=yp[np-2]+r*(yp[np-1]-yp[np-2]);
    };
    yte=(yp[0]+yp[np-1])*0.5;

    for (i=0; i<np; i++) {
      xp[i]-=nosex;
      xp[i]/=(xte-nosex);
      yp[i]-=nosey;
      yp[i]-=(yte-nosey)*xp[i];
    };

    delete splx;
    delete sply;
    delete outline;
    splx = new Spline(np,xp);
    sply = new Spline(np,yp);
    outline = new Spline2D(splx,sply);
  };
}

QString Airfoil::getName()
{
  return(*name);
}

bool Airfoil::isValid()
{
  return(valid);
}

void Airfoil::writeGeometryXML(
        QDomDocument doc,
        QDomElement parent)
{
  QDomElement afElement = doc.createElement("airfoil");
  afElement.setAttribute("np", np);
  afElement.setAttribute("name", *name);
  afElement.setAttribute("index", databaseIndex);
  for (int i=0; i<np; i++)
  {
    QDomElement pointElement = doc.createElement("point");
    pointElement.setAttribute("x", QString("%1").arg(xp[i],9,'f',6));
    pointElement.setAttribute("y", QString("%1").arg(yp[i],9,'f',6));
    afElement.appendChild(pointElement);
  };
  parent.appendChild(afElement);
}

void Airfoil::setIndex(int idx)
{
  databaseIndex = idx;
}

int Airfoil::getIndex()
{
  return(databaseIndex);
}

