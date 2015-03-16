
/************************************************************************/
/*                                                                      */
/*  written by  U.Lehnert                               08/2012         */
/*                                                                      */
/************************************************************************/

#include "flatpanel.h"
#include <iostream>

char   message[80];

FlatPanel::FlatPanel(FlatPanel *p)
{
  NP = p->panelNP();
  CP = new Vector[NP];
  for (int i=0; i<NP; i++) CP[i] = p->panelPoint(i);
  init();
}

FlatPanel::FlatPanel(int n, Vector* pp)
{
  NP = n;
  CP = new Vector[NP];
  for (int i=0; i<NP; i++) CP[i] = pp[i];
  init();
}

FlatPanel::FlatPanel(Vector A, Vector B, Vector C)
{
  NP = 3;
  CP = new Vector[3];
  CP[0] = A;
  CP[1] = B;
  CP[2] = C;
  init();
}

FlatPanel::FlatPanel(Vector A, Vector B, Vector C, Vector D)
{
  NP = 4;
  CP = new Vector[4];
  CP[0] = A;
  CP[1] = B;
  CP[2] = C;
  CP[3] = D;
  init();
}

void FlatPanel::init()
{
  int i;
  double xi1, xi2, eta1, eta2;

  ds = new Vector[NP];
  nue = new Vector[NP];

  for (i=0; i<NP-1; i++) ds[i] = CP[i+1]-CP[i];
  ds[NP-1] = CP[0]-CP[NP-1];
  for (i=0; i<NP; i++) ds[i].normalize();
  Center = CP[0]; for (i=1; i<NP; i++) Center += CP[i];
  Center /= NP;
  // The cross product of panel lines defines the panel normal.
  // We use all panel lines to average the normal for non-flat panels.
  Vector r1 = CP[NP-1] - Center;
  Vector r2 = CP[0] - Center;
  Vector n0 = cross(r1, r2);
  n0.normalize();
  Normal = n0;
  for (i=1; i<NP; i++)
  {
    r1 = CP[i-1] - Center;
    r2 = CP[i] - Center;
    Vector n = cross(r1, r2);
    n.normalize();
    Normal += n;
  };
  Normal.normalize();

  // test if coplanar
  isPlanar=1;
  for (i=0; i<NP; i++)
  {
    if (dot(CP[i]-Center,Normal)>1e-2) isPlanar=0;
  };

  for (i=0; i<NP; i++) nue[i] = cross(ds[i], Normal);

  // define the panel coordinate system
  if (Normal.x*Normal.x < Normal.y*Normal.y) {
    Eta=cross(Normal,Vector(1,0,0));
    Eta.normalize();
    Xi=cross(Eta,Normal);
    Xi.normalize();
  } else {
    Xi=cross(Vector(0,1,0),Normal);
    Xi.normalize();
    Eta=cross(Normal,Xi);
    Eta.normalize();
  };

  // compute Area
  Area = 0.0;
  for (i=0; i<NP-1; i++)
  {
    xi1=dot(CP[i],Xi);
    xi2=dot(CP[i+1],Xi);
    eta1=dot(CP[i],Eta);
    eta2=dot(CP[i+1],Eta);
    Area+=0.5*(xi2+xi1)*(eta2-eta1);
  };
  xi1=dot(CP[NP-1],Xi);
  xi2=dot(CP[0],Xi);
  eta1=dot(CP[NP-1],Eta);
  eta2=dot(CP[0],Eta);
  Area+=0.5*(xi2+xi1)*(eta2-eta1);

  // allocate all memory needed for PIC computation
  R = new Vector[NP];
  r = new double[NP];
  vv = new Vector[NP];
  a = new double[NP];
  aa = new double[NP];
  g = new double[NP];
  l1 = new double[NP];
  l2 = new double[NP];

  F111 = new double[NP];
}

FlatPanel::~FlatPanel()
{
  delete[] R;
  delete[] r;
  delete[] vv;
  delete[] a;
  delete[] aa;
  delete[] g;
  delete[] l1;
  delete[] l2;
  delete[] F111;
  delete[] ds;
  delete[] nue;
  delete[] CP;
}

int FlatPanel::panelNP()
{
  return(NP);
}

Vector FlatPanel::panelPoint(int i)
{
  return(CP[i]);
}

Vector FlatPanel::panelCenter()
{
  return(Center);
}

Vector FlatPanel::panelNormal()
{
  return(Normal);
}

Vector FlatPanel::panelXi()
{
  return(Xi);
}

Vector FlatPanel::panelEta()
{
  return(Eta);
}

int FlatPanel::panelIsPlanar()
{
  return(isPlanar);
}

double FlatPanel::panelArea()
{
  return(Area);
}

Vector FlatPanel::panelNue(int i)
{
  return(nue[i]);
}

int FlatPanel::panelPointIsInterior()
{
  return(isInterior);
}

double FlatPanel::panelPointh()
{
  return(h);
}

double FlatPanel::panelPointa(int i)
{
  return(a[i]);
}

double FlatPanel::panelPointg(int i)
{
  return(g[i]);
}

double FlatPanel::panelH111()
{
  return(H111);
}

double FlatPanel::panelH113()
{
  return(hH113);
}

double FlatPanel::panelF111(int i)
{
  if ((i>=0) && (i<NP))
    return(F111[i]);
  else
    return(0.0);
}

Vector FlatPanel::ConstDoubletVelocity()
{
  return(vConstDoublet);
}

double FlatPanel::ConstDoubletPotential()
{
  return(fConstDoublet);
}

Vector FlatPanel::ConstSourceVelocity()
{
  return(vConstSource);
}

double FlatPanel::ConstSourcePotential()
{
  return(fConstSource);
}

// helper function which computes ArcTan(x/y)
// and gives valid results also for y=0
// (y is known positive or zero)
// derived from a series expansion about y/x=0
// may be identical to the atan2 function -> to be checked
double ATanMod(double x, double y)
{
  double yx, res;
  double ax = fabs(x);
  // cout << "AtanMod x=" << x <<"  y=" <<y;
  if (y >= 0.001 * ax)
  {
    // cout << "   normal";
    res = atan(x/y);
  }
  else
  {
    // cout << "   expansion";
    yx = y/x;
    if (x > 0.0)
      res = 0.5*M_PI + yx*(yx*yx/3.0-1.0);
    else
      res = -0.5*M_PI + yx*(yx*yx/3.0-1.0);
  };
  // cout << "   res=" << res << "   ref=" << atan(x/y) << "\n";
  return(res);
}

void FlatPanel::ComputePIC(Vector P)
{
  int i;
  double I1, I2;

  // cout << "compute PIC\n";
  // compute all geometric parameters
  for (i=0; i<NP; i++)
  {
    R[i] = CP[i] - P;
    r[i]=R[i].norm();
  };

  h = dot(P-Center,Normal);
  ha = fabs(h);
  // cout << "h = " << h << "\n";
  // cout << "ha = " << ha << "\n";

  // g² = a² + h²
  isInterior=1;
  for (i=0; i<NP; i++)
  {
    vv[i] = cross(ds[i],R[i]);
    a[i] = dot(Normal,vv[i]);
    g[i] = vv[i].norm();
    // cout << "a = " << a[i] << "\n";
    // cout << "g = " << g[i] << "\n";
    if (a[i]>0.0) isInterior=0;
    aa[i] = fabs(a[i]);
  };

  for (i=0; i<NP-1; i++)
  {
    l1[i] = dot(R[i],ds[i]);
    l2[i] = dot(R[i+1],ds[i]);
    // cout << "l1 = " << l1[i] << "\n";
    // cout << "l2 = " << l2[i] << "\n";
  };
  l1[NP-1] = dot(R[NP-1],ds[NP-1]);
  l2[NP-1] = dot(R[0],ds[NP-1]);
  // cout << "l1 = " << l1[NP-1] << "\n";
  // cout << "l2 = " << l2[NP-1] << "\n";

  // compute H(1,1,3) integral
  hH113 = 0.0;
  for (i=0; i<NP-1; i++)
  {
    I1 = ATanMod(l2[i], aa[i]) - ATanMod(l1[i], aa[i]);
    I2 = ATanMod(ha*l2[i], aa[i]*r[i+1]) - ATanMod(ha*l1[i], aa[i]*r[i]);
    // cout << "I1 = " << I1 << "\n";
    // cout << "I2 = " << I2 << "\n";
    if (a[i]>=0.0)
      hH113 -= I1-I2;
    else
      hH113 += I1-I2;
  };
  I1 = ATanMod(l2[NP-1], aa[NP-1]) - ATanMod(l1[NP-1], aa[NP-1]);
  I2 = ATanMod(ha*l2[NP-1], aa[NP-1]*r[0]) - ATanMod(ha*l1[NP-1], aa[NP-1]*r[NP-1]);
  // cout << "I1 = " << I1 << "\n";
  // cout << "I2 = " << I2 << "\n";
  if (a[NP-1]>=0.0)
    hH113 -= I1-I2;
  else
    hH113 += I1-I2;

  // compute F(1,1,1) integrals
  // from the definition there is always l2>=l1
  for (i=0; i<NP-1; i++)
  {
    if (l1[i] >= 0.0) {
      F111[i] = log((r[i+1]+l2[i])/(r[i]+l1[i]));
    } else {
      if (l2[i] < 0.0) {
        F111[i] = log((r[i]-l1[i])/(r[i+1]-l2[i]));
      } else {
        F111[i] = log((r[i]-l1[i])*(r[i+1]+l2[i])/(g[i]*g[i]));
      }
    }
  };
  if (l1[NP-1] >= 0.0) {
    F111[NP-1] = log((r[0]+l2[NP-1])/(r[NP-1]+l1[NP-1]));
  } else {
    if (l2[NP-1] < 0.0) {
      F111[NP-1] = log((r[NP-1]-l1[NP-1])/(r[0]-l2[NP-1]));
    } else {
      F111[NP-1] = log((r[NP-1]-l1[NP-1])*(r[0]+l2[NP-1])/(g[NP-1]*g[NP-1]));
    }
  };

  H111 = -ha*hH113;
  for (i=0; i<NP; i++)
    H111 -= a[i]*F111[i];

  // const. strength doublet distribution over the panel area
  // induced velicities and potentials
  vConstDoublet = Vector(0.0, 0.0, 0.0);
  for (i=0; i<NP-1; i++)
  {
    // cout << "g[i] = " << g[i] << "   cos2-cos1 = " <<  l2[i]/r[i+1]-l1[i]/r[i] << "\n";
    if (g[i] > 1e-6)
      vConstDoublet += vv[i]*over4Pi* (l2[i]/r[i+1]-l1[i]/r[i]) / (g[i]*g[i]);
  };
  // cout << "g[i] = " << g[NP-1] << "   cos2-cos1 = " <<  l2[NP-1]/r[0]-l1[NP-1]/r[NP-1] << "\n";
  if (g[NP-1] > 1e-6)
    vConstDoublet += vv[NP-1]*over4Pi* (l2[NP-1]/r[0]-l1[NP-1]/r[NP-1]) / (g[NP-1]*g[NP-1]);
  if (h>=0.0)
    fConstDoublet = over4Pi*hH113;
  else
    fConstDoublet = -over4Pi*hH113;

  // const. strength source distribution over the panel area
  // induced velicities and potentials
  fConstSource = -over4Pi*H111;
  vConstSource = Normal*fConstDoublet;
  for (i=0; i<NP; i++)
    vConstSource += nue[i]*over4Pi*F111[i];

}

