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

#include "wakestripe.h"

WakeStripe::WakeStripe(Streamline *A, Streamline *B)
{
  filamentA = A;
  filamentB = B;
  startA = filamentA->point(0);
  startB = filamentB->point(0);
  Vector firstA = filamentA->point(1);
  Vector firstB = filamentB->point(1);
  Vector d1 = firstA-startA;
  d1.normalize();
  Vector d2 = firstB-startB;
  d2.normalize();
  CP = (startA+startB)*0.5;
  CP += (d1+d2)*1e-4;
  wingIndex=0;		// default
  wingSpan=0.0;		// default
  wingChord=100.0;	// default
  NP = filamentA->segmentN();
  // one should issue a warning if computing a wake stripe
  // between filaments of different length
  if (filamentB->segmentN()<NP) NP=filamentB->segmentN();
  // generate panels connecting the streamlines
  panels = new FlatPanel*[NP];
  for (int ipan=0; ipan<NP; ipan++)
  {
    Vector C = filamentA->point(ipan);
    Vector D = filamentA->point(ipan+1);
    Vector E = filamentB->point(ipan);
    Vector F = filamentB->point(ipan+1);
    panels[ipan] = new FlatPanel(D,C,E,F);
  };
  // the first panel defines the wake normal
  Normal = panels[0]->panelNormal();
}

WakeStripe::~WakeStripe()
{
}

void WakeStripe::setWingIndex(int index)
{
  wingIndex = index;
}

int WakeStripe::getWingIndex()
{
  return wingIndex;
}

void WakeStripe::setWingSpanPos(double s)
{
  wingSpan = s;
}

double WakeStripe::getWingSpanPos()
{
  return wingSpan;
}

void WakeStripe::setWingChord(double c)
{
  wingChord = c;
}

double WakeStripe::getWingChord()
{
  return wingChord;
}

Vector WakeStripe::wakeStart(int i)
{
  Vector v = Vector(0.0, 0.0, 0.0);
  if (i==1) v = startA;
  if (i==2) v = startB;
  return v;
}

Vector WakeStripe::wakeInfinity(int i)
{
  Vector v = Vector(0.0, 0.0, 0.0);
  if (i==1) v = filamentA->point(filamentA->segmentN());
  if (i==2) v = filamentB->point(filamentB->segmentN());
  return v;
}

Vector WakeStripe::wakeCP()
{
  return CP;
}

Vector WakeStripe::wakeNormal()
{
  return Normal;
}

Vector WakeStripe::inducedVelocity(Vector P)
{
  Vector vind = Vector(0,0,0);
  for (int ipan=0; ipan<NP; ipan++)
  {
    // we need to create a copy of the panel for thread-safety
    FlatPanel *p = new FlatPanel(panels[ipan]);
    // FlatPanel *p = panels[ipan];
    p->ComputePIC(P);
    vind += p->ConstDoubletVelocity();
    delete p;
  };
  return vind;
}

double WakeStripe::inducedPotential(Vector P)
{
  double phi = 0.0;
  for (int ipan=0; ipan<NP; ipan++)
  {
    // we need to create a copy of the panel for thread-safety
    FlatPanel *p = new FlatPanel(panels[ipan]);
    // FlatPanel *p = panels[ipan];
    p->ComputePIC(P);
    phi += p->ConstDoubletPotential();
    delete p;
  };
  return phi;
}

