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

#include "geometrymodel.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>

// ******************************************************
// ************* Class Model ****************************
// ******************************************************

Model::Model()
{
  // set default global properties;
  setName("no name");
  setAuthor("no author");
  setMass(1000.0);
  setRefSpan(2000.0);
  setRefArea(50.0);
  setRefChord(100.0);
  setAOA(2.0);
  setXinfinity(10.0);
  setWakePanelNumber(2);
  // a new model initialy has one default wing
  NumberOfWings=1;
  wings = new GeometryWing*[NumberOfWings];
  GeometryWing *newwing = new GeometryWing(3);
  wings[0] = newwing;
  highlighting = NoHighlight;
  airfoilDatabase = new QList<Airfoil *>();
  NumberOfAirfoils=0;
}

Model::Model(QDomElement parent)
{
  QString s;
  QDomNode node, afnode;
  QDomElement elem, afelem;
  // read global properties
  this->setName(parent.attribute("name","no name"));
  this->setAuthor(parent.attribute("author","no author"));
  s = parent.attribute("mass","1000.0");
  this->setMass(s.toDouble());
  s = parent.attribute("refspan","2000.0");
  this->setRefSpan(s.toDouble());
  s = parent.attribute("refarea","50.0");
  this->setRefArea(s.toDouble());
  s = parent.attribute("refchord","100.0");
  this->setRefChord(s.toDouble());
  s = parent.attribute("xinfinity","10.0");
  this->setXinfinity(s.toDouble());
  s = parent.attribute("wakepanelnumber","2");
  this->setWakePanelNumber(s.toInt());
  s = parent.attribute("aoa","2.0");
  this->setAOA(s.toDouble());
  // browse nodes for wing entries
  NumberOfWings=0;
  node=parent.firstChild();
  while( !node.isNull() )
  {
    elem=node.toElement();
    if( !elem.isNull() )
    {
      if( elem.tagName() == "wing" ) NumberOfWings++;
    }
    node=node.nextSibling();
  };
  if (NumberOfWings>0)
  {
    wings = new GeometryWing*[NumberOfWings];
    int nw=0;
    // browse once more and read wings
    node=parent.firstChild();
    while( !node.isNull() )
    {
      elem=node.toElement();
      if( !elem.isNull() )
      {
        if( elem.tagName() == "wing" )
        {
          GeometryWing *newwing = new GeometryWing(elem);
          wings[nw++] = newwing;
	  Globals::MainTextDisplay->append(
	    QString("successfully read wing No. %1").arg(nw));
        }
      }
      node=node.nextSibling();
    };
  } else
  {
    // something's wrong - use default
    Globals::MainTextDisplay->append(
       QString("\n")+
       QString("error reading model data\n")+
       QString("no valid wing found - using default\n\n"));
    NumberOfWings=1;
    wings = new GeometryWing*[NumberOfWings];
    GeometryWing *newwing = new GeometryWing(3);
    wings[0] = newwing;
  };
  airfoilDatabase = new QList<Airfoil *>();
  NumberOfAirfoils=0;
  // start over with the document and search for the airfoil database
  node=parent.firstChild();
  while( !node.isNull() )
  {
    elem=node.toElement();
    if( !elem.isNull() )
    {
      if (elem.tagName() == "airfoildatabase")
      {
        Globals::MainTextDisplay->append("airfoil database found");
        // database found now read all the airfoils
        afnode = elem.firstChild();
        while( !afnode.isNull() )
        {
          afelem=afnode.toElement();
          if( !afelem.isNull() )
          {
            if (afelem.tagName() == "airfoil")
            {
              // read airfoil data
              Airfoil *newaf = new Airfoil(afelem);
              if (newaf->isValid())
              {
                Globals::MainTextDisplay->append(
                    QString("airfoil name= %1 : %2 points : index=%3").
                    arg(newaf->getName()).arg(newaf->numPoints()).arg(newaf->getIndex()));
                int idx = newaf->getIndex();
                // the index is normally set to the database position when storing the airfoil
                addAirfoil(newaf);
                // the index read from the file is restored for the consistency check
                newaf->setIndex(idx);
              } else {
                delete newaf;
              };
            };
          };
          afnode = afnode.nextSibling();
        };
      };
    };
    node=node.nextSibling();
  };
  // check if all references actually point into the database
  checkAirfoilReferences();
}

Model::~Model()
{
  for (int i=1; i<NumberOfWings; i++)
  {
    delete wings[i-1];
  };
  delete wings;
  // delete all airfoils
  int naf = numberOfAirfoils();
  for (int i=1; i<=naf; i++)
    delete airfoilDatabaseAt(i);
  delete airfoilDatabase;
}

void Model::update()
{
  for (int i=0; i<NumberOfWings; i++) wings[i]->update();
}

int Model::numberOfWings()
{
  return(NumberOfWings);
}

GeometryWing* Model::getWing(int n)
{
  if (n<1) return(wings[0]);
  if (n>NumberOfWings) return(wings[NumberOfWings-1]);
  return(wings[n-1]);
}

// handle global properties

QString Model::getName()
{
  return(ModelName);
}

void Model::setName(QString s)
{
  ModelName=s;
}

QString Model::getAuthor()
{
  return(ModelAuthor);
}

void Model::setAuthor(QString s)
{
  ModelAuthor=s;
}

double Model::getMass()
{
  return(ModelMass);
}

void Model::setMass(double m)
{
  ModelMass=m;
}

double Model::getRefSpan()
{
  return(ReferenceSpan);
}

void Model::setRefSpan(double span)
{
  ReferenceSpan=span;
}

double Model::getRefArea()
{
  return(ReferenceArea);
}

void Model::setRefArea(double area)
{
  ReferenceArea=area;
}

double Model::getRefChord()
{
  return(ReferenceChord);
}

void Model::setRefChord(double chord)
{
  ReferenceChord=chord;
}

double Model::getXinfinity()
{
  return(xInfinity);
}

void Model::setXinfinity(double xinf)
{
  xInfinity=xinf;
}

int Model::getWakePanelNumber()
{
  return WakePanelNumber;
}

void Model::setWakePanelNumber(int n)
{
  WakePanelNumber = n;
}

double Model::getAOA()
{
  return(flowAOA);
}

void Model::setAOA(double aoa)
{
  flowAOA=aoa;
}

// manage some wings belonging to the model

void Model::addWing()
{
  // recreate array of wings
  GeometryWing **newarray = new GeometryWing*[NumberOfWings+1];
  for (int i=0; i<NumberOfWings; i++)
    newarray[i] = wings[i];
  NumberOfWings+=1;
  // create a default wing with two stations
  GeometryWing *newwing = new GeometryWing(2);
  newarray[NumberOfWings-1] = newwing;
  // replace the old with the new array
  delete wings;
  wings = newarray;
}

void Model::deleteWing(int n)
{
  if ((NumberOfWings>1) && (n>=1) && (n<=NumberOfWings))
  {
    delete wings[n-1];
    for (int i=n; i<NumberOfWings; i++)
      wings[i-1] = wings[i];
    NumberOfWings-=1;
  };
}

// manage airfoil database belonging to the model

int Model::numberOfAirfoils()
{
  return(NumberOfAirfoils);
}

void Model::addAirfoil(Airfoil *af)
{
  airfoilDatabase->append(af);
  NumberOfAirfoils = airfoilDatabase->count();
  af->setIndex(NumberOfAirfoils);
}

Airfoil* Model::airfoilDatabaseAt(int index)
{
  return airfoilDatabase->at(index-1);
}

bool Model::airfoilIsUsed(int index)
{
  bool used=FALSE;
  for (int w=1; w<=NumberOfWings; w++)
  {
    GeometryWing *wing = getWing(w);
    int ns = wing->numberOfStations();
    for (int st=1; st<=ns; st++)
    {
      GeometryStation *station = wing->getStation(st);
      int id = station->getAirfoilIndex();
      if (id == index) used=TRUE;
    }
  };
  return(used);
}

void Model::checkAirfoilReferences()
{
  // test which indices correspond to indices present in the database
  // reset all others to zero
  // set the valid ones to the sequence number of the airfoil
  for (int w=1; w<=NumberOfWings; w++)
  {
    GeometryWing *wing = getWing(w);
    int ns = wing->numberOfStations();
    for (int st=1; st<=ns; st++)
    {
      GeometryStation *station = wing->getStation(st);
      int index = station->getAirfoilIndex();
      bool found = (index==0);
      for (int af=1; af<=NumberOfAirfoils; af++)
      {
        int readIndex = airfoilDatabaseAt(af)->getIndex();
        if (index==readIndex)
        {
          station->setAirfoilIndex(af);
          found=TRUE;
        };
      };
      if (!found)
      {
        station->setAirfoilIndex(0);
        Globals::MainTextDisplay->append(
          QString("didn't find referenced airfoil for wing No. %1, station No. %2").arg(w).arg(st));
      };
    }
  };
  // set the index value of the airfoil equal to the sequence number
  for (int af=1; af<=NumberOfAirfoils; af++)
  {
    airfoilDatabaseAt(af)->setIndex(af);
  };
}

void Model::deleteAirfoil(int index)
{
  Globals::MainTextDisplay->append(
      QString("removing airfoil No. %1 from database").arg(index));
  if ((index>0) && (index<=NumberOfAirfoils))
  {
    // delete the airfoil
    delete airfoilDatabaseAt(index);
    // remove the database entry
    airfoilDatabase->removeAt(index-1);
    NumberOfAirfoils = airfoilDatabase->count();
    // if one airfoil is removed, all references have to be shifted accordingly
    for (int w=0; w<NumberOfWings; w++)
    {
      GeometryWing *wing = wings[w];
      int ns = wing->numberOfStations();
      for (int st=1; st<=ns; st++)
      {
	GeometryStation *station = wing->getStation(st);
	int id = station->getAirfoilIndex();
	if (id == index) station->setAirfoilIndex(0);
	if (id > index) station->setAirfoilIndex(id-1);
      }
    };
  } else {
    Globals::MainTextDisplay->append(QString("wrong index - doing nothing"));
  };
}

void Model::sourceVTK(
      vtkAppendPolyData *polyData,
      int highlightWing,
      int highlightSegment,
      int highlightStation)
{
  polyData->SetNumberOfInputs(NumberOfWings);
  for (int i=0; i<NumberOfWings; i++)
  {
    if (i+1==highlightWing)
    {
      wings[i]->setHighlighting(highlighting);
    } else {
      wings[i]->setHighlighting(NoHighlight);
    };
    vtkPolyData *wingData = vtkPolyData::New();
    wings[i]->sourceVTK(
      wingData,
      highlightSegment, highlightStation);
    polyData->SetInputDataByNumber(i, wingData);
    wingData->Delete();
  };
}

void Model::setHighlighting (geometryHighlightMode mode)
{
  highlighting = mode;
}
    
void Model::writeGeometryXML(
        QDomDocument doc,
        QDomElement parent)
{
  // write global properties
  parent.setAttribute("name", this->getName());
  parent.setAttribute("author", this->getAuthor());
  parent.setAttribute("mass", QString("%1").arg(this->getMass(),9,'f',1));
  parent.setAttribute("refspan", QString("%1").arg(this->getRefSpan(),12,'f',3));
  parent.setAttribute("refarea", QString("%1").arg(this->getRefArea(),12,'f',4));
  parent.setAttribute("refchord", QString("%1").arg(this->getRefChord(),10,'f',3));
  parent.setAttribute("xinfinity", QString("%1").arg(this->getXinfinity(),12,'f',6));
  parent.setAttribute("wakepanelnumber", QString("%1").arg(this->getWakePanelNumber()));
  parent.setAttribute("aoa", QString("%1").arg(this->getAOA(),12,'f',4));
  // write wing data
  for (int i=0; i<NumberOfWings; i++)
  {
    wings[i]->writeGeometryXML(doc, parent);
  }
  QDomElement afDatabaseElement = doc.createElement("airfoildatabase");
  int naf = airfoilDatabase->size();
  for (int i=0; i<naf; i++)
  {
    Airfoil *af = airfoilDatabase->at(i);
    af->writeGeometryXML(doc, afDatabaseElement);
  };
  parent.appendChild(afDatabaseElement);
}
