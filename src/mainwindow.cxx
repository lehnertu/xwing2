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

#include "mainwindow.h"
#include <QLocale>

#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include "vtkTextProperty.h"
#include "vtkContextScene.h"
#include "vtkChartXY.h"

MainWindow::MainWindow(QMainWindow *parent)
    : QMainWindow(parent)
{

    // setup user interface
    this->setLocale(QLocale::C);
    ui.setupUi(this);
    // this->setGeometry(200,100,1400,900);

    connect( ui.actionAbout, SIGNAL( triggered() ), this, SLOT( HelpAbout() ) );
    connect( ui.actionNew, SIGNAL( triggered() ), this, SLOT( FileNew() ) );
    connect( ui.actionOpen, SIGNAL( triggered() ), this, SLOT( FileOpen() ) );
    connect( ui.actionSave, SIGNAL( triggered() ), this, SLOT( FileSave() ) );
    connect( ui.actionSaveAs, SIGNAL( triggered() ), this, SLOT( FileSaveAs() ) );
    connect( ui.actionQuit, SIGNAL( triggered() ), this, SLOT( close() ) );
    connect( ui.actionSaveTextDisplay, SIGNAL( triggered() ), this, SLOT( SaveTextDisplay() ) );
    connect( ui.actionClearTextDisplay, SIGNAL( triggered() ), ui.TextDisplay, SLOT( clear() ) );
    connect( ui.actionLoadAirfoil, SIGNAL( triggered() ), this, SLOT ( LoadAirfoil() ) );
    connect( ui.actionDeleteAirfoil, SIGNAL( triggered() ), this, SLOT ( DeleteAirfoil() ) );
    connect( ui.actionExportSTL, SIGNAL( triggered() ), this, SLOT ( ExportSTL() ) );
    ui.statusbar->showMessage(QString("this is XWing Version 2.0"));

    // slots mit Namen, die Signalen der enstprechenden UI-Widgets entsprechen,
    // werden automatisch beim Programmstart mit diesen verknüpft

    // create a separate window for rendering the visualization
    graphWindow = new QWidget(this, Qt::Window | Qt::CustomizeWindowHint |
      Qt::WindowTitleHint | Qt::WindowSystemMenuHint );
    graphWindow->setWindowTitle(QString("XWing2 - rendering window"));
    graphWindow->setGeometry(100,50,1050,800);
    QVBoxLayout *graphLayout = new QVBoxLayout;
    graphWidget = new QVTKWidget(graphWindow);
    graphWidget->setMinimumSize(400,300);
    graphWidget->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
    graphLayout->addWidget(graphWidget);
    graphWindow->setLayout(graphLayout);
    graphWindow->show();

    // create an data structure, mapper and actor
    // for the model geometry (skeleton) display
    skeleton = vtkSmartPointer<vtkAppendPolyData>::New();
    skeleton->UserManagedInputsOn();
    skeletonMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    skeletonMapper->SetInputConnection(skeleton->GetOutputPort());
    // tell the mapper to use the color array for the lines
    skeletonMapper->ScalarVisibilityOn();
    skeletonMapper->SetScalarModeToUseCellFieldData();
    skeletonMapper->SelectColorArray("Colors");
    skeletonActor = vtkSmartPointer<vtkActor>::New();
    skeletonActor->SetMapper(skeletonMapper);
    skeletonActor->GetProperty()->SetLineWidth(2.5);
    skeletonActor->GetProperty()->SetDiffuseColor(0.8, 0.8, 0.8);
    skeletonActor->GetProperty()->SetOpacity(1.0);

    // create an data structure, mapper and actor
    // for the triangulated computational model (mesh)
    selectModelGraphics = SelectVLM;
    mesh = vtkSmartPointer<vtkPolyData>::New();
    meshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    meshMapper->SetInputData(mesh);
    meshActor = vtkSmartPointer<vtkActor>::New();
    meshActor->SetMapper(meshMapper);
    meshActor->GetProperty()->SetLineWidth(1.5);
    meshActor->GetProperty()->SetDiffuseColor(0.8, 0.8, 1.0);
    meshActor->GetProperty()->SetOpacity(1.0);
    scalMuMin=-1.0; scalMuMax=1.0;
    scalSigMin=-1.0; scalSigMax=1.0;
    scalCpMin=-1.0; scalCpMax=1.0;


    // creata a color legend for the mesh display
    colorScale = vtkSmartPointer<vtkScalarBarActor>::New();
    colorScale->SetTitle("scale");
    colorScale->SetLookupTable(meshMapper->GetLookupTable());
    colorScale->DrawAnnotationsOn();
    colorScale->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    // rechts oben
    // colorScale->GetPositionCoordinate()->SetValue(0.88, 0.68);
    // links unten
    colorScale->GetPositionCoordinate()->SetValue(0.02, 0.02);
    colorScale->SetWidth(0.10);
    colorScale->SetBarRatio(0.4);
    colorScale->SetHeight(0.3);
    colorScale->SetDrawFrame(0);

    // create an data structure, mapper and actor for the wake
    wake = vtkSmartPointer<vtkPolyData>::New();
    wakeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    wakeMapper->SetInputData(wake);
    wakeActor = vtkSmartPointer<vtkActor>::New();
    wakeActor->SetMapper(wakeMapper);
    wakeActor->GetProperty()->SetLineWidth(1.5);
    wakeActor->GetProperty()->SetDiffuseColor(1.0, 1.0, 0.2);
    wakeActor->GetProperty()->SetOpacity(1.0);

    // create actors for the plots
    chartStripePlot = vtkSmartPointer<vtkChartXY>::New();
    chartStripePlot->SetTitle("pressure distribution");
    plotStripeActor = vtkSmartPointer<vtkContextActor>::New();
    plotStripeActor->GetScene()->AddItem(chartStripePlot);
    chartGammaPlot = vtkSmartPointer<vtkChartXY>::New();
    chartGammaPlot->SetTitle("circulation distribution");
    plotGammaActor = vtkSmartPointer<vtkContextActor>::New();
    plotGammaActor->GetScene()->AddItem(chartGammaPlot);
    chartClCdPlot = vtkSmartPointer<vtkChartXY>::New();
    chartClCdPlot->SetTitle("lift/drag distribution");
    plotClCdActor = vtkSmartPointer<vtkContextActor>::New();
    plotClCdActor->GetScene()->AddItem(chartClCdPlot);

    // create renderers to display the graphics
    // they are selectively added to the graphWidget in updateGraph()
    graphRenderer = vtkRenderer::New();
    graphRenderer->AddActor(skeletonActor);
    graphRenderer->AddActor(meshActor);
    graphRenderer->AddActor(colorScale);
    graphRenderer->AddActor(wakeActor);
    graphRenderer->SetBackground(0.10, 0.20, 0.40);
    plotRenderer = vtkRenderer::New();
    plotRenderer->AddActor(plotStripeActor);
    plotRenderer->AddActor(plotGammaActor);
    plotRenderer->AddActor(plotClCdActor);
    plotRenderer->SetBackground(1.0, 1.0, 1.0);

    camera = graphRenderer->GetActiveCamera();
    camera->SetPosition(100.0, -1000.0, 300.0);
    camera->SetFocalPoint(0.0, 0.0, 0.0);
    camera->SetViewUp(0.0, 0.3, 1.0);
    camera->SetClippingRange(100.0, 2000.0);

    // initialize global data structures
    Globals::MainTextDisplay = ui.TextDisplay;
    Globals::MainTextDisplay->append(
        QString("XWing Version 2.0\n")+
        QString("U.Lehnert 5/2011\n\n")+
        QString("This program is free software; you can redistribute it and/or\n")+
        QString("modify it under the terms of the GNU General Public License\n")+
        QString("as published by the Free Software Foundation; either version 2\n")+
        QString("of the License, or (at your option) any later version.\n\n")+
        QString("This program is distributed in the hope that it will be useful,\n")+
        QString("but WITHOUT ANY WARRANTY - without even the implied warranty of\n")+
        QString("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n")+
        QString("See the GNU General Public License for more details.\n\n"));

    Globals::MainStatusBar = ui.statusbar;
    Globals::lastDataDirectory = QDir::current().absolutePath();
    Globals::lastModelFile = "";

    // setup default model data
    model = new Model();
    Globals::GeometryModel = model;
    flowVLM = NULL;
    flowSPM = NULL;

    // connect user interface to model data
    updateGeometryTab();
    ui.DataTab->setCurrentIndex(0);
    updateGraph();
    ui.GraphicsTab->setCurrentIndex(0);

}

void MainWindow::FileNew()
{
  delete model;
  model = new Model();
  Globals::GeometryModel = model;
  delete flowVLM;
  flowVLM = NULL;
  delete flowSPM;
  flowSPM = NULL;
  Globals::lastModelFile = "";
  updateGeometryTab();
  updateGraph();
  ui.statusbar->showMessage(QString("new model created"));
  ui.selectAirfoil->clear();
  ui.selectAirfoil->addItem(QString("--Strak--"));
  this->setWindowTitle("XWing2");
}

void MainWindow::FileOpen()
{
  QString fileName = QFileDialog::getOpenFileName(this,
      tr("Open XWing2 model file"),
         Globals::lastDataDirectory,
         tr("XWing2 files (*.xw2)\n"
             "all files (*.*)"));
  Globals::MainTextDisplay->clear();
  Globals::MainTextDisplay->append(QString("XWing2 model read : "+fileName+"\n"));
  delete model;
  bool modelValid=FALSE;
  if (!fileName.isEmpty()) {
    Globals::lastDataDirectory = QFileInfo(fileName).absolutePath();
    QApplication::setOverrideCursor(Qt::WaitCursor);
    // read the file
    QDomDocument doc("XWingML");
    QFile file(fileName);
    if( !file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      ui.statusbar->showMessage(QString("cannot read ")+fileName);
    } else {
      ui.statusbar->showMessage(QString("reading model from file"));
      if( !doc.setContent( &file ) )
        ui.statusbar->showMessage(QString("error reading XML file"));
      file.close();
      // interpret data
      QDomElement root = doc.documentElement();
      if( root.tagName() != "xwing2" )
      {
        ui.statusbar->showMessage(QString("wrong XML file type"));
	modelValid=FALSE;
      } else {
        model = new Model(root);
	modelValid=TRUE;
	Globals::lastModelFile = fileName;
	Globals::GeometryModel = model;
      };
    };
    // reading finished
    QApplication::restoreOverrideCursor();
  };
  ui.statusbar->showMessage(QString("done reading."));
  if (modelValid)
  {
    this->setWindowTitle("XWing2 - " + model->getName());
    if (model->numberOfWings() < 1)
    {
      ui.statusbar->showMessage(QString("no valid model read - using default"));
      delete model;
      model = new Model();
      Globals::GeometryModel = model;
    };
  } else {
    ui.statusbar->showMessage(QString("no valid model read - using default"));
    model = new Model();
    Globals::GeometryModel = model;
    modelValid=TRUE;
  };
  // the airfoil database was read, so the ui has to be updated accordingly
  // first clean all old entries and reinsert --Strak--
  ui.selectAirfoil->blockSignals(TRUE);
  ui.selectAirfoil->clear();
  ui.selectAirfoil->addItem(QString("--Strak--"));
  // then insert all available airfoils
  int naf = model->numberOfAirfoils();
  for (int i=1; i<=naf; i++)
  {
    Airfoil *af = model->airfoilDatabaseAt(i);
    ui.selectAirfoil->addItem(QString(af->getName()));
  };
  // finally update the user interface
  updateGeometryTab();
  delete flowVLM;
  flowVLM = NULL;
  delete flowSPM;
  flowSPM = NULL;
  ui.GraphicsTab->setCurrentIndex(0);
  updateGraph();
  ui.selectAirfoil->blockSignals(FALSE);
}

void MainWindow::FileSave()
{
  QString fileName;
  if (Globals::lastModelFile.isEmpty())
  {
    fileName = QFileDialog::getSaveFileName(this,
	tr("Save XWing2 model file"),
	  Globals::lastDataDirectory,
	  tr("XWing2 files (*.xw2)\n"
	      "all files (*.*)"));
  } else {
    fileName = Globals::lastModelFile;
  };
  if (!fileName.isEmpty()) {
    Globals::MainTextDisplay->append(QString("saving XWing2 model : "+fileName+"\n"));
    Globals::lastDataDirectory = QFileInfo(fileName).absolutePath();
    Globals::lastModelFile = fileName;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    // create XML document
    QDomDocument doc("XWingML");
    QDomElement root = doc.createElement("xwing2");
    doc.appendChild(root);
    model->writeGeometryXML(doc, root);
    // write document to file
    QFile file(fileName);
    if( !file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      ui.statusbar->showMessage(QString("cannot write ")+fileName);
    } else {
      ui.statusbar->showMessage(QString("saving model"));
      QTextStream ts( &file );
      ts << doc.toString();
      file.close();
    };
    QApplication::restoreOverrideCursor();
  };
  ui.statusbar->showMessage(QString("done writing."));
}

void MainWindow::FileSaveAs()
{
  QString fileName = QFileDialog::getSaveFileName(this,
      tr("Save XWing2 model file"),
	Globals::lastDataDirectory,
	tr("XWing2 files (*.xw2);;all files (*.*)"));
  if (!fileName.isEmpty()) {
    Globals::MainTextDisplay->append(QString("saving XWing2 model : "+fileName+"\n"));
    Globals::lastDataDirectory = QFileInfo(fileName).absolutePath();
    Globals::lastModelFile = fileName;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    // create XML document
    QDomDocument doc("XWingML");
    QDomElement root = doc.createElement("xwing2");
    doc.appendChild(root);
    model->writeGeometryXML(doc, root);
    // write document to file
    QFile file(fileName);
    if( !file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      ui.statusbar->showMessage(QString("cannot write ")+fileName);
    } else {
      ui.statusbar->showMessage(QString("saving model"));
      QTextStream ts( &file );
      ts << doc.toString();
      file.close();
    };
    QApplication::restoreOverrideCursor();
  };
  ui.statusbar->showMessage(QString("done writing."));
}

void MainWindow::SaveTextDisplay()
{
  QString fileName = QFileDialog::getSaveFileName(this,
      tr("Save XWing2 output text to file"),
	Globals::lastDataDirectory,
	tr("text files (*.txt)\n"
	    "all files (*.*)"));
    if (!fileName.isEmpty())
    {
      QFile file(fileName);
      if( !file.open(QIODevice::WriteOnly | QIODevice::Text))
      {
	ui.statusbar->showMessage(QString("cannot write ")+fileName);
      } else {
	QTextStream ts( &file );
	ts << Globals::MainTextDisplay->toPlainText();
	file.close();
        Globals::MainTextDisplay->append(QString("text saved to "+fileName+"\n"));
      };
    } else {
      ui.statusbar->showMessage(QString("could not write text file"));
    }
}

void MainWindow::LoadAirfoil()
{
  QString fileName = QFileDialog::getOpenFileName(this,
      tr("Open airfoil file"),
      Globals::lastDataDirectory,
      tr("airfoil files (*.dat)\n"
         "all files (*.*)"));
  if (!fileName.isEmpty()) {
    Globals::lastDataDirectory = QFileInfo(fileName).absolutePath();
    Airfoil *af = new Airfoil(fileName);
    if (af->isValid())
    {
      model->addAirfoil(af);
      ui.selectAirfoil->addItem(QString(af->getName()));
    } else {
      Globals::MainTextDisplay->append(QString("failed reading airfoil"));
      delete af;
    };
  };
}

void MainWindow::DeleteAirfoil()
{
  // we use a standard dialog for choosing the airfoil
  int i, id;
  QStringList items;
  int naf = model->numberOfAirfoils();
  for (i=1; i<=naf; i++)
  {
    items.append(model->airfoilDatabaseAt(i)->getName());
  };
  bool ok;
  QString item = QInputDialog::getItem(this, "select airfoil", "", items, 0, false, &ok);
  if (ok && !item.isEmpty())
  {
    Globals::MainTextDisplay->append(QString("delete airfoil : ")+item);
    id=0;
    for (i=1; i<=naf; i++)
    {
      if (model->airfoilDatabaseAt(i)->getName() == item) id = i;
    }
    if ((id>=1) && (id<=naf))
    {
      // delete the airfoil
      model->deleteAirfoil(id);
      // recreate the selection list from the database
      // (this is safer, than just deleting the one entry
      // first clean all old entries and reinsert --Strak--
      ui.selectAirfoil->clear();
      ui.selectAirfoil->addItem(QString("--Strak--"));
      // then reinsert all available airfoils
      int naf = model->numberOfAirfoils();
      for (i=1; i<=naf; i++)
      {
        Airfoil *af = model->airfoilDatabaseAt(i);
        ui.selectAirfoil->addItem(QString(af->getName()));
      };
      updateGeometryTab();
      updateGraph();
    };
  };
}

void MainWindow::ExportXWing1()
{
  QString fileName = QFileDialog::getSaveFileName(this,
      tr("export XWing 1.0"),
	Globals::lastDataDirectory,
	tr("XWing 1.0 files (*.xw);;all files (*.*)"));
  if (!fileName.isEmpty()) {
    QFileInfo fileinfo = QFileInfo(fileName);
    QString path = fileinfo.absolutePath();
    Globals::lastDataDirectory = path;
    QString name = fileinfo.baseName();
    // force file suffix to ".stl"
    // fileName = path+QString("/")+name+QString(".stl");
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QFile file(fileName);
    if( !file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
      ui.statusbar->showMessage(QString("cannot write ")+fileName);
    } else {
      ui.statusbar->showMessage(QString("exporting to ")+fileName);
      QTextStream ts( &file );
      // here goes the actual export routine
      file.close();
    };
    QApplication::restoreOverrideCursor();
  };
}

void MainWindow::ExportSTL()
{
  if (flowSPM == 0)
  {
    ui.statusbar->showMessage(QString("no valid surface mesh model available - cannot export."));
    return;
  }
  if (flowSPM->isValid())
  {
    QString fileName = QFileDialog::getSaveFileName(this,
	tr("export mesh file"),
	  Globals::lastDataDirectory,
	  tr("STL files (*.stl);;all files (*.*)"));
    if (!fileName.isEmpty()) {
      QFileInfo fileinfo = QFileInfo(fileName);
      QString path = fileinfo.absolutePath();
      Globals::lastDataDirectory = path;
      QString name = fileinfo.baseName();
      // force file suffix to ".stl"
      // fileName = path+QString("/")+name+QString(".stl");
      QApplication::setOverrideCursor(Qt::WaitCursor);
      QFile file(fileName);
      if( !file.open(QIODevice::WriteOnly | QIODevice::Text))
      {
	ui.statusbar->showMessage(QString("cannot write ")+fileName);
      } else {
	ui.statusbar->showMessage(QString("exporting ")+fileName);
	QTextStream ts( &file );
	ts << QString("solid ")+name << endl;
	flowSPM->exportSTL(&ts);
	ts << QString("endsolid ")+name << endl;
	file.close();
      };
      QApplication::restoreOverrideCursor();
    };
  }
  else
  {
    ui.statusbar->showMessage(QString("no valid surface mesh model available - cannot export."));
  }
}

void MainWindow::HelpAbout()
{
  QMessageBox::about(this,"XWing2 - about",
    "XWing Version 2.0\n"
    "U.Lehnert 11/2013\n\n"
    "Computation of the potential flow\n"
    "around wings and complete airplanes\n"
    "based on a panel method.\n");
}

/******************************************************/
/****** private slots *********************************/
/******************************************************/

// Geometry Wing Tab
#include "mainwindow.geometry.cxx"

// Global Settings Tab
#include "mainwindow.global.cxx"

// Flight Condition Tab
#include "mainwindow.modeling.cxx"

// Graphics Tab
#include "mainwindow.graphics.cxx"
