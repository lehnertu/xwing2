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

#ifndef MAINWINDOWFORM_H
#define MAINWINDOWFORM_H

#include <QtGui>
#include <QtXml/QDomDocument>

// these defines are necessary for VTK only when building with qmake
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include "QVTKWidget.h"
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include "vtkChartXY.h"
#include "vtkContextActor.h"
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkScalarBarActor.h"
#include <vtkSmartPointer.h>

#include "ui_mainwindow.h"
#include "global.h"
#include "vector.h"
#include "airfoil.h"
#include "guiextend.h"
#include "geometrystation.h"
#include "geometrysegment.h"
#include "geometrywing.h"
#include "geometrymodel.h"
#include "sourcedoubletmodel.h"
#include "vortexlatticemodel.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    MainWindow(QMainWindow *parent = 0);

public slots:

    void FileNew();
    void FileOpen();
    void FileSave();
    void FileSaveAs();
    void SaveTextDisplay();
    void LoadAirfoil();
    void DeleteAirfoil();
    void ExportXWing1();
    void ExportSTL();
    void HelpAbout();

private slots:

    // slots mit Namen, die Signalen der enstprechenden UI-Widgets entsprechen,
    // werden automatisch beim Programmstart mit diesen verknüpft

    // Geometry Wing Tab
    void on_selectWing_valueChanged();
    void on_geometryWingName_editingFinished();
    void on_addWingButton_pressed();
    void on_deleteWingButton_pressed();
    void on_geometryWingLeftClosure_stateChanged();
    void on_geometryWingRightClosure_stateChanged();

    void on_selectStation_valueChanged();
    void on_addStationButton_pressed();
    void on_deleteStationButton_pressed();
    void on_doubleStationX_editingFinished();
    void on_doubleStationY_editingFinished();
    void on_doubleStationZ_editingFinished();
    void on_doubleStationChord_editingFinished();
    void on_doubleStationAlfa_editingFinished();
    void on_selectAirfoil_currentIndexChanged(int index);

    void on_selectSegment_valueChanged();
    void on_segmentSpanwisePanels_valueChanged(int value);
    void on_segmentChordwisePanels_valueChanged(int value);
    void on_segmentChordwisePanelsAllButton_pressed();
    void on_segmentFlapPanels_valueChanged(int value);
    void on_segmentFlapPanelsAllButton_pressed();
    void on_segmentFlapGroup_valueChanged(int value);
    void on_doubleSegmentFlapDepthLeft_editingFinished();
    void on_doubleSegmentFlapDepthRight_editingFinished();

    void on_DataTab_currentChanged(int value);
    void on_groupBoxGeometryStation_FocusChange(int value);
    void on_groupBoxGeometrySegment_FocusChange(int value);
    void on_groupBoxGeometryWing_FocusChange(int value);

    // Global Settings Tab
    void on_globalSettingsModelName_editingFinished();
    void on_globalSettingsModelAuthor_editingFinished();
    void on_globalSettingsMass_editingFinished();
    void on_globalSettingsReferenceSpan_editingFinished();
    void on_globalSettingsReferenceArea_editingFinished();
    void on_globalSettingsReferenceChord_editingFinished();
    void on_globalSettingsNumberWakePanels_valueChanged(int value);
    void on_globalSettingsXinfinity_editingFinished();
    void on_globalSettingsAoA_editingFinished();

    // vortex lattice model tab
    void on_VLM_button_Create_pressed();
    void on_VLM_button_Run_pressed();
    void on_VLM_button_Wake_pressed();
    void on_VLM_button_PrintPaneling_pressed();
    void on_VLM_button_PrintFlowSolution_pressed();
    void on_VLM_button_PrintCirculation_pressed();

    // surface panel model tab
    void on_SPM_button_Create_pressed();
    void on_SPM_button_Run_pressed();
    void on_SPM_button_Wake_pressed();
    void on_SPM_button_PrintPaneling_pressed();
    void on_SPM_button_PrintFlowSolution_pressed();
    void on_SPM_button_PrintCirculation_pressed();

    // Graphics Tab
    void on_GraphicsTab_currentChanged();
    void on_graphicsFlowShowWake_stateChanged();
    void on_graphicsFlowShowPanelLines_stateChanged();
    void on_graphicsFlowShowPanelNormals_stateChanged();
    void on_graphicsFlowShowFlowDirection_stateChanged();
    void on_graphicsFlowRenderWireFrame_toggled();
    void on_graphicsFlowRenderCP_toggled();
    void on_graphicsFlowRenderSource_toggled();
    void on_graphicsFlowRenderDoublet_toggled();
    void on_graphicsFlowSelectModel_pressed();
    void on_doubleGraphicsRenderScaleMin_editingFinished();
    void on_doubleGraphicsRenderScaleMax_editingFinished();
    void on_graphicsSectionSelectModel_pressed();
    void on_graphicsSectionSelectStripe_valueChanged(int value);
    void on_buttonGraphicsSectionPrintStripe_pressed();
    void on_graphGammaShowVLM_toggled();
    void on_graphGammaShowSPM_toggled();
    void on_graphLiftDragShowVLM_toggled();
    void on_graphLiftDragShowSPM_toggled();

private:

    Ui::MainWindow ui;

    QWidget *graphWindow;
    QVTKWidget *graphWidget;

    Model *model;
    GeometryWing *wing;
    GeometryStation *station;
    GeometrySegment *segment;
    SourceDoubletModel *flowSPM;
    VortexLatticeModel *flowVLM;
    graphicsSelectionType selectModelGraphics;

    // rendering of the geometry graphics
    vtkSmartPointer<vtkAppendPolyData> skeleton;
    vtkSmartPointer<vtkPolyDataMapper> skeletonMapper;
    vtkSmartPointer<vtkActor> skeletonActor;

    // rendering of the paneled computational model
    vtkSmartPointer<vtkPolyData> mesh;
    vtkSmartPointer<vtkPolyDataMapper> meshMapper;
    vtkSmartPointer<vtkActor> meshActor;
    vtkSmartPointer<vtkScalarBarActor> colorScale;
    double scalMuMin, scalMuMax;
    double scalSigMin, scalSigMax;
    double scalCpMin, scalCpMax;

    // rendering of the wake
    vtkSmartPointer<vtkPolyData> wake;
    vtkSmartPointer<vtkPolyDataMapper> wakeMapper;
    vtkSmartPointer<vtkActor> wakeActor;

    // rendering of 2D plots
    vtkSmartPointer<vtkChartXY> chartStripePlot;
    vtkSmartPointer<vtkContextActor> plotStripeActor;
    vtkSmartPointer<vtkChartXY> chartGammaPlot;
    vtkSmartPointer<vtkContextActor> plotGammaActor;
    vtkSmartPointer<vtkChartXY> chartClCdPlot;
    vtkSmartPointer<vtkContextActor> plotClCdActor;

    // rendering display port
    vtkRenderer *graphRenderer;
    vtkRenderer *plotRenderer;
    vtkCamera *camera;

    // update all elements of the geometry UI
    // for instance after change of selection
    // this makes sure that *wing *station and *segment
    // point to the currently selected items
    // and all informations are updated/diplayed properly
    // this has to be executed after all changes of selections etc.
    void updateGeometryTab();

    // update all elements of the graphics control tab
    // for instance changes of the model selection
    // updateGraphicsTab() should usually be called before updateGraph()
    void updateGraphicsTab();

    // regenerate display objects for graphics display
    // depending on the selection made in the GraphicsTab
    // the choosen actors are activated/deactivated
    // and the display in rendered
    // updateGraphicsTab() should usually be called before updateGraph()
    void updateGraph();

};

#endif