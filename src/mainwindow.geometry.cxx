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

void MainWindow::updateGeometryTab()
{
  Vector v;
  // alle diese Zuweisungen erzeugen valueChanged() Signale,
  // die dann nochmal die entsprechenden Slots ausführen und diese
  // rufen rekursiv updateGeometryTab() auf
  // deswegen werden die Signale geblockt

  // elements of the Geometry Wing Tab : Wing
  ui.selectWing->setMaximum(model->numberOfWings());
  wing = model->getWing(ui.selectWing->value());
  ui.geometryWingName->blockSignals(TRUE);
  ui.geometryWingName->setText(wing->getName());
  ui.geometryWingName->blockSignals(FALSE);
  ui.geometryWingLeftClosure->blockSignals(TRUE);
  ui.geometryWingLeftClosure->setChecked(wing->getLeftClosure());
  ui.geometryWingLeftClosure->blockSignals(FALSE);
  ui.geometryWingRightClosure->blockSignals(TRUE);
  ui.geometryWingRightClosure->setChecked(wing->getRightClosure());
  ui.geometryWingRightClosure->blockSignals(FALSE);
  ui.labelWingSpan->setText(QString("span = %1 mm").arg(wing->wingSpan(),0,'f',2));
  ui.labelWingArea->setText(QString("total area = %1").arg(wing->wingArea(),0,'f',3) +
    QString::fromUtf8(" dm²"));
  ui.labelWingLambda->setText(QString("aspect ratio = %1").arg(wing->wingLambda(),0,'f',2));
  ui.labelWingLiftingArea->setText(QString("lifting area = %1").arg(wing->wingLiftingArea(),0,'f',3) +
    QString::fromUtf8(" dm²"));
  ui.labelWingXC->setText(QString("xC = %1 mm").arg(wing->wingXC(),0,'f',2));
  ui.labelWingLmu->setText(QString("ref. chord = %1 mm").arg(wing->wingLmu(),0,'f',2));
  ui.labelWingSideArea->setText(QString("side area = %1").arg(wing->wingSideArea(),0,'f',3) +
    QString::fromUtf8(" dm²"));
  ui.labelWingXS->setText(QString("xS = %1 mm").arg(wing->wingXS(),0,'f',2));
  // elements of the Geometry Wing Tab : Station
  ui.selectStation->setMaximum(wing->numberOfStations());
  station = wing->getStation(ui.selectStation->value());
  ui.selectAirfoil->blockSignals(TRUE);
  ui.selectAirfoil->setCurrentIndex(station->getAirfoilIndex());
  ui.selectAirfoil->blockSignals(FALSE);
  ui.doubleStationX->blockSignals(TRUE);
  ui.doubleStationX->setValue(station->getX());
  ui.doubleStationX->blockSignals(FALSE);
  ui.doubleStationY->blockSignals(TRUE);
  ui.doubleStationY->setValue(station->getY());
  ui.doubleStationY->blockSignals(FALSE);
  ui.doubleStationZ->blockSignals(TRUE);
  ui.doubleStationZ->setValue(station->getZ());
  ui.doubleStationZ->blockSignals(FALSE);
  ui.doubleStationChord->blockSignals(TRUE);
  ui.doubleStationChord->setValue(station->getChord());
  ui.doubleStationChord->blockSignals(FALSE);
  ui.doubleStationAlfa->blockSignals(TRUE);
  ui.doubleStationAlfa->setValue(station->getAlfa());
  ui.doubleStationAlfa->blockSignals(FALSE);
  ui.labelStationDihedral->setText(QString("dihedral angle = %1 deg").arg(station->getDihedral()*180.0/M_PI,0,'f',2));
  ui.labelStationStretchZ->setText(QString("stretch factor = %1").arg(station->getStretchZ(),0,'f',4));
  v = station->getChordVector();
  ui.labelStationChordVec->setText(QString("chord direction = (%1, %2, %3)").arg(v.x,0,'f',3).arg(v.y,0,'f',3).arg(v.z,0,'f',3));
  v = station->getUpVector();
  ui.labelStationUpVec->setText(QString("up direction = (%1, %2, %3)").arg(v.x,0,'f',3).arg(v.y,0,'f',3).arg(v.z,0,'f',3));
  // elements of the Geometry Wing Tab : Segment
  ui.selectSegment->setMaximum(wing->numberOfSegments());
  segment = wing->getSegment(ui.selectSegment->value());

  ui.segmentSpanwisePanels->blockSignals(TRUE);
  ui.segmentSpanwisePanels->setValue(segment->getSpanN());
  ui.segmentSpanwisePanels->blockSignals(FALSE);
  ui.segmentChordwisePanels->blockSignals(TRUE);
  ui.segmentChordwisePanels->setValue(segment->getChordN());
  ui.segmentChordwisePanels->blockSignals(FALSE);
  ui.segmentFlapPanels->blockSignals(TRUE);
  ui.segmentFlapPanels->setValue(segment->getFlapN());
  ui.segmentFlapPanels->blockSignals(FALSE);
  ui.segmentFlapGroup->blockSignals(TRUE);
  ui.segmentFlapGroup->setValue(segment->getFlapGroup());
  ui.segmentFlapGroup->blockSignals(FALSE);
  ui.doubleSegmentFlapDepthRight->blockSignals(TRUE);
  ui.doubleSegmentFlapDepthRight->setValue(segment->getFlapLeft());
  ui.doubleSegmentFlapDepthRight->blockSignals(FALSE);
  ui.doubleSegmentFlapDepthLeft->blockSignals(TRUE);
  ui.doubleSegmentFlapDepthLeft->setValue(segment->getFlapRight());
  ui.doubleSegmentFlapDepthLeft->blockSignals(FALSE);
  ui.labelSegmentSpan->setText(QString("span = %1 mm").arg(segment->segmentSpan(),0,'f',2));
  ui.labelSegmentArea->setText(QString("area = %1").arg(segment->segmentArea(),0,'f',3) +
    QString::fromUtf8(" dm²"));
  ui.labelSegmentXC->setText(QString("xC = %1 mm").arg(segment->segmentXC(),0,'f',2));
  ui.labelSegmentDihedral->setText(QString("dihedral angle = %1 deg").arg(segment->segmentDihedral()*180.0/M_PI,0,'f',2));

  // elements of the  Global Settings Tab
  ui.globalSettingsModelName->blockSignals(TRUE);
  ui.globalSettingsModelName->setText(model->getName());
  ui.globalSettingsModelName->blockSignals(FALSE);
  ui.globalSettingsModelAuthor->blockSignals(TRUE);
  ui.globalSettingsModelAuthor->setText(model->getAuthor());
  ui.globalSettingsModelAuthor->blockSignals(FALSE);
  ui.globalSettingsMass->blockSignals(TRUE);
  ui.globalSettingsMass->setValue(model->getMass());
  ui.globalSettingsMass->blockSignals(FALSE);
  ui.globalSettingsReferenceSpan->blockSignals(TRUE);
  ui.globalSettingsReferenceSpan->setValue(model->getRefSpan());
  ui.globalSettingsReferenceSpan->blockSignals(FALSE);
  ui.globalSettingsReferenceArea->blockSignals(TRUE);
  ui.globalSettingsReferenceArea->setValue(model->getRefArea());
  ui.globalSettingsReferenceArea->blockSignals(FALSE);
  ui.globalSettingsReferenceChord->blockSignals(TRUE);
  ui.globalSettingsReferenceChord->setValue(model->getRefChord());
  ui.globalSettingsReferenceChord->blockSignals(FALSE);
  ui.globalSettingsXinfinity->blockSignals(TRUE);
  ui.globalSettingsXinfinity->setValue(model->getXinfinity());
  ui.globalSettingsXinfinity->blockSignals(FALSE);
  ui.globalSettingsNumberWakePanels->blockSignals(TRUE);
  ui.globalSettingsNumberWakePanels->setValue(model->getWakePanelNumber());
  ui.globalSettingsNumberWakePanels->blockSignals(FALSE);
  ui.globalSettingsAoA->blockSignals(TRUE);
  ui.globalSettingsAoA->setValue(model->getAOA());
  ui.globalSettingsAoA->blockSignals(FALSE);
}

void MainWindow::on_selectWing_valueChanged()
{
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_geometryWingName_editingFinished()
{
  wing->setName(ui.geometryWingName->text());
  updateGeometryTab();
}

void MainWindow::on_addWingButton_pressed()
{
  model->addWing();
  model->update();
  updateGeometryTab();
  updateGraph();
  ui.statusbar->showMessage(QString("new wing added"));
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
}

void MainWindow::on_deleteWingButton_pressed()
{
  int n=ui.selectWing->value();
  if (model->numberOfWings()>1)
  {
    model->deleteWing(n);
    model->update();
    updateGeometryTab();
    updateGraph();
    ui.statusbar->showMessage(QString("wing No. %1 deleted").arg(n,3));
    if (flowVLM != NULL) delete flowVLM;
    flowVLM = NULL;
    if (flowSPM != NULL) delete flowSPM;
    flowSPM = NULL;
  } else
  {
    ui.statusbar->showMessage(QString("can't delete wing No. %1").arg(n,3));
  }
}

void MainWindow::on_geometryWingLeftClosure_stateChanged()
{
  if (ui.geometryWingLeftClosure->isChecked())
    wing->setLeftClosure(TRUE);
  else
    wing->setLeftClosure(FALSE);
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
}

void MainWindow::on_geometryWingRightClosure_stateChanged()
{
  if (ui.geometryWingRightClosure->isChecked())
    wing->setRightClosure(TRUE);
  else
    wing->setRightClosure(FALSE);
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
}

void MainWindow::on_selectStation_valueChanged()
{
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_addStationButton_pressed()
{
  // n is the logic number of the newly created station
  int n=ui.selectStation->value()+1;
  if (n<2) n=2;
  if (n>wing->numberOfStations()) n=wing->numberOfStations();
  wing->addStation(n);
  model->update();
  ui.selectStation->setValue(n);
  updateGeometryTab();
  updateGraph();
  ui.statusbar->showMessage(QString("new station No. %1 added").arg(n,3));
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
}

void MainWindow::on_deleteStationButton_pressed()
{
  int n=ui.selectStation->value();
  if (wing->numberOfStations()>2)
  {
    wing->deleteStation(n);
    model->update();
    updateGeometryTab();
    updateGraph();
    ui.statusbar->showMessage(QString("station No. %1 deleted").arg(n,3));
    if (flowVLM != NULL) delete flowVLM;
    flowVLM = NULL;
    if (flowSPM != NULL) delete flowSPM;
    flowSPM = NULL;
  } else
  {
    ui.statusbar->showMessage(QString("can't delete station No. %1").arg(n,3));
  }
}

void MainWindow::on_doubleStationX_editingFinished()
{
  station->setX(ui.doubleStationX->value());
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_doubleStationY_editingFinished()
{
  station->setY(ui.doubleStationY->value());
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_doubleStationZ_editingFinished()
{
  station->setZ(ui.doubleStationZ->value());
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_doubleStationChord_editingFinished()
{
  station->setChord(ui.doubleStationChord->value());
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_doubleStationAlfa_editingFinished()
{
  station->setAlfa(ui.doubleStationAlfa->value());
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_selectAirfoil_currentIndexChanged(int index)
{
  station->setAirfoilIndex(index);
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_selectSegment_valueChanged()
{
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_segmentSpanwisePanels_valueChanged(int value)
{
  segment->setSpanN(value);
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_segmentChordwisePanels_valueChanged(int value)
{
  segment->setChordN(value);
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_segmentChordwisePanelsAllButton_pressed()
{
  int value = segment->getChordN();
  int nseg = wing->numberOfSegments();
  // n is the logical number running from 1 to NumberOfSegments
  for (int n=1; n<=nseg; n++)
  {
    GeometrySegment *seg = wing->getSegment(n);
    seg->setChordN(value);
  }
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_segmentFlapPanels_valueChanged(int value)
{
  segment->setFlapN(value);
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_segmentFlapPanelsAllButton_pressed()
{
  int value = segment->getFlapN();
  int nseg = wing->numberOfSegments();
  // n is the logical number running from 1 to NumberOfSegments
  for (int n=1; n<=nseg; n++)
  {
    GeometrySegment *seg = wing->getSegment(n);
    seg->setFlapN(value);
  }
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_segmentFlapGroup_valueChanged(int value)
{
  segment->setFlapGroup(value);
  model->update();
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_doubleSegmentFlapDepthLeft_editingFinished()
{
  segment->setFlapLeft(ui.doubleSegmentFlapDepthLeft->value());
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_doubleSegmentFlapDepthRight_editingFinished()
{
  segment->setFlapRight(ui.doubleSegmentFlapDepthRight->value());
  model->update();
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = NULL;
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = NULL;
  updateGeometryTab();
  updateGraph();
}

void MainWindow::on_groupBoxGeometryWing_FocusChange(int value)
{
  if (value==1) model->setHighlighting(WingHighlight);
  ui.GraphicsTab->setCurrentIndex(0);
  updateGraph();
}

void MainWindow::on_groupBoxGeometrySegment_FocusChange(int value)
{
  if (value==1) model->setHighlighting(SegmentHighlight);
  ui.GraphicsTab->setCurrentIndex(0);
  updateGraph();
}

void MainWindow::on_groupBoxGeometryStation_FocusChange(int value)
{
  if (value==1) model->setHighlighting(StationHighlight);
  ui.GraphicsTab->setCurrentIndex(0);
  updateGraph();
}

void MainWindow::on_DataTab_currentChanged(int value)
{
  if (value!=0)
  {
    model->setHighlighting(NoHighlight);
    updateGraph();
  };
}
