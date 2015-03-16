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

//  ********************
//  vortex lattice model
//  ********************

void MainWindow::on_VLM_button_Create_pressed()
{
  if (flowVLM != NULL) delete flowVLM;
  flowVLM = new VortexLatticeModel(model);
  // set graphics display to show the triangulation
  if (flowVLM->isValid())
  {
    ui.GraphicsTab->setCurrentIndex(1);
  } else {
    ui.GraphicsTab->setCurrentIndex(0);
  }
  selectModelGraphics=SelectVLM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_VLM_button_Run_pressed()
{
  if (flowVLM != NULL)
  {
    if (flowVLM->isValid())
    {
      flowVLM->runModelAOA(model->getAOA());
      if (flowVLM->isSolved())
      {
	flowVLM->analyzeWake();
	flowVLM->analyzeCirculation();
	// preset the color range for visualization
	flowVLM->getMuRange(&scalMuMin, &scalMuMax);
      };
    } else {
      Globals::MainTextDisplay->append(QString("Cannot start computation - you have to create a flow model first."));
    }
  } else {
    Globals::MainTextDisplay->append(QString("Cannot start computation - you have to create a flow model first."));
  }
  selectModelGraphics=SelectVLM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_VLM_button_Wake_pressed()
{
  if (flowVLM != NULL)
  {
    if (flowVLM->isValid() && flowVLM->isSolved())
    {
      flowVLM->relaxWake();
    } else {
      Globals::MainTextDisplay->append(QString("Cannot recompute wake - you have to create and solve the flow model first."));
    }
  } else {
    Globals::MainTextDisplay->append(QString("Cannot recompute wake - you have to create and solve the flow model first."));
  }
  selectModelGraphics=SelectVLM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_VLM_button_PrintPaneling_pressed()
{
  if (flowVLM != NULL) flowVLM->printPaneling();
}

void MainWindow::on_VLM_button_PrintFlowSolution_pressed()
{
  if (flowVLM != NULL) flowVLM->printSolution();
}

void MainWindow::on_VLM_button_PrintCirculation_pressed()
{
  if (flowVLM != NULL) flowVLM->printCirculation();
}

//  *******************
//  surface panel model
//  *******************

void MainWindow::on_SPM_button_Create_pressed()
{
  if (flowSPM != NULL) delete flowSPM;
  flowSPM = new SourceDoubletModel(model);
  // set graphics display to show the triangulation
  if (flowSPM->isValid())
  {
    ui.GraphicsTab->setCurrentIndex(1);
  } else {
    ui.GraphicsTab->setCurrentIndex(0);
  }
  selectModelGraphics=SelectSPM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_SPM_button_Run_pressed()
{
  if (flowSPM != NULL)
  {
    if (flowSPM->isValid())
    {
      flowSPM->runModelAOA(model->getAOA());
      if (flowSPM->isSolved())
      {
	flowSPM->analyzeWake();
	flowSPM->analyzePressure();
	// preset the color range for visualization
	flowSPM->getSigmaRange(&scalSigMin, &scalSigMax);
	flowSPM->getMuRange(&scalMuMin, &scalMuMax);
	flowSPM->getCpRange(&scalCpMin, &scalCpMax);
      };
    } else {
      Globals::MainTextDisplay->append(QString("Cannot start computation - you have to create a flow model first."));
    }
  } else {
    Globals::MainTextDisplay->append(QString("Cannot start computation - you have to create a flow model first."));
  }
  selectModelGraphics=SelectSPM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_SPM_button_Wake_pressed()
{
  if (flowSPM != NULL)
  {
    if (flowSPM->isValid() && flowSPM->isSolved())
    {
      flowSPM->relaxWake();
    } else {
      Globals::MainTextDisplay->append(QString("Cannot recompute wake - you have to create and solve the flow model first."));
    }
  } else {
    Globals::MainTextDisplay->append(QString("Cannot recompute wake - you have to create and solve the flow model first."));
  }
  selectModelGraphics=SelectSPM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_SPM_button_PrintPaneling_pressed()
{
  if (flowSPM != NULL) flowSPM->printPaneling();
}

void MainWindow::on_SPM_button_PrintFlowSolution_pressed()
{
  if (flowSPM != NULL) flowSPM->printSolution();
}

void MainWindow::on_SPM_button_PrintCirculation_pressed()
{
  if (flowSPM != NULL) flowSPM->printCirculation();
}
