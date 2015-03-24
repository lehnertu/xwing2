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

void MainWindow::updateGraph()
{
  // switch off all actors
  // just the needed actors will be switched visible in the following
  skeletonActor->VisibilityOff();
  meshActor->VisibilityOff();
  colorScale->VisibilityOff();
  wakeActor->VisibilityOff();
  plotGammaActor->VisibilityOff();
  plotClCdActor->VisibilityOff();
  graphWidget->GetRenderWindow()->RemoveRenderer(graphRenderer);
  graphWidget->GetRenderWindow()->RemoveRenderer(plotRenderer);

  switch(ui.GraphicsTab->currentIndex())
  {
    case 0 :
    // geometry display
      {
	// recompute the graphics data
	model->sourceVTK(skeleton,
	  ui.selectWing->value(),
	  ui.selectSegment->value(),
	  ui.selectStation->value()); // wing, segment station - logical number
	skeletonActor->VisibilityOn();
	graphWidget->GetRenderWindow()->AddRenderer(graphRenderer);
	break;
      };
    case 1 :
    // flow model display
      {
	if ( (selectModelGraphics==SelectVLM) && (flowVLM != NULL))
	{
	  if (ui.graphicsFlowRenderWireFrame->isChecked())
	  {
	    flowVLM->sourcePanelsVTK(mesh,
				      ui.graphicsFlowShowPanelLines->isChecked(),
				      ui.graphicsFlowShowPanelNormals->isChecked(),
				      ui.graphicsFlowShowFlowDirection->isChecked(),
				      WireFrame
				      );
	    // set the scaling range for color display
	    meshMapper->SetScalarRange(-1.0,1.0);
	    meshActor->VisibilityOn();
	    if (ui.graphicsFlowShowWake->isChecked())
	    {
	      flowVLM->sourceWakeVTK(wake);
	      wakeActor->VisibilityOn();
	    };
	  };
	  if (ui.graphicsFlowRenderDoublet->isChecked() && flowVLM->isSolved())
	  {
	    flowVLM->sourcePanelsVTK(mesh,
				      ui.graphicsFlowShowPanelLines->isChecked(),
				      ui.graphicsFlowShowPanelNormals->isChecked(),
				      ui.graphicsFlowShowFlowDirection->isChecked(),
				      RenderDoublet
				      );
	    meshActor->VisibilityOn();
	    if (ui.graphicsFlowShowWake->isChecked())
	    {
	      flowVLM->sourceWakeVTK(wake);
	      wakeActor->VisibilityOn();
	    };
	    // set the scaling range for color display
	    meshMapper->SetScalarRange(scalMuMin,scalMuMax);
	    // retrieve and display the full range of values
            colorScale->SetTitle("doublet");
	    colorScale->VisibilityOn();
	  };
	}
	if ( (selectModelGraphics==SelectSPM) && (flowSPM != NULL))
	{
	  if (ui.graphicsFlowRenderWireFrame->isChecked())
	  {
	    flowSPM->sourcePanelsVTK(mesh,
				      ui.graphicsFlowShowPanelLines->isChecked(),
				      ui.graphicsFlowShowPanelNormals->isChecked(),
				      ui.graphicsFlowShowFlowDirection->isChecked(),
				      WireFrame
				      );
	    // set the scaling range for color display
	    meshMapper->SetScalarRange(-1.0,1.0);
	    meshActor->VisibilityOn();
	    if (ui.graphicsFlowShowWake->isChecked())
	    {
	      flowSPM->sourceWakeVTK(wake);
	      wakeActor->VisibilityOn();
	    };
	  };
	  if (ui.graphicsFlowRenderCP->isChecked() && flowSPM->isSolved())
	  {
	    flowSPM->sourcePanelsVTK(mesh,
				      ui.graphicsFlowShowPanelLines->isChecked(),
				      ui.graphicsFlowShowPanelNormals->isChecked(),
				      ui.graphicsFlowShowFlowDirection->isChecked(),
				      RenderCP
				      );
	    meshActor->VisibilityOn();
	    if (ui.graphicsFlowShowWake->isChecked())
	    {
	      flowSPM->sourceWakeVTK(wake);
	      wakeActor->VisibilityOn();
	    };
	    // set the scaling range for color display
	    meshMapper->SetScalarRange(scalCpMin,scalCpMax);
	    colorScale->SetTitle("pressure");
	    colorScale->VisibilityOn();
	  };
	  if (ui.graphicsFlowRenderSource->isChecked() && flowSPM->isSolved())
	  {
	    flowSPM->sourcePanelsVTK(mesh,
				      ui.graphicsFlowShowPanelLines->isChecked(),
				      ui.graphicsFlowShowPanelNormals->isChecked(),
				      ui.graphicsFlowShowFlowDirection->isChecked(),
				      RenderSource
				      );
	    meshActor->VisibilityOn();
	    if (ui.graphicsFlowShowWake->isChecked())
	    {
	      flowSPM->sourceWakeVTK(wake);
	      wakeActor->VisibilityOn();
	    };
	    // set the scaling range for color display
	    meshMapper->SetScalarRange(scalSigMin,scalSigMax);
	    // retrieve and display the full range of values
            colorScale->SetTitle("source");
	    colorScale->VisibilityOn();
	  };
	  if (ui.graphicsFlowRenderDoublet->isChecked() && flowSPM->isSolved())
	  {
	    flowSPM->sourcePanelsVTK(mesh,
				      ui.graphicsFlowShowPanelLines->isChecked(),
				      ui.graphicsFlowShowPanelNormals->isChecked(),
				      ui.graphicsFlowShowFlowDirection->isChecked(),
				      RenderDoublet
				      );
	    meshActor->VisibilityOn();
	    if (ui.graphicsFlowShowWake->isChecked())
	    {
	      flowSPM->sourceWakeVTK(wake);
	      wakeActor->VisibilityOn();
	    };
	    // set the scaling range for color display
	    meshMapper->SetScalarRange(scalMuMin,scalMuMax);
	    // retrieve and display the full range of values
            colorScale->SetTitle("doublet");
	    colorScale->VisibilityOn();
	  };
	}
	graphWidget->GetRenderWindow()->AddRenderer(graphRenderer);
	break;
      };
    case 2 :
    // section display
    {
      break;
    };
    case 3 :
    // gamma 2D display
      {
	  // remove all previous plots
	chartGammaPlot->ClearPlots();
	if ( ui.graphGammaShowVLM->isChecked() && (flowVLM != NULL) )
	{
	  flowVLM->sourceGammaPlot(chartGammaPlot);
	};
	if ( ui.graphGammaShowSPM->isChecked() && (flowSPM != NULL) )
	{
	  flowSPM->sourceGammaPlot(chartGammaPlot);
	};
	plotGammaActor->VisibilityOn();
	graphWidget->GetRenderWindow()->AddRenderer(plotRenderer);
	break;
      };
    case 4 :
    // cl cd plot
      {
	chartClCdPlot->ClearPlots();
	if ( ui.graphLiftDragShowVLM->isChecked() && (flowVLM != NULL) )
	{
	  flowVLM->sourceClCdPlot(chartClCdPlot);
	};
	if ( ui.graphLiftDragShowSPM->isChecked() && (flowSPM != NULL) )
	{
	  flowSPM->sourceClCdPlot(chartClCdPlot);
	};
	plotClCdActor->VisibilityOn();
	graphWidget->GetRenderWindow()->AddRenderer(plotRenderer);
	break;
      };
  };
  // render graphics
  graphWidget->GetRenderWindow()->Render();
}

void MainWindow::updateGraphicsTab()
{
  double min, max;
  // update the panel content
  switch(ui.GraphicsTab->currentIndex())
  {
    case 0 :
    // geometry display
      {
	break;
      };
    case 1 :
    // flow model display
      {
	ui.labelGraphicsFlowRenderRangeMin->setVisible(FALSE);
	ui.labelGraphicsFlowRenderRangeMax->setVisible(FALSE);
	ui.labelGraphicsFlowRenderScaleMin->setVisible(FALSE);
	ui.labelGraphicsFlowRenderScaleMax->setVisible(FALSE);
	ui.graphicsFlowRenderRangeMin->setVisible(FALSE);
	ui.graphicsFlowRenderRangeMax->setVisible(FALSE);
	ui.doubleGraphicsRenderScaleMin->setVisible(FALSE);
	ui.doubleGraphicsRenderScaleMax->setVisible(FALSE);
	// vortex lattice model
	if (selectModelGraphics==SelectVLM)
	{
	  ui.graphicsFlowSelectModel->setText("VLM");
	  ui.graphicsFlowSelectModel->update();
	  ui.graphicsFlowRenderCP->blockSignals(TRUE);
	  ui.graphicsFlowRenderCP->setChecked(FALSE);
	  ui.graphicsFlowRenderCP->setCheckable(FALSE);
	  ui.graphicsFlowRenderCP->setEnabled(FALSE);
	  ui.graphicsFlowRenderCP->blockSignals(FALSE);
	  ui.graphicsFlowRenderSource->blockSignals(TRUE);
	  ui.graphicsFlowRenderSource->setChecked(FALSE);
	  ui.graphicsFlowRenderSource->setCheckable(FALSE);
	  ui.graphicsFlowRenderSource->setEnabled(FALSE);
	  ui.graphicsFlowRenderSource->blockSignals(FALSE);
	  if (flowVLM==NULL)
	  {
	    ui.graphicsFlowRenderWireFrame->blockSignals(TRUE);
	    ui.graphicsFlowRenderWireFrame->setChecked(TRUE);
	    ui.graphicsFlowRenderWireFrame->blockSignals(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(TRUE);
	    ui.graphicsFlowRenderDoublet->setCheckable(FALSE);
	    ui.graphicsFlowRenderDoublet->setChecked(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(FALSE);
	  };
	};
	if ( (selectModelGraphics==SelectVLM) && (flowVLM != NULL))
	{
	  if (flowVLM->isSolved())
	  {
	    // enable all display options
	    ui.graphicsFlowRenderWireFrame->setCheckable(TRUE);
	    ui.graphicsFlowRenderDoublet->setCheckable(TRUE);
	  }
	  else
	  {
	    // if we dont have a valid solution for the flow field
	    // we only can render the wire-frame model
	    ui.graphicsFlowRenderWireFrame->blockSignals(TRUE);
	    ui.graphicsFlowRenderWireFrame->setChecked(TRUE);
	    ui.graphicsFlowRenderWireFrame->blockSignals(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(TRUE);
	    ui.graphicsFlowRenderDoublet->setChecked(FALSE);
	    ui.graphicsFlowRenderDoublet->setCheckable(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(FALSE);
	  };
	  if (ui.graphicsFlowRenderDoublet->isChecked() && flowVLM->isSolved())
	  {
	    ui.labelGraphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMax->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setVisible(TRUE);
	    // set the input box for the scale range to the previously defined values
	    ui.doubleGraphicsRenderScaleMin->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setValue(scalMuMin);
	    ui.doubleGraphicsRenderScaleMin->blockSignals(FALSE);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setValue(scalMuMax);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(FALSE);
	    // retrieve and display the full range of values
	    flowVLM->getMuRange(&min, &max);
	    ui.graphicsFlowRenderRangeMin->setText(QString("%1").arg(min,0,'f',2));
	    ui.graphicsFlowRenderRangeMax->setText(QString("%1").arg(max,0,'f',2));
	  };
	};
	// surface panel model
	if (selectModelGraphics==SelectSPM)
	{
	  ui.graphicsFlowSelectModel->setText("SPM");
	  ui.graphicsFlowSelectModel->update();
	  ui.graphicsFlowRenderWireFrame->setEnabled(TRUE);
	  ui.graphicsFlowRenderCP->setEnabled(TRUE);
	  ui.graphicsFlowRenderSource->setEnabled(TRUE);
	  ui.graphicsFlowRenderDoublet->setEnabled(TRUE);
	  if (flowSPM==NULL)
	  {
	    ui.graphicsFlowRenderWireFrame->blockSignals(TRUE);
	    ui.graphicsFlowRenderWireFrame->setChecked(TRUE);
	    ui.graphicsFlowRenderWireFrame->blockSignals(FALSE);
	    ui.graphicsFlowRenderCP->blockSignals(TRUE);
	    ui.graphicsFlowRenderCP->setChecked(FALSE);
	    ui.graphicsFlowRenderCP->setCheckable(FALSE);
	    ui.graphicsFlowRenderCP->blockSignals(FALSE);
	    ui.graphicsFlowRenderSource->blockSignals(TRUE);
	    ui.graphicsFlowRenderSource->setChecked(FALSE);
	    ui.graphicsFlowRenderSource->setCheckable(FALSE);
	    ui.graphicsFlowRenderSource->blockSignals(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(TRUE);
	    ui.graphicsFlowRenderDoublet->setChecked(FALSE);
	    ui.graphicsFlowRenderDoublet->setCheckable(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(FALSE);
	  };
	};
	if ( (selectModelGraphics==SelectSPM) && (flowSPM != NULL))
	{
	  if (flowSPM->isSolved())
	  {
	    // enable all display options
	    ui.graphicsFlowRenderWireFrame->setCheckable(TRUE);
	    ui.graphicsFlowRenderCP->setCheckable(TRUE);
	    ui.graphicsFlowRenderSource->setCheckable(TRUE);
	    ui.graphicsFlowRenderDoublet->setCheckable(TRUE);
	  }
	  else
	  {
	    // if we dont have a valid solution for the flow field
	    // we only can render the wire-frame model
	    ui.graphicsFlowRenderWireFrame->blockSignals(TRUE);
	    ui.graphicsFlowRenderWireFrame->setChecked(TRUE);
	    ui.graphicsFlowRenderWireFrame->blockSignals(FALSE);
	    ui.graphicsFlowRenderCP->blockSignals(TRUE);
	    ui.graphicsFlowRenderCP->setChecked(FALSE);
	    ui.graphicsFlowRenderCP->setCheckable(FALSE);
	    ui.graphicsFlowRenderCP->blockSignals(FALSE);
	    ui.graphicsFlowRenderSource->blockSignals(TRUE);
	    ui.graphicsFlowRenderSource->setChecked(FALSE);
	    ui.graphicsFlowRenderSource->setCheckable(FALSE);
	    ui.graphicsFlowRenderSource->blockSignals(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(TRUE);
	    ui.graphicsFlowRenderDoublet->setChecked(FALSE);
	    ui.graphicsFlowRenderDoublet->setCheckable(FALSE);
	    ui.graphicsFlowRenderDoublet->blockSignals(FALSE);
	  };
	  if (ui.graphicsFlowRenderCP->isChecked() && flowSPM->isSolved())
	  {
	    ui.labelGraphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMax->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setVisible(TRUE);
	    // set the input box for the scale range to the previously defined values
	    ui.doubleGraphicsRenderScaleMin->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setValue(scalCpMin);
	    ui.doubleGraphicsRenderScaleMin->blockSignals(FALSE);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setValue(scalCpMax);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(FALSE);
	    // retrieve and display the full range of values
	    flowSPM->getCpRange(&min, &max);
	    ui.graphicsFlowRenderRangeMin->setText(QString("%1").arg(min,0,'f',2));
	    ui.graphicsFlowRenderRangeMax->setText(QString("%1").arg(max,0,'f',2));
	  };
	  if (ui.graphicsFlowRenderSource->isChecked() && flowSPM->isSolved())
	  {
	    ui.labelGraphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMax->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setVisible(TRUE);
	    // set the input box for the scale range to the previously defined values
	    ui.doubleGraphicsRenderScaleMin->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setValue(scalSigMin);
	    ui.doubleGraphicsRenderScaleMin->blockSignals(FALSE);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setValue(scalSigMax);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(FALSE);
	    // retrieve and display the full range of values
	    flowSPM->getSigmaRange(&min, &max);
	    ui.graphicsFlowRenderRangeMin->setText(QString("%1").arg(min,0,'f',2));
	    ui.graphicsFlowRenderRangeMax->setText(QString("%1").arg(max,0,'f',2));
	  };
	  if (ui.graphicsFlowRenderDoublet->isChecked() && flowSPM->isSolved())
	  {
	    ui.labelGraphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMin->setVisible(TRUE);
	    ui.labelGraphicsFlowRenderScaleMax->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMin->setVisible(TRUE);
	    ui.graphicsFlowRenderRangeMax->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setVisible(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setVisible(TRUE);
	    // set the input box for the scale range to the previously defined values
	    ui.doubleGraphicsRenderScaleMin->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMin->setValue(scalMuMin);
	    ui.doubleGraphicsRenderScaleMin->blockSignals(FALSE);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(TRUE);
	    ui.doubleGraphicsRenderScaleMax->setValue(scalMuMax);
	    ui.doubleGraphicsRenderScaleMax->blockSignals(FALSE);
	    // retrieve and display the full range of values
	    flowSPM->getMuRange(&min, &max);
	    ui.graphicsFlowRenderRangeMin->setText(QString("%1").arg(min,0,'f',2));
	    ui.graphicsFlowRenderRangeMax->setText(QString("%1").arg(max,0,'f',2));
	  };
	}
	break;
      };
    case 2 :
    // section display
      {
	if (selectModelGraphics==SelectVLM)
	{
	  ui.graphicsSectionSelectModel->setText("VLM");
          if (flowVLM != NULL)
          {
            if (flowVLM->isValid())
            {
              ui.graphicsSectionSelectStripe->setMaximum(0);
            }
            else
            {
              ui.graphicsSectionSelectStripe->setMaximum(0);
            };
          }
          else
          {
            ui.graphicsSectionSelectStripe->setMaximum(0);
          };
	}
	if (selectModelGraphics==SelectSPM)
	{
	  ui.graphicsSectionSelectModel->setText("SPM");
          if (flowSPM != NULL)
          {
            if (flowSPM->isValid())
            {
              ui.graphicsSectionSelectStripe->setMaximum(flowSPM->numberWakes()-1);
            }
            else
            {
              ui.graphicsSectionSelectStripe->setMaximum(0);
            };
          }
          else
          {
            ui.graphicsSectionSelectStripe->setMaximum(0);
          };
	}
	break;
      };
    case 3 :
    // gamma 2D display
      {
	break;
      };
    case 4 :
    // cl cd plot
      {
	break;
      };
  };
}

void MainWindow::on_GraphicsTab_currentChanged()
{
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphicsFlowShowWake_stateChanged()
{
  updateGraph();
}

void MainWindow::on_graphicsFlowShowPanelLines_stateChanged()
{
  updateGraph();
}

void MainWindow::on_graphicsFlowShowPanelNormals_stateChanged()
{
  updateGraph();
}

void MainWindow::on_graphicsFlowShowFlowDirection_stateChanged()
{
  updateGraph();
}

void MainWindow::on_graphicsFlowRenderWireFrame_toggled()
{
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphicsFlowRenderCP_toggled()
{
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphicsFlowRenderSource_toggled()
{
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphicsFlowRenderDoublet_toggled()
{
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphicsFlowSelectModel_pressed()
{
  if (selectModelGraphics==SelectVLM)
    selectModelGraphics=SelectSPM;
  else
    selectModelGraphics=SelectVLM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_doubleGraphicsRenderScaleMin_editingFinished()
{
  if (ui.graphicsFlowRenderCP->isChecked())
    scalCpMin=ui.doubleGraphicsRenderScaleMin->value();
  if (ui.graphicsFlowRenderSource->isChecked())
    scalSigMin=ui.doubleGraphicsRenderScaleMin->value();
  if (ui.graphicsFlowRenderDoublet->isChecked())
    scalMuMin=ui.doubleGraphicsRenderScaleMin->value();
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_doubleGraphicsRenderScaleMax_editingFinished()
{
  if (ui.graphicsFlowRenderCP->isChecked())
    scalCpMax=ui.doubleGraphicsRenderScaleMax->value();
  if (ui.graphicsFlowRenderSource->isChecked())
    scalSigMax=ui.doubleGraphicsRenderScaleMax->value();
  if (ui.graphicsFlowRenderDoublet->isChecked())
    scalMuMax=ui.doubleGraphicsRenderScaleMax->value();
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphicsSectionSelectModel_pressed()
{
  if (selectModelGraphics==SelectVLM)
    selectModelGraphics=SelectSPM;
  else
    selectModelGraphics=SelectVLM;
  updateGraphicsTab();
  updateGraph();
}

void MainWindow::on_graphGammaShowVLM_toggled()
{
  updateGraph();
}

void MainWindow::on_graphGammaShowSPM_toggled()
{
  updateGraph();
}

void MainWindow::on_graphLiftDragShowVLM_toggled()
{
  updateGraph();
}

void MainWindow::on_graphLiftDragShowSPM_toggled()
{
  updateGraph();
}

