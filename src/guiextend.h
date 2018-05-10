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

#ifndef GUIEXTEND_H
#define GUIEXTEND_H

#include <QtGui>
#include <QEvent>
#include <QtWidgets/QGroupBox>

// Es wird eine eigene Klasse von der QGroupBox abgeleitet,
// um ein zusätzliches Signal FocusChange() zu generieren.
// Dazu müssen enterEvent() und leaveEvent() entsprechend überschrieben werden.

class FocusGroupBox : public QGroupBox
{
  Q_OBJECT

  public:
  
    FocusGroupBox(QWidget *parent = 0);
    
  signals :
    void FocusChange(int);
    
  public slots:
    
  protected:
    void enterEvent(QEvent *e);
    void leaveEvent(QEvent *e);
};

#endif