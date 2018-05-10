######################################################################
# Automatically generated by qmake (2.01a) Mo. Aug 19 18:40:46 2013
# manualy modified - do not regenerate !!
######################################################################

TEMPLATE = app
TARGET = XWing2


DEPENDPATH += . ./src
INCLUDEPATH += . ./src
SOURCES_DIR = ./src

# OPEN_BLAS_INCLUDE_PATH = /usr/include
OPEN_BLAS_INCLUDE_PATH = /home/ulf/Programming/OpenBLAS
OPEN_BLAS_LIB_PATH = /home/ulf/Programming/OpenBLAS/lib

VTK_INCLUDE_PATH = /usr/include/vtk-6.1
VTK_LIB_PATH = /usr/lib/x86_64-linux-gnu

INCLUDEPATH += $$OPEN_BLAS_INCLUDE_PATH
INCLUDEPATH += $$VTK_INCLUDE_PATH

MOC_DIR = ./obj
UI_SOURCES_DIR = ./src
UI_HEADERS_DIR = ./obj
OBJECTS_DIR = ./obj

# Input
HEADERS += src/airfoil.h \
           src/flatpanel.h \
           src/geometrymodel.h \
           src/geometrysegment.h \
           src/geometrystation.h \
           src/geometrywing.h \
           src/global.h \
           src/guiextend.h \
           src/mainwindow.h \
           src/sourcedoubletmodel.h \
           src/spline.h \
           src/streamline.h \
           src/vector.h \
           src/vortexlatticemodel.h \
           src/wakestripe.h

FORMS += src/mainwindow.ui

SOURCES += src/airfoil.cxx \
           src/flatpanel.cxx \
           src/geometrymodel.cxx \
           src/geometrysegment.cxx \
           src/geometrystation.cxx \
           src/geometrywing.cxx \
           src/global.cxx \
           src/guiextend.cxx \
           src/main.cxx \
           src/mainwindow.cxx \
           src/mainwindow.geometry.cxx \
           src/mainwindow.global.cxx \
           src/mainwindow.graphics.cxx \
           src/mainwindow.modeling.cxx \
           src/sourcedoubletmodel.cxx \
           src/spline.cxx \
           src/streamline.cxx \
           src/vector.cxx \
           src/vortexlatticemodel.cxx \
           src/wakestripe.cxx

LIBS += -L$$VTK_LIB_PATH \
	-lvtkChartsCore-6.1 \
	-lvtkCommonCore-6.1 \
	-lvtkCommonDataModel-6.1 \
	-lvtkCommonExecutionModel-6.1 \
	-lvtkFiltersCore-6.1 \
	-lvtkFiltersExtraction-6.1 \
	-lvtkFiltersStatistics-6.1 \
	-lvtkImagingCore-6.1 \
	-lvtkImagingHybrid-6.1 \
	-lvtkInteractionStyle-6.1 \
	-lvtkRenderingAnnotation-6.1 \
	-lvtkRenderingCore-6.1 \
	-lvtkRenderingContext2D-6.1 \
	-lvtkRenderingFreeType-6.1 \
	-lvtkRenderingOpenGL-6.1 \
	-lvtkRenderingFreeTypeOpenGL-6.1 \
	-lvtkRenderingVolumeOpenGL-6.1 \
	-lvtkViewsContext2D-6.1 \
	-lvtkGUISupportQt-6.1 \
	-lvtkIOImage-6.1 \
	-lvtksys-6.1 \
	$$OPEN_BLAS_LIB_PATH/libopenblas.a \
	/usr/lib/liblapacke.a \
	-lgomp

QMAKE_CXXFLAGS += -Wno-deprecated \
		  -fopenmp

QMAKE_LFLAGS += -Wl,-rpath,$$VTK_LIB_PATH

CONFIG += debug
# CONFIG += console

QT += xml
