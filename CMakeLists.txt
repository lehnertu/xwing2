cmake_minimum_required (VERSION 2.6)
project (XWing2)

set(CMAKE_VERBOSE_MAKEFILE ON)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "CMAKE_BUILD_TYPE not given; setting to 'Release'.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# VTK library
find_package(VTK REQUIRED)
if (VTK_FOUND)
    message(STATUS "VTK includes ${VTK_INCLUDE_DIRS}")
    message(STATUS "VTK libraries ${VTK_LIBRARIES}")
endif (VTK_FOUND)
set(MY_VTK_LIBRARIES vtkChartsCore
	vtkCommonCore
	vtkCommonDataModel
	vtkCommonExecutionModel
	vtkFiltersCore
	vtkFiltersExtraction
	vtkFiltersStatistics
	vtkImagingCore
	vtkImagingHybrid
	vtkInteractionStyle
	vtkRenderingAnnotation
	vtkRenderingCore
	vtkRenderingContext2D
	vtkRenderingFreeType
	vtkRenderingOpenGL
	vtkRenderingFreeTypeOpenGL
	vtkRenderingVolumeOpenGL
	vtkViewsContext2D
	vtkGUISupportQt
	vtkIOImage
	vtksys)

# OpenBLAS library (linear algebra)
# the library uses cmake build system, so there exists a CMakeLists.txt in there
set(OpenBLAS_INCLUDE_DIRS /usr/include/openblas)
# if (OpenBLAS_FOUND)
#     message(STATUS "OpenBLAS includes ${OpenBLAS_INCLUDE_DIRS}")
#     message(STATUS "OpenBLAS libraries ${OpenBLAS_LIBRARIES}")
# endif (OpenBLAS_FOUND)

# declaring files in our project
include_directories(${PROJECT_SOURCE_DIR}/src
                    ${VTK_INCLUDE_DIRS}
                    ${OpenBLAS_INCLUDE_DIRS})

set ( SOURCES
    src/airfoil.cxx
    src/flatpanel.cxx
    src/geometrymodel.cxx
    src/geometrysegment.cxx
    src/geometrystation.cxx
    src/geometrywing.cxx
    src/global.cxx
    src/guiextend.cxx
    src/main.cxx
    src/mainwindow.cxx
    src/sourcedoubletmodel.cxx
    src/spline.cxx
    src/streamline.cxx
    src/vector.cxx
    src/vortexlatticemodel.cxx
    src/wakestripe.cxx
)

set ( MOC_HEADERS
  src/mainwindow.h
)

set ( UIS
  src/mainwindow.ui
)

# set ( RESOURCES
#  main.qrc
# )

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
# set(CMAKE_AUTORCC ON)
# set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Next lines needed for building all Qt projects
find_package( Qt4 REQUIRED )
include( ${QT_USE_FILE} )
add_definitions( ${QT_DEFINITIONS} )
include_directories( ${CMAKE_BINARY_DIR} )

# Next, using precompiler, compiler and linker

# using Qt meta-system (precompiler)
QT4_ADD_RESOURCES( RES_SOURCES ${RESOURCES} )
QT4_WRAP_UI( UI_HEADERS ${UIS} )
QT4_WRAP_CPP( MOC_SRCS ${MOC_HEADERS} )

# compile
add_executable( PROJECT_NAME
            ${SOURCES}
            # ${MOC_SRCS}
            # ${RES_SOURCES}
            ${UI_HEADERS} )

# or use line below instead, if you using Windows ™ Operating System.
#add_executable( PROJECT_NAME WIN32 ${SOURCES} ${MOC_SRCS} ${RES_SOURCES} ${UI_HEADERS} )

# build it (link libraries)
target_link_libraries( PROJECT_NAME ${QT4_LIBRARIES} ${VTK_LIBRARIES} )
