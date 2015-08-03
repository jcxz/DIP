#
# @author: Matus Fedorko
#

QT += core
QT += gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DIP

#CONFIG += console
#CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -g

# OpenCL
INCLUDEPATH += "$$(MY_LIB_PATH)\opencl_khronos_headers\opencl_1_2"
#INCLUDEPATH += "$$(AMDAPPSDKROOT)\include"
LIBS += -L"$$(AMDAPPSDKROOT)\lib\x86_64"
LIBS += -lOpenCL

# Boost
INCLUDEPATH += "$$(MY_LIB_PATH)\boost\boost_1_57_0"

# Boost.Compute
INCLUDEPATH += "$$(MY_LIB_PATH)\boost_compute\compute\include"

#Other includes
INCLUDEPATH += src
#INCLUDEPATH += src/utils
#INCLUDEPATH += src/gui
#INCLUDEPATH += src/core

#release: BUILD_DIR = build/release
#debug:   BUILD_DIR = build/debug

#OBJECTS_DIR = $$BUILD_DIR/obj
#MOC_DIR = $$BUILD_DIR/moc
#RCC_DIR = $$BUILD_DIR/qrc
#UI_DIR = $$BUILD_DIR/ui
#DESTDIR = ./

HEADERS += \
    src/utils/macros.h \
    src/utils/debug.h \
    src/utils/misc.h \
    src/utils/stats.h \
    src/utils/maths.h \
    src/utils/timer.h \
    src/utils/gllogger.h \
    src/utils/trackball.h \
    src/utils/geom.h \
    src/utils/ogl.h \
    src/utils/ocl.h \
    \
    src/core/base_renderer.h \
    src/core/instancing_renderer.h \
    src/core/point_sprite_renderer.h \
    src/core/curvature_flow_renderer.h \
    src/core/ray_cast_renderer.h \
    src/core/particle_system.h \
    src/core/particle_system_internal.h \
    src/core/debug_particle_system.h \
    src/core/fluid_particle_system.h \
    src/core/sph_ocl_common.h \
    src/core/sph_params.h \
    src/core/sph_naive.h \
    src/core/sph_optimized.h \
    src/core/sph_uniform_grid.h \
    \
    src/io/voxel_mesh.h \
    \
    src/test/benchmark.h \
    src/test/sph_test.h \
    src/test/sph_test_adapter.h \
    src/test/sph_naive_test_adapter.h \
    src/test/sph_optimized_test_adapter.h \
    src/test/sph_uniform_grid_test_adapter.h \
    \
    src/gui/display_widget.h

SOURCES += \
    src/main.cpp \
    \
    src/utils/timer.cpp \
    src/utils/gllogger.cpp \
    src/utils/trackball.cpp \
    src/utils/geom.cpp \
    src/utils/ogl.cpp \
    src/utils/ocl.cpp \
    \
    src/core/base_renderer.cpp \
    src/core/instancing_renderer.cpp \
    src/core/point_sprite_renderer.cpp \
    src/core/curvature_flow_renderer.cpp \
    src/core/ray_cast_renderer.cpp \
    src/core/particle_system.cpp \
    src/core/debug_particle_system.cpp \
    src/core/fluid_particle_system.cpp \
    src/core/sph_params.cpp \
    src/core/sph_naive.cpp \
    src/core/sph_optimized.cpp \
    src/core/sph_uniform_grid.cpp \
    \
    src/io/voxel_mesh.cpp \
    \
    src/test/benchmark.cpp \
    src/test/sph_test.cpp \
    src/test/sph_naive_test_adapter.cpp \
    src/test/sph_optimized_test_adapter.cpp \
    src/test/sph_uniform_grid_test_adapter.cpp \
    \
    src/gui/display_widget.cpp

FORMS +=

RESOURCES += \
    resources_clgl.qrc \
    resources_data.qrc

OTHER_FILES +=
