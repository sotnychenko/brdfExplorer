TEMPLATE = app
CONFIG += qt4  #debug

isEmpty(prefix) {
        prefix = /local
}
isEmpty(prefix) {
	error("$prefix is undefined. Please pass prefix=<path> to qmake")
}

DEST = $$prefix
isEmpty(LIBDIR) {
	LIBDIR = $$system(pf-makevar lib 2>/dev/null)
}
isEmpty(LIBDIR) {
	LIBDIR = lib
}

TARGET = brdf
target.path = $$DEST/bin

HEADERS = *.h \
    projectToPCA.h \
    redsvd.hpp \
    redsvdFile.hpp \
    redsvdIncr.hpp \
    util.hpp \
    global.h
SOURCES = \
    BRDFAnalytic.cpp \
    BRDFBase.cpp \
    BRDFImageSlice.cpp \
    BRDFMeasuredAniso.cpp \
    BRDFMeasuredMERL.cpp \
    ColorVarWidget.cpp \
    FloatVarWidget.cpp \
    DGLFrameBuffer.cpp \
    DGLShader.cpp \
    IBLWidget.cpp \
    IBLWindow.cpp \
    ImageSliceWidget.cpp \
    ImageSliceWindow.cpp \
    LitSphereWindow.cpp \
    main.cpp \
    MainWindow.cpp \
    ParameterGroupWidget.cpp \
    ParameterWindow.cpp \
    SharedContextGLWidget.cpp \
    ShowingDockWidget.cpp \
    PlotCartesianWindow.cpp \
    PlotCartesianWidget.cpp \
    PlotPolarWidget.cpp \
    Plot3DWidget.cpp \
    LitSphereWidget.cpp \
    SimpleModel.cpp \
    Paths.cpp \
    ptex/PtexReader.cpp \
    ptex/PtexUtils.cpp \
    ptex/PtexCache.cpp \
    ptex/PtexHalf.cpp \
    projectToPCA.cpp \
    cnpy.cpp \
    redsvdFile.cpp \
    util.cpp


QT   += opengl
DEFINES += PTEX_STATIC NOMINMAX

macx {
	INCLUDEPATH += /usr/X11/include
}

brdfs.path = $$DEST/share/brdf/brdfs
brdfs.files = ../brdfs/*

data.path = $$DEST/share/brdf/data
data.files = ../data/*

images.path = $$DEST/share/brdf/images
images.files = ../images/*

probes.path = $$DEST/share/brdf/probes
probes.files = ../probes/*

shaderTemplates.path = $$DEST/share/brdf/shaderTemplates
shaderTemplates.files = ../shaderTemplates/*

pkgconfig.path = $$DEST/$$LIBDIR/pkgconfig
pkgconfig.files = brdf.pc

INSTALLS = target brdfs data images probes shaderTemplates pkgconfig




# Windows cross compile at disney
linux-mingw32-custom{
    DEFINES += GLEW_STATIC
    WINDOWS_BUILD=/jobs2/soft/users/aselle/windows-build
    INCLUDEPATH += $$WINDOWS_BUILD/glew-1.9.0/include/
    INCLUDEPATH += $$WINDOWS_BUILD/glut-3.7.6-bin/
    LIBS += -L$$WINDOWS_BUILD/glut-3.7.6-bin/
    LIBS += -L$$WINDOWS_BUILD/glew-1.9.0/lib/
    LIBS += -static-libgcc
    LIBS += -lglew32s

}

win32 {
DEFINES += GLEW_STATIC
 INCLUDEPATH += C:\sotnychenko\mpii\dependencies\include\
 INCLUDEPATH += C:\sotnychenko\mpii\dependencies\bin\
 INCLUDEPATH += C:\sotnychenko\mpii\dependencies\Eigen\include"
 LIBS +=C:\Qt\4.8.6\lib\libQtOpenGL4.a
 LIBS +=C:\Qt\4.8.6\lib\libQtGui4.a
 LIBS +=C:\Qt\4.8.6\lib\libQtCore4.a
 LIBS += C:\sotnychenko\mpii\dependencies\lib\glut32.lib
 LIBS += C:\sotnychenko\mpii\dependencies\lib\libglew32.dll.a



}
