QT += core
QT -= gui

TARGET = project
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    cnpy.cpp \
    redsvdFile.cpp \
    util.cpp

HEADERS += \
    cnpy.h \
    redsvd.hpp \
    redsvdFile.hpp \
    util.hpp \
    redsvdIncr.hpp

INCLUDEPATH += "C:\Program Files (x86)\Eigen\include"
