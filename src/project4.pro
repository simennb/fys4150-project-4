TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS -= -std=c++0x

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    ising.cpp \
    functions.cpp

HEADERS += \
    ising.h \
    functions.h

