INCLUDEPATH += /c/PFM/compile/include
LIBS += -L /c/PFM/compile/lib -lbag -lchrtr2 -lnvutility -lbag -lbeecrypt -lhdf5 -lgdal -lxml2 -lpoppler -lwsock32 -liconv
DEFINES += WIN32 NVWIN3X
CONFIG += console
QMAKE_LFLAGS += 
######################################################################
# Automatically generated by qmake (2.01a) Fri Dec 11 16:04:03 2020
######################################################################

TEMPLATE = app
TARGET = bag_residual
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += version.h
SOURCES += main.c
