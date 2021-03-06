# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'zmwanalysiswidget.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from ImageViewZMW import *


class Ui_ZmwAnalysisWidget(object):
    def setupUi(self, ZmwAnalysisWidget):
        ZmwAnalysisWidget.setObjectName("ZmwAnalysisWidget")
        ZmwAnalysisWidget.resize(930, 739)
        self.centralwidget = QtWidgets.QWidget(ZmwAnalysisWidget)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.imagePlot = ImageViewZMW(self.centralwidget)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(8)
        font.setBold(False)
        font.setWeight(50)
        self.imagePlot.setFont(font)
        self.imagePlot.setObjectName("imagePlot")
        self.gridLayout.addWidget(self.imagePlot, 0, 0, 1, 1)
        ZmwAnalysisWidget.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(ZmwAnalysisWidget)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 930, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtWidgets.QMenu(self.menubar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuView = QtWidgets.QMenu(self.menubar)
        self.menuView.setObjectName("menuView")
        self.menuAnalyze = QtWidgets.QMenu(self.menubar)
        self.menuAnalyze.setObjectName("menuAnalyze")
        ZmwAnalysisWidget.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(ZmwAnalysisWidget)
        self.statusbar.setObjectName("statusbar")
        ZmwAnalysisWidget.setStatusBar(self.statusbar)
        self.actionLoad = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionLoad.setObjectName("actionLoad")
        self.actionCrop = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionCrop.setObjectName("actionCrop")
        self.actionFilter = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionFilter.setObjectName("actionFilter")
        self.actionZ_Project = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionZ_Project.setObjectName("actionZ_Project")
        self.actionView_Stack = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionView_Stack.setObjectName("actionView_Stack")
        self.actionBase_Call = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionBase_Call.setObjectName("actionBase_Call")
        self.actionCheck_Controls = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionCheck_Controls.setObjectName("actionCheck_Controls")
        self.actionSet_Threshold = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionSet_Threshold.setObjectName("actionSet_Threshold")
        self.actionClear_Settings = QtWidgets.QAction(ZmwAnalysisWidget)
        self.actionClear_Settings.setObjectName("actionClear_Settings")
        self.menuFile.addAction(self.actionLoad)
        self.menuEdit.addAction(self.actionCrop)
        self.menuEdit.addAction(self.actionFilter)
        self.menuEdit.addAction(self.actionClear_Settings)
        self.menuView.addAction(self.actionZ_Project)
        self.menuView.addAction(self.actionView_Stack)
        self.menuAnalyze.addAction(self.actionBase_Call)
        self.menuAnalyze.addAction(self.actionCheck_Controls)
        self.menuAnalyze.addAction(self.actionSet_Threshold)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.menubar.addAction(self.menuAnalyze.menuAction())

        self.retranslateUi(ZmwAnalysisWidget)
        QtCore.QMetaObject.connectSlotsByName(ZmwAnalysisWidget)

    def retranslateUi(self, ZmwAnalysisWidget):
        _translate = QtCore.QCoreApplication.translate
        ZmwAnalysisWidget.setWindowTitle(_translate("ZmwAnalysisWidget", "MainWindow"))
        self.menuFile.setTitle(_translate("ZmwAnalysisWidget", "File"))
        self.menuEdit.setTitle(_translate("ZmwAnalysisWidget", "Edit"))
        self.menuView.setTitle(_translate("ZmwAnalysisWidget", "View"))
        self.menuAnalyze.setTitle(_translate("ZmwAnalysisWidget", "Analyze"))
        self.actionLoad.setText(_translate("ZmwAnalysisWidget", "load"))
        self.actionCrop.setText(_translate("ZmwAnalysisWidget", "Crop"))
        self.actionFilter.setText(_translate("ZmwAnalysisWidget", "Filter"))
        self.actionZ_Project.setText(_translate("ZmwAnalysisWidget", "Z-Project"))
        self.actionView_Stack.setText(_translate("ZmwAnalysisWidget", "View Stack"))
        self.actionBase_Call.setText(_translate("ZmwAnalysisWidget", "Base Call"))
        self.actionCheck_Controls.setText(_translate("ZmwAnalysisWidget", "Check Controls"))
        self.actionSet_Threshold.setText(_translate("ZmwAnalysisWidget", "Set Threshold"))
        self.actionClear_Settings.setText(_translate("ZmwAnalysisWidget", "Clear Settings"))
