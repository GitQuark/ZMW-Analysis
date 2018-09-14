# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'filterdialog.ui'
#
# Created: Mon Nov 14 18:58:08 2016
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(400, 249)
        self.filterBox = QtGui.QComboBox(Dialog)
        self.filterBox.setGeometry(QtCore.QRect(81, 10, 97, 26))
        self.filterBox.setObjectName(_fromUtf8("filterBox"))
        self.filterBox.addItem(_fromUtf8(""))
        self.filterBox.addItem(_fromUtf8(""))
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(12, 12, 64, 16))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(12, 42, 87, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(12, 94, 87, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(12, 146, 86, 16))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.xRangeEntry = QtGui.QLineEdit(Dialog)
        self.xRangeEntry.setGeometry(QtCore.QRect(84, 63, 127, 21))
        self.xRangeEntry.setAlignment(QtCore.Qt.AlignCenter)
        self.xRangeEntry.setObjectName(_fromUtf8("xRangeEntry"))
        self.yRangeEntry = QtGui.QLineEdit(Dialog)
        self.yRangeEntry.setGeometry(QtCore.QRect(84, 115, 127, 21))
        self.yRangeEntry.setAlignment(QtCore.Qt.AlignCenter)
        self.yRangeEntry.setObjectName(_fromUtf8("yRangeEntry"))
        self.zRangeEntry = QtGui.QLineEdit(Dialog)
        self.zRangeEntry.setGeometry(QtCore.QRect(84, 167, 127, 21))
        self.zRangeEntry.setAlignment(QtCore.Qt.AlignCenter)
        self.zRangeEntry.setObjectName(_fromUtf8("zRangeEntry"))
        self.okButton = QtGui.QPushButton(Dialog)
        self.okButton.setGeometry(QtCore.QRect(280, 210, 110, 32))
        self.okButton.setObjectName(_fromUtf8("okButton"))
        self.cancelButton = QtGui.QPushButton(Dialog)
        self.cancelButton.setGeometry(QtCore.QRect(160, 210, 110, 32))
        self.cancelButton.setObjectName(_fromUtf8("cancelButton"))

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.filterBox.setItemText(0, _translate("Dialog", "Median", None))
        self.filterBox.setItemText(1, _translate("Dialog", "Mean", None))
        self.label.setText(_translate("Dialog", "Filter Type:", None))
        self.label_2.setText(_translate("Dialog", "X Filter Range:", None))
        self.label_3.setText(_translate("Dialog", "Y Filter Range:", None))
        self.label_4.setText(_translate("Dialog", "Z Filter Range:", None))
        self.xRangeEntry.setText(_translate("Dialog", "3", None))
        self.yRangeEntry.setText(_translate("Dialog", "3", None))
        self.zRangeEntry.setText(_translate("Dialog", "3", None))
        self.okButton.setText(_translate("Dialog", "OK", None))
        self.cancelButton.setText(_translate("Dialog", "Cancel", None))

