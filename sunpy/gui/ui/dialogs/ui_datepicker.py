# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/matt/programming/sunpy/sunpy/gui/ui/dialogs/datepicker.ui'
#
# Created: Tue Sep 13 14:11:26 2011
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_DatePickerDlg(object):
    def setupUi(self, DatePickerDlg):
        DatePickerDlg.setObjectName(_fromUtf8("DatePickerDlg"))
        DatePickerDlg.resize(322, 223)
        DatePickerDlg.setWindowTitle(QtGui.QApplication.translate("DatePickerDlg", "Select date", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(DatePickerDlg)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.calendarWidget = QtGui.QCalendarWidget(DatePickerDlg)
        self.calendarWidget.setObjectName(_fromUtf8("calendarWidget"))
        self.verticalLayout.addWidget(self.calendarWidget)
        self.buttonBox = QtGui.QDialogButtonBox(DatePickerDlg)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(DatePickerDlg)
        QtCore.QMetaObject.connectSlotsByName(DatePickerDlg)

    def retranslateUi(self, DatePickerDlg):
        pass

