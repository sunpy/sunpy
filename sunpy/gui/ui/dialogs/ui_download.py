# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/matt/programming/sunpy/sunpy/gui/ui/dialogs/download.ui'
#
# Created: Tue Sep 13 14:11:25 2011
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_DownloadDlg(object):
    def setupUi(self, DownloadDlg):
        DownloadDlg.setObjectName(_fromUtf8("DownloadDlg"))
        DownloadDlg.resize(545, 437)
        DownloadDlg.setWindowTitle(QtGui.QApplication.translate("DownloadDlg", "Retrieve Data", None, QtGui.QApplication.UnicodeUTF8))
        self.widget = QtGui.QWidget(DownloadDlg)
        self.widget.setGeometry(QtCore.QRect(20, 20, 331, 66))
        self.widget.setObjectName(_fromUtf8("widget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.widget)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.timeEdit = QtGui.QTimeEdit(self.widget)
        self.timeEdit.setObjectName(_fromUtf8("timeEdit"))
        self.horizontalLayout_2.addWidget(self.timeEdit)
        self.dateEdit = QtGui.QDateEdit(self.widget)
        self.dateEdit.setObjectName(_fromUtf8("dateEdit"))
        self.horizontalLayout_2.addWidget(self.dateEdit)
        self.pushButton = QtGui.QPushButton(self.widget)
        self.pushButton.setText(QtGui.QApplication.translate("DownloadDlg", "PushButton", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.horizontalLayout_2.addWidget(self.pushButton)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.timeEdit_2 = QtGui.QTimeEdit(self.widget)
        self.timeEdit_2.setObjectName(_fromUtf8("timeEdit_2"))
        self.horizontalLayout.addWidget(self.timeEdit_2)
        self.dateEdit_2 = QtGui.QDateEdit(self.widget)
        self.dateEdit_2.setObjectName(_fromUtf8("dateEdit_2"))
        self.horizontalLayout.addWidget(self.dateEdit_2)
        self.pushButton_2 = QtGui.QPushButton(self.widget)
        self.pushButton_2.setText(QtGui.QApplication.translate("DownloadDlg", "...", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_2.setObjectName(_fromUtf8("pushButton_2"))
        self.horizontalLayout.addWidget(self.pushButton_2)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(DownloadDlg)
        QtCore.QMetaObject.connectSlotsByName(DownloadDlg)

    def retranslateUi(self, DownloadDlg):
        pass

