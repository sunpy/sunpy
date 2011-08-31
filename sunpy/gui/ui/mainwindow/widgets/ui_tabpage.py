# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/matt/programming/sunpy/sunpy/gui/ui/mainwindow/widgets/tabpage.ui'
#
# Created: Wed Aug 31 18:27:08 2011
#      by: PyQt4 UI code generator 4.8.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_TabPage(object):
    def setupUi(self, TabPage):
        TabPage.setObjectName(_fromUtf8("TabPage"))
        TabPage.setWindowModality(QtCore.Qt.NonModal)
        TabPage.resize(540, 401)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(TabPage.sizePolicy().hasHeightForWidth())
        TabPage.setSizePolicy(sizePolicy)
        self.horizontalLayout_2 = QtGui.QHBoxLayout(TabPage)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.FigureCanvas = QtGui.QWidget(TabPage)
        self.FigureCanvas.setObjectName(_fromUtf8("FigureCanvas"))
        self.pushButton = QtGui.QPushButton(self.FigureCanvas)
        self.pushButton.setGeometry(QtCore.QRect(110, 110, 97, 27))
        self.pushButton.setCheckable(True)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.verticalLayout.addWidget(self.FigureCanvas)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.plotToolBar = QtGui.QWidget(TabPage)
        self.plotToolBar.setObjectName(_fromUtf8("plotToolBar"))
        self.horizontalLayout.addWidget(self.plotToolBar)
        self.cmComboBox = QtGui.QComboBox(TabPage)
        self.cmComboBox.setObjectName(_fromUtf8("cmComboBox"))
        self.horizontalLayout.addWidget(self.cmComboBox)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(TabPage)
        QtCore.QMetaObject.connectSlotsByName(TabPage)

    def retranslateUi(self, TabPage):
        TabPage.setWindowTitle(QtGui.QApplication.translate("TabPage", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("TabPage", "PushButton", None, QtGui.QApplication.UnicodeUTF8))

