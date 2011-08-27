# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/matt/programming/sunpy/sunpy/gui/ui/dialogs/cmselector.ui'
#
# Created: Sat Aug 27 23:54:41 2011
#      by: PyQt4 UI code generator 4.8.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_CMSelectorDlg(object):
    def setupUi(self, CMSelectorDlg):
        CMSelectorDlg.setObjectName(_fromUtf8("CMSelectorDlg"))
        CMSelectorDlg.resize(364, 177)
        CMSelectorDlg.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedKingdom))
        self.verticalLayout_2 = QtGui.QVBoxLayout(CMSelectorDlg)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label_2 = QtGui.QLabel(CMSelectorDlg)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.plotComboBox = QtGui.QComboBox(CMSelectorDlg)
        self.plotComboBox.setObjectName(_fromUtf8("plotComboBox"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.plotComboBox)
        self.label = QtGui.QLabel(CMSelectorDlg)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label)
        self.subplotComboBox = QtGui.QComboBox(CMSelectorDlg)
        self.subplotComboBox.setObjectName(_fromUtf8("subplotComboBox"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.subplotComboBox)
        self.label_3 = QtGui.QLabel(CMSelectorDlg)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_3)
        self.cmComboBox = QtGui.QComboBox(CMSelectorDlg)
        self.cmComboBox.setObjectName(_fromUtf8("cmComboBox"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.cmComboBox)
        self.verticalLayout.addLayout(self.formLayout)
        self.line = QtGui.QFrame(CMSelectorDlg)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.verticalLayout.addWidget(self.line)
        self.buttonBox = QtGui.QDialogButtonBox(CMSelectorDlg)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.label_2.setBuddy(self.plotComboBox)
        self.label.setBuddy(self.subplotComboBox)
        self.label_3.setBuddy(self.cmComboBox)

        self.retranslateUi(CMSelectorDlg)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), CMSelectorDlg.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), CMSelectorDlg.reject)
        QtCore.QMetaObject.connectSlotsByName(CMSelectorDlg)

    def retranslateUi(self, CMSelectorDlg):
        CMSelectorDlg.setWindowTitle(QtGui.QApplication.translate("CMSelectorDlg", "Color Map Selector", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("CMSelectorDlg", "&Plot: ", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("CMSelectorDlg", "&Subplot: ", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("CMSelectorDlg", "&Color map:", None, QtGui.QApplication.UnicodeUTF8))

