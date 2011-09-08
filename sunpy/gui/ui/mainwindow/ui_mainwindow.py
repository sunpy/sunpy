# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/matt/programming/sunpy/sunpy/gui/ui/mainwindow/mainwindow.ui'
#
# Created: Thu Sep  8 17:25:24 2011
#      by: PyQt4 UI code generator 4.8.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(662, 615)
        MainWindow.setCursor(QtCore.Qt.ArrowCursor)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/main.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        MainWindow.setWindowOpacity(1.0)
        MainWindow.setLayoutDirection(QtCore.Qt.LeftToRight)
        MainWindow.setToolButtonStyle(QtCore.Qt.ToolButtonIconOnly)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setTabPosition(QtGui.QTabWidget.South)
        self.tabWidget.setElideMode(QtCore.Qt.ElideNone)
        self.tabWidget.setTabsClosable(True)
        self.tabWidget.setMovable(True)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.horizontalLayout.addWidget(self.tabWidget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 662, 23))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.fileToolBar = QtGui.QToolBar(MainWindow)
        self.fileToolBar.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.fileToolBar.setIconSize(QtCore.QSize(22, 22))
        self.fileToolBar.setObjectName(_fromUtf8("fileToolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.fileToolBar)
        self.dockWidget = QtGui.QDockWidget(MainWindow)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dockWidget.sizePolicy().hasHeightForWidth())
        self.dockWidget.setSizePolicy(sizePolicy)
        self.dockWidget.setFloating(False)
        self.dockWidget.setFeatures(QtGui.QDockWidget.DockWidgetFloatable|QtGui.QDockWidget.DockWidgetMovable)
        self.dockWidget.setObjectName(_fromUtf8("dockWidget"))
        self.dockWidgetContents = QtGui.QWidget()
        self.dockWidgetContents.setObjectName(_fromUtf8("dockWidgetContents"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.dockWidgetContents)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.widget = ColorOptions(self.dockWidgetContents)
        self.widget.setAutoFillBackground(False)
        self.widget.setObjectName(_fromUtf8("widget"))
        self.horizontalLayout_2.addWidget(self.widget)
        self.dockWidget.setWidget(self.dockWidgetContents)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(2), self.dockWidget)
        self.actionOpen_file = QtGui.QAction(MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/open_plot.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionOpen_file.setIcon(icon1)
        self.actionOpen_file.setObjectName(_fromUtf8("actionOpen_file"))
        self.actionExit_PlotMan = QtGui.QAction(MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/exit.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionExit_PlotMan.setIcon(icon2)
        self.actionExit_PlotMan.setObjectName(_fromUtf8("actionExit_PlotMan"))
        self.menuFile.addAction(self.actionOpen_file)
        self.menuFile.addAction(self.actionExit_PlotMan)
        self.menubar.addAction(self.menuFile.menuAction())
        self.fileToolBar.addAction(self.actionOpen_file)
        self.fileToolBar.addAction(self.actionExit_PlotMan)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(-1)
        QtCore.QObject.connect(self.actionExit_PlotMan, QtCore.SIGNAL(_fromUtf8("triggered()")), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "SunPy PlotMan", None, QtGui.QApplication.UnicodeUTF8))
        self.menuFile.setTitle(QtGui.QApplication.translate("MainWindow", "&File", None, QtGui.QApplication.UnicodeUTF8))
        self.fileToolBar.setWindowTitle(QtGui.QApplication.translate("MainWindow", "File", None, QtGui.QApplication.UnicodeUTF8))
        self.dockWidget.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Color options", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen_file.setText(QtGui.QApplication.translate("MainWindow", "Open file...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen_file.setToolTip(QtGui.QApplication.translate("MainWindow", "Open a FITS file for plotting...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen_file.setStatusTip(QtGui.QApplication.translate("MainWindow", "Open a FITS file for plotting...", None, QtGui.QApplication.UnicodeUTF8))
        self.actionOpen_file.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit_PlotMan.setText(QtGui.QApplication.translate("MainWindow", "Exit PlotMan", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit_PlotMan.setToolTip(QtGui.QApplication.translate("MainWindow", "Close the application.", None, QtGui.QApplication.UnicodeUTF8))
        self.actionExit_PlotMan.setStatusTip(QtGui.QApplication.translate("MainWindow", "Close the application.", None, QtGui.QApplication.UnicodeUTF8))

from widgets.coloroptions import ColorOptions
from resources import qrc_resources
