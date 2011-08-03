#!/usr/bin/env python
#-*- coding:utf-8 -*-
import sys
from PyQt4 import QtGui
import Ui_RGBComposite

def main():
    app = QtGui.QApplication(sys.argv)
    win = RGBCompositeImageApp()
    win.show()
    sys.exit(app.exec_())

#        super(RGBCompositeImageApp, self).__init__()
class RGBCompositeImageApp(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = Ui_RGBComposite.Ui_RGBComposite()
        self.ui.setupUi(self)

if __name__ == '__main__':
    sys.exit(main())