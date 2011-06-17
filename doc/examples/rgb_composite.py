#!/usr/bin/env python
#-*- coding:utf-8 -*-
import sys
from PyQt4 import QtGui

def main():
    app = QtGui.QApplication(sys.argv)
    win = RGBCompositeImageApp()
    win.show()
    sys.exit(app.exec_())

class RGBCompositeImageApp(QtGui.QWidget):
    def __init__(self):
        super(RGBCompositeImageApp, self).__init__()
        self.setupUI()
        
    def setupUI(self):        
        okButton = QtGui.QPushButton("OK")
        cancelButton = QtGui.QPushButton("Cancel")

        hbox = QtGui.QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(okButton)
        hbox.addWidget(cancelButton)

        vbox = QtGui.QVBoxLayout()
        vbox.addStretch(1)
        vbox.addLayout(hbox)
        
        self.setLayout(vbox)
        
        self.setWindowTitle('SunPy Demo - RGB Composite Image Creator')
        self.resize(1024, 768)

if __name__ == '__main__':
    sys.exit(main())