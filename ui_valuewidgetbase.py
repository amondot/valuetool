# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_valuewidgetbase.ui'
#
# Created: Sun Oct 28 17:48:41 2012
#      by: PyQt4 UI code generator 4.8.5
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_ValueWidgetBase(object):
    def setupUi(self, ValueWidgetBase):
        ValueWidgetBase.setObjectName(_fromUtf8("ValueWidgetBase"))
        ValueWidgetBase.resize(336, 153)
        ValueWidgetBase.setWindowTitle(QtGui.QApplication.translate("ValueWidgetBase", "Form", None, QtGui.QApplication.UnicodeUTF8))
        self.verticalLayout = QtGui.QVBoxLayout(ValueWidgetBase)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.cbxActive = QtGui.QCheckBox(ValueWidgetBase)
        self.cbxActive.setToolTip(QtGui.QApplication.translate("ValueWidgetBase", "(Shift+A) to toggle", None, QtGui.QApplication.UnicodeUTF8))
        self.cbxActive.setText(QtGui.QApplication.translate("ValueWidgetBase", "Active", None, QtGui.QApplication.UnicodeUTF8))
        self.cbxActive.setObjectName(_fromUtf8("cbxActive"))
        self.horizontalLayout.addWidget(self.cbxActive)
        self.line_2 = QtGui.QFrame(ValueWidgetBase)
        self.line_2.setFrameShape(QtGui.QFrame.VLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.horizontalLayout.addWidget(self.line_2)
        self.cbxDigits = QtGui.QCheckBox(ValueWidgetBase)
        self.cbxDigits.setText(QtGui.QApplication.translate("ValueWidgetBase", "Decimals", None, QtGui.QApplication.UnicodeUTF8))
        self.cbxDigits.setObjectName(_fromUtf8("cbxDigits"))
        self.horizontalLayout.addWidget(self.cbxDigits)
        self.spinBox = QtGui.QSpinBox(ValueWidgetBase)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.spinBox.sizePolicy().hasHeightForWidth())
        self.spinBox.setSizePolicy(sizePolicy)
        self.spinBox.setMaximum(99)
        self.spinBox.setProperty("value", 2)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.horizontalLayout.addWidget(self.spinBox)
        self.line = QtGui.QFrame(ValueWidgetBase)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.horizontalLayout.addWidget(self.line)
        self.cbxGraph = QtGui.QCheckBox(ValueWidgetBase)
        self.cbxGraph.setText(QtGui.QApplication.translate("ValueWidgetBase", "Graph", None, QtGui.QApplication.UnicodeUTF8))
        self.cbxGraph.setObjectName(_fromUtf8("cbxGraph"))
        self.horizontalLayout.addWidget(self.cbxGraph)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.graphControls = QtGui.QWidget(ValueWidgetBase)
        self.graphControls.setObjectName(_fromUtf8("graphControls"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.graphControls)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.cbxStats = QtGui.QCheckBox(self.graphControls)
        self.cbxStats.setToolTip(QtGui.QApplication.translate("ValueWidgetBase", "Compute min/max when layers are loaded", None, QtGui.QApplication.UnicodeUTF8))
        self.cbxStats.setText(QtGui.QApplication.translate("ValueWidgetBase", "Stats", None, QtGui.QApplication.UnicodeUTF8))
        self.cbxStats.setObjectName(_fromUtf8("cbxStats"))
        self.horizontalLayout_2.addWidget(self.cbxStats)
        self.label = QtGui.QLabel(self.graphControls)
        self.label.setText(QtGui.QApplication.translate("ValueWidgetBase", "Y min", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_2.addWidget(self.label)
        self.leYMin = QtGui.QLineEdit(self.graphControls)
        self.leYMin.setObjectName(_fromUtf8("leYMin"))
        self.horizontalLayout_2.addWidget(self.leYMin)
        self.label_2 = QtGui.QLabel(self.graphControls)
        self.label_2.setText(QtGui.QApplication.translate("ValueWidgetBase", "Y max", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_2.addWidget(self.label_2)
        self.leYMax = QtGui.QLineEdit(self.graphControls)
        self.leYMax.setObjectName(_fromUtf8("leYMax"))
        self.horizontalLayout_2.addWidget(self.leYMax)
        self.plotSelector = QtGui.QComboBox(self.graphControls)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotSelector.sizePolicy().hasHeightForWidth())
        self.plotSelector.setSizePolicy(sizePolicy)
        self.plotSelector.setToolTip(QtGui.QApplication.translate("ValueWidgetBase", "Select plotting toolkit", None, QtGui.QApplication.UnicodeUTF8))
        self.plotSelector.setObjectName(_fromUtf8("plotSelector"))
        self.horizontalLayout_2.addWidget(self.plotSelector)
        spacerItem1 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout.addWidget(self.graphControls)
        self.stackedWidget = QtGui.QStackedWidget(ValueWidgetBase)
        self.stackedWidget.setObjectName(_fromUtf8("stackedWidget"))
        self.page = QtGui.QWidget()
        self.page.setObjectName(_fromUtf8("page"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.page)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.tableWidget = QtGui.QTableWidget(self.page)
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(0)
        item = QtGui.QTableWidgetItem()
        item.setText(QtGui.QApplication.translate("ValueWidgetBase", "Layer", None, QtGui.QApplication.UnicodeUTF8))
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        item.setText(QtGui.QApplication.translate("ValueWidgetBase", "Value", None, QtGui.QApplication.UnicodeUTF8))
        self.tableWidget.setHorizontalHeaderItem(1, item)
        self.verticalLayout_2.addWidget(self.tableWidget)
        self.stackedWidget.addWidget(self.page)
        self.verticalLayout.addWidget(self.stackedWidget)

        self.retranslateUi(ValueWidgetBase)
        QtCore.QMetaObject.connectSlotsByName(ValueWidgetBase)

    def retranslateUi(self, ValueWidgetBase):
        item = self.tableWidget.horizontalHeaderItem(0)
        item = self.tableWidget.horizontalHeaderItem(1)

