# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cpxt/gui/ui/main_ui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainUI(object):
    def setupUi(self, MainUI):
        MainUI.setObjectName("MainUI")
        MainUI.setWindowModality(QtCore.Qt.WindowModal)
        MainUI.resize(538, 442)
        MainUI.setWindowTitle("ClusterPyXT (cpxt)")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("cpxt/assets/cpxt.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainUI.setWindowIcon(icon)
        self.centralwidget = QtWidgets.QWidget(MainUI)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.cluster_list_box = QtWidgets.QListWidget(self.centralwidget)
        self.cluster_list_box.setObjectName("cluster_list_box")
        self.gridLayout.addWidget(self.cluster_list_box, 0, 0, 1, 1)
        self.continue_button = QtWidgets.QPushButton(self.centralwidget)
        self.continue_button.setObjectName("continue_button")
        self.gridLayout.addWidget(self.continue_button, 1, 0, 1, 1)
        self.new_cluster_button = QtWidgets.QPushButton(self.centralwidget)
        self.new_cluster_button.setObjectName("new_cluster_button")
        self.gridLayout.addWidget(self.new_cluster_button, 2, 0, 1, 1)
        MainUI.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainUI)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 538, 20))
        self.menubar.setObjectName("menubar")
        self.menu_File = QtWidgets.QMenu(self.menubar)
        self.menu_File.setObjectName("menu_File")
        MainUI.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainUI)
        self.statusbar.setObjectName("statusbar")
        MainUI.setStatusBar(self.statusbar)
        self.action_Settings = QtWidgets.QAction(MainUI)
        self.action_Settings.setObjectName("action_Settings")
        self.actionE_xit = QtWidgets.QAction(MainUI)
        self.actionE_xit.setObjectName("actionE_xit")
        self.menu_File.addAction(self.action_Settings)
        self.menu_File.addSeparator()
        self.menu_File.addAction(self.actionE_xit)
        self.menubar.addAction(self.menu_File.menuAction())

        self.retranslateUi(MainUI)
        QtCore.QMetaObject.connectSlotsByName(MainUI)

    def retranslateUi(self, MainUI):
        _translate = QtCore.QCoreApplication.translate
        self.continue_button.setText(_translate("MainUI", "Continue Cluster"))
        self.new_cluster_button.setText(_translate("MainUI", "New Cluster"))
        self.menu_File.setTitle(_translate("MainUI", "&File"))
        self.action_Settings.setText(_translate("MainUI", "&Settings"))
        self.actionE_xit.setText(_translate("MainUI", "E&xit"))

