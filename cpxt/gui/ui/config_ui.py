# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cpxt/gui/ui/config_ui.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ConfigUI(object):
    def setupUi(self, ConfigUI):
        ConfigUI.setObjectName("ConfigUI")
        ConfigUI.resize(612, 232)
        self.directory_text = QtWidgets.QLineEdit(ConfigUI)
        self.directory_text.setGeometry(QtCore.QRect(100, 10, 401, 25))
        self.directory_text.setObjectName("directory_text")
        self.directory_label = QtWidgets.QLabel(ConfigUI)
        self.directory_label.setGeometry(QtCore.QRect(9, 14, 87, 17))
        self.directory_label.setObjectName("directory_label")
        self.find_directory_button = QtWidgets.QPushButton(ConfigUI)
        self.find_directory_button.setGeometry(QtCore.QRect(510, 10, 88, 25))
        self.find_directory_button.setObjectName("find_directory_button")
        self.astro_query_check = QtWidgets.QCheckBox(ConfigUI)
        self.astro_query_check.setGeometry(QtCore.QRect(9, 50, 357, 23))
        self.astro_query_check.setChecked(True)
        self.astro_query_check.setObjectName("astro_query_check")
        self.obs_download_label = QtWidgets.QLabel(ConfigUI)
        self.obs_download_label.setGeometry(QtCore.QRect(9, 84, 391, 17))
        self.obs_download_label.setObjectName("obs_download_label")
        self.obs_download_combo = QtWidgets.QComboBox(ConfigUI)
        self.obs_download_combo.setGeometry(QtCore.QRect(419, 84, 181, 25))
        self.obs_download_combo.setObjectName("obs_download_combo")
        self.obs_download_combo.addItem("")
        self.obs_download_combo.addItem("")
        self.obs_download_combo.addItem("")
        self.obs_download_combo.addItem("")
        self.obs_download_combo.addItem("")
        self.num_cores_label = QtWidgets.QLabel(ConfigUI)
        self.num_cores_label.setGeometry(QtCore.QRect(9, 120, 404, 17))
        self.num_cores_label.setObjectName("num_cores_label")
        self.num_cores_text = QtWidgets.QLineEdit(ConfigUI)
        self.num_cores_text.setGeometry(QtCore.QRect(419, 120, 181, 25))
        self.num_cores_text.setObjectName("num_cores_text")
        self.log_level_label = QtWidgets.QLabel(ConfigUI)
        self.log_level_label.setGeometry(QtCore.QRect(9, 156, 241, 17))
        self.log_level_label.setObjectName("log_level_label")
        self.log_level_combo = QtWidgets.QComboBox(ConfigUI)
        self.log_level_combo.setGeometry(QtCore.QRect(419, 156, 181, 25))
        self.log_level_combo.setObjectName("log_level_combo")
        self.log_level_combo.addItem("")
        self.log_level_combo.addItem("")
        self.log_level_combo.addItem("")
        self.log_level_combo.addItem("")
        self.log_level_combo.addItem("")
        self.save_button = QtWidgets.QPushButton(ConfigUI)
        self.save_button.setGeometry(QtCore.QRect(9, 192, 591, 25))
        self.save_button.setObjectName("save_button")

        self.retranslateUi(ConfigUI)
        QtCore.QMetaObject.connectSlotsByName(ConfigUI)

    def retranslateUi(self, ConfigUI):
        _translate = QtCore.QCoreApplication.translate
        ConfigUI.setWindowTitle(_translate("ConfigUI", "ClusterPyXT Configuration"))
        self.directory_label.setText(_translate("ConfigUI", "Data directory: "))
        self.find_directory_button.setText(_translate("ConfigUI", "Find Directory"))
        self.astro_query_check.setText(_translate("ConfigUI", "Automatically query NASA HEASARC and NED for nH and z?"))
        self.obs_download_label.setText(_translate("ConfigUI", "Number of observations to download at one time? (1-5)"))
        self.obs_download_combo.setItemText(0, _translate("ConfigUI", "1"))
        self.obs_download_combo.setItemText(1, _translate("ConfigUI", "2"))
        self.obs_download_combo.setItemText(2, _translate("ConfigUI", "3"))
        self.obs_download_combo.setItemText(3, _translate("ConfigUI", "4"))
        self.obs_download_combo.setItemText(4, _translate("ConfigUI", "5"))
        self.num_cores_label.setText(_translate("ConfigUI", "Number of CPU cores to use for parallel processing? (Max: {cpu_count}):"))
        self.log_level_label.setText(_translate("ConfigUI", "Log level (from most to least verbose):"))
        self.log_level_combo.setItemText(0, _translate("ConfigUI", "10"))
        self.log_level_combo.setItemText(1, _translate("ConfigUI", "20"))
        self.log_level_combo.setItemText(2, _translate("ConfigUI", "30"))
        self.log_level_combo.setItemText(3, _translate("ConfigUI", "40"))
        self.log_level_combo.setItemText(4, _translate("ConfigUI", "50"))
        self.save_button.setText(_translate("ConfigUI", "Save Configuration + Continue"))

