
# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017/2018. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

from PyQt5 import QtCore, QtGui, QtWidgets
import ui_main
from piicr_analysis import *
from piicr_cross_checks import *
import ellipse_selector
import os
import sys

# pyuic5 -x piicr_analysis_gui.ui -o ui_main.py

class MainWindow(QtWidgets.QMainWindow, ui_main.Ui_MainWindow):
	def __init__ (self):
		QtWidgets.QMainWindow.__init__(self)
		self.setupUi(self)
		# buttons
		self.exit_bttn.clicked.connect(self.my_exit)
		self.load.clicked.connect(self.my_load)
		self.load_5.clicked.connect(self.my_loads)
		self.load_2.clicked.connect(self.my_load_singles)
		self.load_3.clicked.connect(self.my_load_parameter_fix)
		self.run.clicked.connect(self.run_analysis)
		self.run_2.clicked.connect(self.run_cross_checks)
		self.run_3.clicked.connect(self.run_single_fit)
		self.run_8.clicked.connect(self.run_combined_cross_checks)
		self.groupBox.clicked.connect(self.clear_position)
		self.groupBox_2.clicked.connect(self.clear_tof)
		self.groupBox_3.clicked.connect(self.clear_z_class)
		self.groupBox_4.clicked.connect(self.clear_parameter_fix)
		self.groupBox_5.clicked.connect(self.clear_ground_or_isomer)
		# text inputs
		self.isotope1.textChanged.connect(self.set_isotope1)
		self.isotope2.textChanged.connect(self.set_isotope2)
		self.isotope3.textChanged.connect(self.set_isotope3)
		self.isotope4.textChanged.connect(self.set_isotope4)
		self.isotope5.textChanged.connect(self.set_isotope5)
		self.isotope1_2.textChanged.connect(self.set_isotope1_2)
		self.x_min.textChanged.connect(self.set_x_min)
		self.x_max.textChanged.connect(self.set_x_max)
		self.y_min.textChanged.connect(self.set_y_min)
		self.y_max.textChanged.connect(self.set_y_max)
		self.z_class_min.textChanged.connect(self.set_z_class_min)
		self.z_class_max.textChanged.connect(self.set_z_class_max)
		self.tof_min.textChanged.connect(self.set_tof_min)
		self.tof_max.textChanged.connect(self.set_tof_max)
		self.ref.textChanged.connect(self.set_ref)
		self.meas.textChanged.connect(self.set_meas)
		self.folder_path.textChanged.connect(self.set_folder_path)
		self.folder_path_2.textChanged.connect(self.set_parameter_fix)
		self.folder_paths.textChanged.connect(self.set_folder_paths)
		self.file_names.textChanged.connect(self.set_folder_path_singles)
		self.ground_or_isomer.textChanged.connect(self.set_ground_or_isomer)
		self.isomer_mass_excess.textChanged.connect(self.set_isomer_mass_excess)
		self.isomer_mass_excess_unc.textChanged.connect(self.set_isomer_mass_excess_unc)

		self.isomer_mass_excess.setText('-62910')
		self.isomer_mass_excess_unc.setText('150')

		self.isotope1_name = ''
		self.isotope2_name = ''
		self.isotope3_name = ''
		self.isotope4_name = ''
		self.isotope5_name = ''
		self.isotope1_2_name = ''
		self.x_min_name = ''
		self.x_max_name = ''
		self.y_min_name = ''
		self.y_max_name = ''
		self.z_class_min_name = ''
		self.z_class_max_name = ''
		self.tof_min_name = ''
		self.tof_max_name = ''
		self.ref_name = ''
		self.meas_name = ''
		self.folder_names = ''
		self.file_path_singles = []
		self.FWHM_parameter_fix = {'p1': {}, 'p2': {}, 'c': {}}
		self.all_peak_pos_c_file_path = ''
		self.isomer_analysis = False
		self.ellipse_selection = False
		self.check_status = {0: False,
							 2: True}
		self.cross_checks_isomer_analysis_str = ''
		self.isomer_state__empty_str_if_not_isotope_analysis = ''
		self.isomer_ME = 0
		self.isomer_ME_unc = 0


	def my_exit(self):
		sys.exit()

	def my_load(self):
		self.folder_name = str(QtWidgets.QFileDialog.getExistingDirectory(QtWidgets.QFileDialog(), 'Open File'))
		self.folder_path.setText(self.folder_name)

	def my_load_parameter_fix(self):
		self.all_peak_pos_c_file_path = str(QtWidgets.QFileDialog.getOpenFileName(QtWidgets.QFileDialog(), 'Open File', '', 'Comma separated value file (*.csv)')[0])
		self.folder_path_2.setText(self.all_peak_pos_c_file_path)

	def my_loads(self):
		file_dialog = QtWidgets.QFileDialog()
		file_dialog.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
		file_dialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog, True)
		file_view = file_dialog.findChild(QtWidgets.QListView, 'listView')

		# make it possible to select multiple directories:
		if file_view:
		    file_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
		f_tree_view = file_dialog.findChild(QtWidgets.QTreeView)
		if f_tree_view:
		    f_tree_view.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)

		if file_dialog.exec_():
			hilf_folder_names = file_dialog.selectedFiles()
			self.folder_names = [str(x) for x in hilf_folder_names]
			hilf_folder_names_2 = ''
			for fns in self.folder_names:
				if hilf_folder_names_2 == '':
					hilf_folder_names_2 += fns
				else:
					hilf_folder_names_2 += ('\n' + fns)
			self.folder_paths.setText(hilf_folder_names_2)

	def my_load_singles(self):
		self.file_path_singles = []
		self.file_names.clear()
		self.file_name = list(QtWidgets.QFileDialog.getOpenFileNames(QtWidgets.QFileDialog(), 'Open File', '', 'Text files (*.txt)'))
		for i in self.file_name[0]:
			self.file_path_singles.append(str(i))

		for i in self.file_path_singles:
			self.file_names.append(i)

	def set_isotope1(self):
		self.isotope1_name = str(self.isotope1.text())
		self.ref.setText(self.isotope1.text())

	def set_isotope2(self):
		self.isotope2_name = str(self.isotope2.text())
		self.meas.setText(self.isotope2.text())

	def set_isotope3(self):
		self.isotope3_name = str(self.isotope3.text())

	def set_isotope4(self):
		self.isotope4_name = str(self.isotope4.text())

	def set_isotope5(self):
		self.isotope5_name = str(self.isotope5.text())

	def set_isotope1_2(self):
		self.isotope1_2_name = str(self.isotope1_2.text())

	def set_x_min(self):
		self.x_min_name = str(self.x_min.text())

	def set_x_max(self):
		self.x_max_name = str(self.x_max.text())

	def set_y_min(self):
		self.y_min_name = str(self.y_min.text())

	def set_y_max(self):
		self.y_max_name = str(self.y_max.text())

	def set_z_class_min(self):
		self.z_class_min_name = str(self.z_class_min.text())

	def set_z_class_max(self):
		self.z_class_max_name = str(self.z_class_max.text())

	def set_tof_min(self):
		self.tof_min_name = str(self.tof_min.text())

	def set_tof_max(self):
		self.tof_max_name = str(self.tof_max.text())

	def set_ref(self):
		self.ref_name = str(self.ref.text())

	def set_meas(self):
		self.meas_name = str(self.meas.text())

	def set_folder_path(self):
		self.folder_name = str(self.folder_path.text())

	def set_parameter_fix(self):
		self.all_peak_pos_c_file_path = str(self.folder_path_2.text())

	def set_folder_paths(self):
		hilf_folder_names = (self.folder_paths.toPlainText())
		self.folder_names = [str(x) for x in hilf_folder_names.splitlines()]

	def set_folder_path_singles(self):	#testtesttest
		hilf_folder_names = (self.file_names.toPlainText())
		self.file_path_singles = [str(x) for x in hilf_folder_names.splitlines()]

	def set_file_path(self):
		self.file_name = str(self.file_path.text())

	def set_ground_or_isomer(self):
		self.cross_checks_isomer_analysis_str = str(self.ground_or_isomer.text())
		if self.cross_checks_isomer_analysis_str == 'ground':
			self.isomer_state__empty_str_if_not_isotope_analysis = '_ground'
		elif self.cross_checks_isomer_analysis_str == 'isomer':
			self.isomer_state__empty_str_if_not_isotope_analysis = '_isomer'
		else:
			self.isomer_state__empty_str_if_not_isotope_analysis = ''

	def set_isomer_mass_excess(self):
		self.isomer_ME = str(self.isomer_mass_excess.text())

	def set_isomer_mass_excess_unc(self):
		self.isomer_ME_unc = str(self.isomer_mass_excess_unc.text())

	def clear_tof(self):
		self.tof_min.setText('')
		self.tof_max.setText('')

	def clear_z_class(self):
		self.z_class_min.setText('')
		self.z_class_max.setText('')

	def clear_position(self):
		self.x_min.setText('')
		self.x_max.setText('')
		self.y_min.setText('')
		self.y_max.setText('')

	def clear_parameter_fix(self):
		self.folder_path_2.setText('')
		self.all_peak_pos_c_file_path = ''

	def clear_ground_or_isomer(self):
		self.ground_or_isomer.setText('')
		self.cross_checks_isomer_analysis_str = ''

		self.isomer_mass_excess.setText('')
		self.isomer_mass_excess_unc.setText('')
		self.isomer_ME = 0
		self.isomer_ME_unc = 0


	def run_analysis(self):
		isotopes = []
		if not self.isotope1_name == '' and type(self.isotope1_name) == str:
			isotopes.append(self.isotope1_name)
		if not self.isotope2_name == '' and type(self.isotope2_name) == str:
			isotopes.append(self.isotope2_name)
		if not self.isotope3_name == '' and type(self.isotope3_name) == str:
			isotopes.append(self.isotope3_name)
		if not self.isotope4_name == '' and type(self.isotope4_name) == str:
			isotopes.append(self.isotope4_name)
		if not self.isotope5_name == '' and type(self.isotope5_name) == str:
			isotopes.append(self.isotope5_name)
		print 'Folder name              :: ', self.folder_name

		self.isomer_analysis = self.check_status[self.checkBox.checkState()]
		self.ellipse_selection = self.check_status[self.checkBox_ellipse.checkState()]

		folder_creation_analysis(self.folder_name, isotopes)
		if self.all_peak_pos_c_file_path != '':
			self.FWHM_parameter_fix = get_FWHM_parameter(self.all_peak_pos_c_file_path)

		for i in isotopes:
			analysis_main(self.folder_name+'/'+i, [self.x_min_name, self.x_max_name, self.y_min_name, self.y_max_name], [self.tof_min_name, self.tof_max_name], [self.z_class_min_name, self.z_class_max_name], self.FWHM_parameter_fix, self.isomer_analysis, self.ellipse_selection)


	def run_single_fit(self):
		print 'Currently in beta.'

		folder_creation_analysis(os.path.dirname(self.file_path_singles[0]), [self.isotope1_2_name])
		for fps in self.file_path_singles:
			print fps

	def run_cross_checks(self):
		folder_creation_cross_checks(self.folder_name, self.ref_name, self.meas_name)
		cross_checks(self.folder_name, self.ref_name, self.meas_name, self.isomer_state__empty_str_if_not_isotope_analysis, float(self.isomer_ME), float(self.isomer_ME_unc))

	def run_combined_cross_checks(self):
		print self.folder_names
		plot_all_cross_checks(self.folder_names)

if __name__ == "__main__":
    print '\nPI-ICR analysis software :: jonas.karthein@cern.ch\n'
    import sys
    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())

