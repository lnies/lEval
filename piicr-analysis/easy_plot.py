from PyQt5 import QtCore, QtGui, QtWidgets
import ui_easy_plot
# from analysis_functions import *
# from cross_checks_functions import *
import data_class

# python2 -m PyQt5.uic.pyuic easy_plot.ui -o ui_easy_plot.py -x

class MainWindow(QtWidgets.QMainWindow, ui_easy_plot.Ui_MainWindow):
	def __init__ (self):
		QtWidgets.QMainWindow.__init__(self)
		self.setupUi(self)
		# buttons
		self.bttn_exit.clicked.connect(self.my_exit)
		self.bttn_load.clicked.connect(self.my_load)
		self.bttn_plot.clicked.connect(self.run_plot)
		self.groupBox_2.clicked.connect(self.clear_isotope)
		self.groupBox_3.clicked.connect(self.clear_fit)
		self.groupBox_5.clicked.connect(self.clear_bounds)
		# text inputs
		self.isotope.textChanged.connect(self.set_isotope)
		self.delimiter.textChanged.connect(self.set_delimiter)
		self.folder_path.textChanged.connect(self.set_folder_path)
		self.a_min.textChanged.connect(self.set_fit_parameter)
		self.a_max.textChanged.connect(self.set_fit_parameter)
		self.b_min.textChanged.connect(self.set_fit_parameter)
		self.b_max.textChanged.connect(self.set_fit_parameter)
		self.c_min.textChanged.connect(self.set_fit_parameter)
		self.c_max.textChanged.connect(self.set_fit_parameter)
		# var
		self.isotope_name = ''
		self.delimiter_name = ';'
		self.folder_name = ''
		self.fit_method = 'x^2'
		self.column_dict = {'x-y': ['x', 'y'], 'x-y-yerr': ['x', 'y', 'yerr'], 'y-x': ['y', 'x'], 'x-xerr-y-yerr': ['x', 'xerr', 'y', 'yerr']}
		self.column = ['x', 'y']
		self.fit = False
		self.fit_bounds = []
		self.a_min_name = ''
		self.a_max_name = ''
		self.b_min_name = ''
		self.b_max_name = ''
		self.c_min_name = ''
		self.c_max_name = ''


	def my_exit(self):
		sys.exit()

	def my_load(self):
		self.folder_name = str(QtWidgets.QFileDialog.getOpenFileName(QtWidgets.QFileDialog(), 'Open File', filter='*.txt *.csv *.xlsx')[0])
		self.folder_path.setText(self.folder_name)

	def set_folder_path(self):
		self.a_min_name = str(self.a_min.text())
		self.a_max_name = str(self.a_max.text())
		self.b_min_name = str(self.b_min.text())
		self.b_max_name = str(self.b_max.text())
		self.c_min_name = str(self.c_min.text())
		self.c_max_name = str(self.c_max.text())

	def set_fit_parameter(self):
		self.folder_name = str(self.folder_path.text())

	def set_delimiter(self):
		self.delimiter_name = str(self.delimiter.text())
		if self.delimiter_name == '':
			self.delimiter_name = ';'

	def set_isotope(self):
		self.isotope_name = '$^{{{}}}${}'.format(filter(str.isdigit, str(self.isotope.text())), filter(str.isalpha, str(self.isotope.text())))
		if self.isotope_name == '$^{}$':
			self.isotope_name = ''

	def clear_isotope(self):
		self.isotope.setText('')

	def clear_bounds(self):
		self.a_min.setText('')
		self.a_max.setText('')
		self.b_min.setText('')
		self.b_max.setText('')
		self.c_min.setText('')
		self.c_max.setText('')

	def clear_fit(self):
		if self.groupBox_3.isChecked():
			self.fit = True
		else:
			self.fit = False

	def run_plot(self):
		print '\nFit                 ::', self.fit
		if self.isotope_name != '':
			print 'Isotope             ::', self.isotope_name
		jonas = data_class.Jonas()
		if self.a_min_name != '' and self.a_max_name != '' and self.b_min_name != '' and self.b_max_name != '':
			fit_parameter = [float(self.a_min_name),float(self.a_max_name),float(self.b_min_name),float(self.b_max_name)]
		elif self.a_min_name != '' and self.a_max_name != '' and self.b_min_name != '' and self.b_max_name != '' and self.c_min_name != '' and self.c_max_name != '':
			fit_parameter = [float(self.a_min_name),float(self.a_max_name),float(self.b_min_name),float(self.b_max_name),float(self.c_min_name),float(self.c_max_name)]
		else:
			fit_parameter = []

		if self.fit == True and self.isotope_name != '':
			self.fit_method = str(self.drop_fit.currentText())
			self.column = self.column_dict[str(self.drop_column.currentText())]
			jonas.plot(self.folder_name, self.column, fit_function=self.fit_method, isotope=self.isotope_name, fit_bounds=self.fit_bounds)
		elif self.fit == True:
			self.fit_method = str(self.drop_fit.currentText())
			self.column = self.column_dict[str(self.drop_column.currentText())]
			jonas.plot(self.folder_name, self.column, fit_function=self.fit_method)
		else:
			self.column = self.column_dict[str(self.drop_column.currentText())]
			jonas.plot(self.folder_name, self.column)


if __name__ == "__main__":
    try:
        app
    except:
        import sys
        app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())

