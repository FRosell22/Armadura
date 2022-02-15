import sys
from PyQt5 import QtWidgets
from UI_Armadura import Ui_Armadura

class Ui_MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.ui = Ui_Armadura()
        self.ui.setupUi(self)




if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    my_app = Ui_MainWindow()
    my_app.show()
    sys.exit(app.exec_())
