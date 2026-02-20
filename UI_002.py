#Multiple windows for ell, dtl and pillbox or arbitrary polygon
import sys, os
from typing import Optional

from PySide6.QtCore import Qt, QTimer, QByteArray, QProcess
from PySide6.QtGui import QIcon, QPixmap, QPalette, QColor
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QLabel,
    QLineEdit,
    QPushButton,
    QHBoxLayout,
    QVBoxLayout,
    QComboBox,
    QCheckBox,
    QProgressBar,
    QStatusBar,
    QFileDialog,
    QMessageBox,
)

# -----------------------------------------------------------------------------------------------------------------------------------------
FENICS_PY = "/opt/anaconda3/envs/fenics-py312/bin/python"   # path to the fenics env python
FENICS_SCRIPT = "/Users/straniero/Documents/Dphil/hyperfish/run.py"  # path to the script fenics must run
DTL_SCRIPT = "/Users/straniero/Documents/Dphil/hyperfish/run-DTL.py"  # path to the script fenics must run
# -----------------------------------------------------------------------------------------------------------------------------------------


class Color(QWidget):

    def __init__(self, color):
        super(Color, self).__init__()
        self.setAutoFillBackground(True)

        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor(color))
        self.setPalette(palette) 




#THIS NEXT PART DEFINES THE GUI
class Ell_window(QWidget):
    
    def __init__(self):
        super().__init__()
        self.resize(600, 400)

        self.interp_edit = QLineEdit(FENICS_PY)
    

        self.setWindowTitle('Mesh Generator')
        layout = QVBoxLayout()
        layout0 =QHBoxLayout()
        layout1 = QHBoxLayout()
        layout2=QHBoxLayout()
        layout3=QHBoxLayout()
        layout4=QHBoxLayout()
        layout5 = QHBoxLayout()
        layout6=QHBoxLayout()
        layout7=QHBoxLayout()
        layout8=QHBoxLayout()
        layout9=QHBoxLayout()
        layout10=QHBoxLayout()
        
        layout.addLayout(layout0) #control description
        layout.addLayout(layout1) #Input frequency
        layout.addLayout(layout6) #power entry
        layout.addLayout(layout7)
        layout.addLayout(layout9) # extra row 1 (two fields)
        layout.addLayout(layout10) # extra row 2 (two fields)
        layout.addLayout(layout3) 
        layout.addLayout(layout2) #start/stop buttons


        label1 =QLabel("Mesh Generator for ellipse (when pressing start, make sure all entries are filled)")
        label1.setStyleSheet("background-color: green")
        layout0.addWidget(label1)
        

        label = QLabel("Bore radius (m)")               
        layout1.addWidget(label)
        self.bore_rad = QLineEdit()
        layout1.addWidget(self.bore_rad)

        label3 =QLabel("equator radius (m)")             
        layout6.addWidget(label3)
        self.eq_radius = QLineEdit()                              
        layout6.addWidget(self.eq_radius)

        label4 =QLabel("wall angle (deg)")
        layout7.addWidget(label4)
        self.wall_angle = QLineEdit()                                
        layout7.addWidget(self.wall_angle)

        # --- Four additional QLineEdits, arranged two-per-row ---
        label5 = QLabel("Dome b (m)")
        self.dome_b= QLineEdit()
        label6 = QLabel("Dome a/b")
        self.dome_a_b = QLineEdit()
        label7 = QLabel("Iris b (m)")
        self.iris_b = QLineEdit()
        label8 = QLabel("Iris a/b")
        self.iris_a_b = QLineEdit()
        label9= QLabel('Lenght (m)')
        self.length = QLineEdit()




        # first extra row: two fields
        layout9.addWidget(label5)
        layout9.addWidget(self.dome_b)
        layout9.addWidget(label6)
        layout9.addWidget(self.dome_a_b)

        # second extra row: two fields
        layout10.addWidget(label7)
        layout10.addWidget(self.iris_b)
        layout10.addWidget(label8)
        layout10.addWidget(self.iris_a_b)

        layout3.addWidget(label9)
        layout3.addWidget(self.length)

        self.start_button = QPushButton("generate mesh")
        layout2.addWidget(self.start_button)
        self.start_button.clicked.connect(self.send_func)

        # prepare a QProcess to run the fenics script using the specified interpreter
        self.proc = QProcess(self)
        self.proc.finished.connect(self.on_fenics_finished)

        # set the layout on this QWidget so its contents are shown
        self.setLayout(layout)


    def send_func(self):
        if not all([self.dome_b.text(), self.dome_a_b.text(), self.iris_b.text(), self.iris_a_b.text(), self.bore_rad.text(), self.eq_radius.text(), self.wall_angle.text(), self.length.text()]):
            print("Error: All fields must be filled.")
            return
        dome_b = float(self.dome_b.text())
        dome_a_b = float(self.dome_a_b.text())
        iris_b = float(self.iris_b.text())
        iris_a_b = float(self.iris_a_b.text())
        bore_rad = float(self.bore_rad.text())
        eq_radius = float(self.eq_radius.text())
        wall_angle = float(self.wall_angle.text())
        length = float(self.length.text())

        
        print(f"Received parameters: dome_b={dome_b}, dome_a/b={dome_a_b}, iris_b={iris_b}, iris_a/b={iris_a_b}, bore_rad={bore_rad}, eq_radius={eq_radius}, wall_angle={wall_angle}, length={length}")

        # Verify interpreter and script exist
        if not os.path.exists(FENICS_PY):
            QMessageBox.critical(self, "Interpreter not found", f"Interpreter not found: {FENICS_PY}")
            return
        if not os.path.exists(FENICS_SCRIPT):
            QMessageBox.critical(self, "Script not found", f"Script not found: {FENICS_SCRIPT}")
            return

        # start the external process using the specified python executable
        if self.proc.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(self, "Already running", "Fenics process is already running")
            return

        self.start_button.setEnabled(False)

        # Pass GUI inputs as command-line arguments to the fenics script
        args = [
            FENICS_SCRIPT,
            str(dome_b),
            str(dome_a_b),
            str(iris_b),
            str(iris_a_b),
            str(bore_rad),
            str(eq_radius),
            str(wall_angle),
            str(length),
        ]

        # connect output readers for basic logging
        try:
            self.proc.readyReadStandardOutput.connect(lambda: print(bytes(self.proc.readAllStandardOutput()).decode()))
            self.proc.readyReadStandardError.connect(lambda: print(bytes(self.proc.readAllStandardError()).decode(), file=sys.stderr))
        except Exception:
            pass

        self.proc.start(FENICS_PY, args)
        if not self.proc.waitForStarted(1000):
            QMessageBox.critical(self, "Failed to start", "Could not start Fenics process")
            self.start_button.setEnabled(True)
            return
        print("Fenics process started with args:", args)

    def on_fenics_finished(self, exitCode, exitStatus):
        print(f"Fenics finished: code={exitCode}, status={exitStatus}")
        self.start_button.setEnabled(True)
#====================================================================================================================================


class DTL_window(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(800, 400)
        self.setWindowTitle('DTL Mesh Generator')
        
        layout = QVBoxLayout()
        layout1 = QVBoxLayout()
        layout2 = QVBoxLayout()
        layout1a = QHBoxLayout()
        layout1b = QHBoxLayout()
        layout1c = QHBoxLayout()
        layout1d = QHBoxLayout()
        layout1e = QHBoxLayout()



        layout.addLayout(layout1)
        #layout.addLayout(layout2)
        layout1.addLayout(layout1a)
        layout1.addLayout(layout1b)
        layout1.addLayout(layout1c)
        layout1.addLayout(layout1d)
        layout1.addLayout(layout1e)
        


        self.title = QLabel('DTL cavity mesh generator, when pressing start, make sure all entries are filled')

        self.label1 = QLabel('enter radius lenght (m)')
        self.label2 = QLabel('enter cavity lenght L (m)')
        self.label3 = QLabel('enter dtl angle (deg)')
        self.label4 = QLabel('enter DTL g (m)')
        self.label5 = QLabel('enter DTL d (m)')
        self.label6 = QLabel('enter bore radius (m)')

        self.R_dim = QLineEdit()
        self.Z_dim = QLineEdit()
        self.angle = QLineEdit()
        self.g_dtl = QLineEdit()
        self.d_dtl = QLineEdit()
        self.bore_radius = QLineEdit()

        layout1a.addWidget(self.title)
        layout1b.addWidget(self.label1)
        layout1b.addWidget(self.R_dim)
        layout1b.addWidget(self.label2)
        layout1b.addWidget(self.Z_dim)
        layout1c.addWidget(self.label4)
        layout1c.addWidget(self.g_dtl)
        layout1c.addWidget(self.label5)
        layout1c.addWidget(self.d_dtl)
        layout1d.addWidget(self.label3)
        layout1d.addWidget(self.angle)
        layout1d.addWidget(self.label6)
        layout1d.addWidget(self.bore_radius)

        self.start_button = QPushButton('generate mesh')
        layout1e.addWidget(self.start_button)
        self.start_button.clicked.connect(self.send_func)

        # prepare a QProcess to run the fenics script using the specified interpreter
        self.proc = QProcess(self)
        self.proc.finished.connect(self.on_fenics_finished)

        self.setLayout(layout)

    def send_func(self):
        if not all([self.R_dim.text(), self.Z_dim.text(), self.angle.text(), self.g_dtl.text(), self.d_dtl.text(), self.bore_radius.text()]):
            print("Error: All fields must be filled.")
            return
        R_dim = float(self.R_dim.text())
        Z_dim = float(self.Z_dim.text())
        angle = float(self.angle.text())
        g_dtl = float(self.g_dtl.text())
        d_dtl = float(self.d_dtl.text())
        bore_radius = float(self.bore_radius.text())

        
        print(f"Received parameters: R_dim={R_dim}, Z_dim={Z_dim}, angle={angle}, g_dtl={g_dtl}, d_dtl={d_dtl}, bore_radius={bore_radius}")

        # Verify interpreter and script exist
        if not os.path.exists(FENICS_PY):
            QMessageBox.critical(self, "Interpreter not found", f"Interpreter not found: {FENICS_PY}")
            return
        if not os.path.exists(DTL_SCRIPT):
            QMessageBox.critical(self, "Script not found", f"Script not found: {DTL_SCRIPT}")
            return

        # start the external process using the specified python executable
        if self.proc.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(self, "Already running", "Fenics process is already running")
            return

        self.start_button.setEnabled(False)

        # Pass GUI inputs as command-line arguments to the fenics script
        args = [
            DTL_SCRIPT,
            str(R_dim),
            str(Z_dim),
            str(angle),
            str(g_dtl),
            str(d_dtl),
            str(bore_radius),
        ]

        # connect output readers for basic logging
        try:
            self.proc.readyReadStandardOutput.connect(lambda: print(bytes(self.proc.readAllStandardOutput()).decode()))
            self.proc.readyReadStandardError.connect(lambda: print(bytes(self.proc.readAllStandardError()).decode(), file=sys.stderr))
        except Exception:
            pass

        self.proc.start(FENICS_PY, args)
        if not self.proc.waitForStarted(1000):
            QMessageBox.critical(self, "Failed to start", "Could not start Fenics process")
            self.start_button.setEnabled(True)
            return
        print("Fenics process started with args:", args)


    def on_fenics_finished(self, exitCode, exitStatus):
        print(f"Fenics finished: code={exitCode}, status={exitStatus}")
        self.start_button.setEnabled(True)

        
    

#====================================================================================================================================


class Pillbox_window(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(600, 400)
        self.setWindowTitle('Pillbox Mesh Generator')
        layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout2 = QHBoxLayout()

        layout.addLayout(layout1)
        layout.addLayout(layout2)

        self.label1 = QLabel('enter radius lenght (m)')
        self.label2 = QLabel('enter cavity lenght L (m)')

        self.R_dim = QLineEdit()
        self.Z_dim = QLineEdit()

        layout2.addWidget(self.label1)
        layout2.addWidget(self.R_dim)
        layout2.addWidget(self.label2)
        layout2.addWidget(self.Z_dim)

        self.setLayout(layout)

#====================================================================================================================================


class Polygon_window(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(600, 400)
        self.setWindowTitle('Polygon Mesh Generator')
        #layout


#====================================================================================================================================



class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(600, 400)

        self.setWindowTitle('Mesh Generator')
        layout_m = QVBoxLayout()
        layout0_m =QHBoxLayout()
        layout1_m = QHBoxLayout()

        layout_m.addLayout(layout0_m)
        layout_m.addLayout(layout1_m)

        label0_m = QLabel('Choose what type of cavity you want to generate:',)
        layout0_m.addWidget(label0_m)


        self.pillbox_button = QPushButton('Pillbox cavity')
        self.dtl_button = QPushButton('DTL cavity')
        self.ell_button = QPushButton('Elliptical cavity')
        self.polygon_button = QPushButton('Arbitrary polygon')

        layout1_m.addWidget(self.pillbox_button)
        layout1_m.addWidget(self.dtl_button)
        layout1_m.addWidget(self.ell_button)
        layout1_m.addWidget(self.polygon_button)

        self.ell_button.clicked.connect(self.ell_func)
        self.pillbox_button.clicked.connect(self.pillbox_func)
        self.dtl_button.clicked.connect(self.dtl_func)

        widget = QWidget()
        widget.setLayout(layout_m)
        self.setCentralWidget(widget)

    
    def ell_func(self, checked=False):
        print('clicked')
        # keep a reference so the window isn't garbage-collected
        self.ell_window = Ell_window()
        self.ell_window.show()
        print('done')
    
    def pillbox_func(self, checked=False):
        # keep a reference so the window isn't garbage-collected
        self.pill_window = Pillbox_window()
        self.pill_window.show()
    
    def dtl_func(self, checked=False):
        # keep a reference so the window isn't garbage-collected
        self.dtl_window = DTL_window()
        self.dtl_window.show()
       
    



app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()