import sys, os
from typing import Optional

from PySide6.QtCore import Qt, QTimer, QByteArray
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
FENICS_SCRIPT = "/Users/straniero/Documents/Dphil/hyperfish/ellipse.py"  # path to the script fenics must run
# -----------------------------------------------------------------------------------------------------------------------------------------


class Color(QWidget):

    def __init__(self, color):
        super(Color, self).__init__()
        self.setAutoFillBackground(True)

        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor(color))
        self.setPalette(palette) 

#THIS NEXT PART DEFINES THE GUI
class MainWindow(QMainWindow):
    
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
        layout.addLayout(layout2) #start/stop buttons
        layout.addLayout(layout3) 
        layout.addLayout(layout5) #monitor label
        layout.addLayout(layout8)
        layout.addLayout(layout4) #monitor graph

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

        self.start_button = QPushButton("generate mesh")
        layout2.addWidget(self.start_button)
        self.start_button.clicked.connect(self.send_func)

        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)


    def send_func(self):
        if not all([self.dome_b.text(), self.dome_a_b.text(), self.iris_b.text(), self.iris_a_b.text(), self.bore_rad.text(), self.eq_radius.text(), self.wall_angle.text()]):
            print("Error: All fields must be filled.")
            return
        dome_b = float(self.dome_b.text())
        dome_a_b = float(self.dome_a_b.text())
        iris_b = float(self.iris_b.text())
        iris_a_b = float(self.iris_a_b.text())
        bore_rad = float(self.bore_rad.text())
        eq_radius = float(self.eq_radius.text())
        wall_angle = float(self.wall_angle.text())

        
        print(f"Received parameters: dome_b={dome_b}, dome_a/b={dome_a_b}, iris_b={iris_b}, iris_a/b={iris_a_b}, bore_rad={bore_rad}, eq_radius={eq_radius}, wall_angle={wall_angle}")




app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()