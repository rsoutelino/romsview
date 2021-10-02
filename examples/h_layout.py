# Filename: h_layout.py

"""Horizontal layout example."""

import sys

from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QWidget

app = QApplication(sys.argv)
window = QWidget()
window.setWindowTitle("Future")
layout = QHBoxLayout()
layout.addWidget(QPushButton("Tati"))
layout.addWidget(QPushButton("Rafa"))
layout.addWidget(QPushButton("Noosa"))
window.setLayout(layout)
window.show()
sys.exit(app.exec_())
