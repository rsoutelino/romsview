# Filename: hello.py

"""Simple Hello World example with PyQt5."""

import sys

# 1. Import `QApplication` and all the required widgets
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QWidget


app = QApplication(sys.argv)

# 3. Create an instance of your application's GUI
window = QWidget()
window.setWindowTitle("Tati App")
window.setGeometry(100, 100, 1280, 180)
window.move(60, 15)
helloMsg = QLabel("<h1>Hello Tati!</h1>", parent=window)
helloMsg.move(600, 90)

# 4. Show your application's GUI
window.show()

# 5. Run your application's event loop (or main loop)
sys.exit(app.exec_())
