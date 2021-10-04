#! /usr/bin/env python3
from settings import *
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT as NavigationToolbar,
)
from PyQt5.QtWidgets import (
    QMainWindow,
    QApplication,
    QWidget,
    QComboBox,
    QLineEdit,
    QDialog,
    QSlider,
    QToolBar,
    QGroupBox,
    QHBoxLayout,
    QGridLayout,
    QStatusBar,
    QVBoxLayout,
    QMessageBox,
    QLabel,
    QFileDialog,
)
from PyQt5.QtGui import QPalette, QColor
from PyQt5.QtCore import Qt
import sys
from os.path import basename
from functools import partial

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr

# import cartopy.crs as ccrs

mpl.use("Qt5Agg")


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=10, height=10, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


class Ui(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setWindowTitle("ROMSView")
        self._state = AppState()
        self.generalLayout = QHBoxLayout()
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.centralWidget.setLayout(self.generalLayout)

        self._createMenu()
        self._createToolBar()
        self._createSideBar()
        self._createMplCanvas()
        self._createStatusBar()

    def _createMenu(self):
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction("&Load roms_grd.nc", not_found_dialog)
        self.fileMenu.addAction("&Load roms_clm.nc", not_found_dialog)
        self.fileMenu.addAction("&Load roms_ini.nc", not_found_dialog)
        self.fileMenu.addAction("&Load roms_his.nc", not_found_dialog)
        self.fileMenu.addAction("&Quit", self.close)

        self.toolsMenu = self.menuBar().addMenu("&Plot")
        self.toolsMenu.addAction("&Hslice", not_found_dialog)
        self.toolsMenu.addAction("&Vslice", not_found_dialog)

    def _createToolBar(self):
        tools = QToolBar()
        self.addToolBar(tools)
        for key in RomsNCFiles.__dataclass_fields__.keys():
            tools.addAction(key.upper(), partial(
                self.openFile, f"*_{key}*.nc"))

    def _createSideBar(self):
        self.sideBarLayout = QVBoxLayout()

        self._createPlotSelector()
        self._createVarSelector()
        self._createTimeSelector()
        self._createLevSelector()
        self._createCbarSelector()
        self._createAlphaSelector()
        self._createRangeBox()

        widget = QWidget()
        widget.setLayout(self.sideBarLayout)
        widget.setFixedWidth(185)
        self.generalLayout.addWidget(widget)

    def _createPlotSelector(self):
        self.plotSelector = QComboBox()
        self.plotSelector.setToolTip(
            "What to plot when clicking horizontal slice point (s)")
        self.plotSelector.addItems(["Timeseries on click", "Vslice on click"])
        self.plotSelector.setDisabled(True)
        self.sideBarLayout.addWidget(self.plotSelector)

    def _createVarSelector(self):
        self.varSelector = QComboBox()
        self.varSelector.setToolTip("Variables")
        self.varSelector.addItem("Variables")
        self.varSelector.setDisabled(True)
        self.varSelector.activated[str].connect(self.toggle_var)
        self.sideBarLayout.addWidget(self.varSelector)

    def _createTimeSelector(self):
        self.timeSelector = QComboBox()
        self.timeSelector.setToolTip("Times")
        self.timeSelector.addItem("Times")
        self.timeSelector.setDisabled(True)
        self.timeSelector.activated[str].connect(self.toggle_time)
        self.sideBarLayout.addWidget(self.timeSelector)

    def _createLevSelector(self):
        self.levSelector = QComboBox()
        self.levSelector.setToolTip("Levels")
        self.levSelector.addItem("Levels")
        self.levSelector.setDisabled(True)
        self.levSelector.activated[str].connect(self.toggle_lev)
        self.sideBarLayout.addWidget(self.levSelector)

    def _createCbarSelector(self):
        self.cbarSelector = QComboBox()
        self.cbarSelector.setToolTip("Colorbars")
        self.cbarSelector.addItems(["viridis", "jet", "RdBu_r"])
        self.cbarSelector.activated[str].connect(self.set_colorbar)
        self.cbarSelector.setDisabled(True)
        self.sideBarLayout.addWidget(self.cbarSelector)

    def _createAlphaSelector(self):
        alpha = QSlider(Qt.Horizontal)
        alpha.setValue(100)
        alpha.valueChanged[int].connect(self.set_alpha)
        self.sideBarLayout.addWidget(alpha)

    def _createRangeBox(self):
        self.rangeBox = QGroupBox()
        vmin_label = QLabel("Vmin")
        vmax_label = QLabel("Vmax")
        self.vmin = QLineEdit()
        self.vmax = QLineEdit()
        rangeLayout = QGridLayout()
        rangeLayout.addWidget(vmax_label, 0, 0, 1, 1)
        rangeLayout.addWidget(self.vmax, 0, 1, 1, 1)
        rangeLayout.addWidget(vmin_label, 1, 0, 1, 1)
        rangeLayout.addWidget(self.vmin, 1, 1, 1, 1)
        self.rangeBox.setLayout(rangeLayout)
        self.rangeBox.setFixedHeight(100)
        self.vmin.returnPressed.connect(self.set_range)
        self.vmax.returnPressed.connect(self.set_range)
        self.sideBarLayout.addWidget(self.rangeBox)
        self.rangeBox.setDisabled(True)

    def _createMplCanvas(self):
        self.mplcanvas = MplCanvas(self, width=5, height=4, dpi=100)
        cid = self.mplcanvas.mpl_connect(
            'button_press_event', self.timeseries_or_vslice)
        self.init_plot()

        # Create toolbar, passing canvas as first parameter, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(self.mplcanvas, self)

        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(self.mplcanvas)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QWidget()
        widget.setLayout(layout)
        self.generalLayout.addWidget(widget)

    def _createStatusBar(self):
        self.status = QStatusBar()
        self.status.showMessage("Ready...")
        self.setStatusBar(self.status)

    def openFile(self, pattern="*.nc"):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "QFileDialog.getOpenFileName()",
            "/source/roms-py/tests",
            f"NetCDF Files ({pattern});;All Files (*)",
            options=options,
        )
        if filename:
            self.onOpenFile(filename)

    def onOpenFile(self, filename):
        self.status.showMessage(f"Current file: {filename}")
        self._load_dataset(filename)
        self._state.filetype = detect_roms_file(filename)
        # getting a representative var based on settings.rep_var
        rep_var = getattr(REP_VAR, self._state.filetype)
        self._state.da = last2d(self._state.ds[rep_var])
        self.hslice(var_changed=True)

    def _reset_mpl_axes(self):
        for ax in self.mplcanvas.figure.axes:
            ax.remove()

        self.mplcanvas.axes = self.mplcanvas.figure.add_subplot(111)

    def init_plot(self):
        img = plt.imread("./icons/welcome.png")
        self._plot = self.mplcanvas.axes.imshow(img)
        self.mplcanvas.axes.set_axis_off()
        self.mplcanvas.draw()

    def _load_dataset(self, filename):
        self._state.ds = xr.open_dataset(filename)

    def hslice(self, var_changed=False):
        self._reset_mpl_axes()

        self._plot = self._state.da.plot(ax=self.mplcanvas.axes)
        if hasattr(self._plot, "set_cmap"):
            self.cbarSelector.setEnabled(True)
            self.rangeBox.setEnabled(True)
            self.plotSelector.setEnabled(True)

            if (
                hasattr(self._state, "vmin") and hasattr(self._state, "vmax")
            ) and not var_changed:
                self._plot.set_norm(
                    mpl.colors.Normalize(self._state.vmin, self._state.vmax)
                )
            else:
                self._reset_range(np.nanmin(self._state.da),
                                  np.nanmax(self._state.da))

            self.set_colorbar(cbar=self.cbarSelector.currentText())
        else:
            self.cbarSelector.setDisabled(True)
            self.rangeBox.setDisabled(True)

        self.mplcanvas.draw()

        self._update_vars()
        self._update_times()
        self._update_levels()

    def timeseries_or_vslice(self, evt):
        if evt.inaxes != self.mplcanvas.axes:
            return

        if 'Vslice' in self.plotSelector.currentText():
            dialog = VsliceDialog(title="Vertical Slice")

        if 'Tseries' in self.plotSelector.currentText():
            dialog = TseriesDialog(title="Time Series")

        dialog.setGeometry(2000, 60, 800, 400)
        dialog.show()

        try:
            self.dialogs.append(dialog)
        except AttributeError:
            self.dialogs = [dialog]

    def _reset_range(self, vmin, vmax):
        self._state.vmin = vmin
        self._state.vmax = vmax
        self.vmin.setText(f"{vmin:0.6f}")
        self.vmax.setText(f"{vmax:0.6f}")

    def _update_vars(self):
        self.varSelector.setEnabled(True)
        self.varSelector.clear()
        self.varSelector.addItems(self._state.ds.data_vars.keys())
        self.varSelector.setCurrentText(self._state.da.name)

    def _update_times(self):
        for dim in self._state.ds.dims.keys():
            if (
                "time" in dim
                and self._state.filetype not in ["grd"]
                and dim in self._state.da.coords.keys()
            ):
                self.timeSelector.setEnabled(True)
                self.timeSelector.clear()
                times = [numpydatetime2str(t)
                         for t in self._state.ds[dim].values]
                self.timeSelector.addItems(times)
                current = numpydatetime2str(self._state.da[dim].values)
                self.timeSelector.setCurrentText(current)
                break

            self.timeSelector.setDisabled(True)

    def _update_levels(self):
        for dim, val in self._state.ds.dims.items():
            if (
                "s_rho" in dim
                and self._state.filetype not in ["grd", "bry"]
                and dim in self._state.da.coords.keys()
            ):
                self.levSelector.setEnabled(True)
                self.levSelector.clear()
                levels = [str(l) for l in self._state.ds[dim].values]
                self.levSelector.addItems(levels)
                self.levSelector.setCurrentText(
                    str(self._state.da[dim].values))
                break

            self.levSelector.setDisabled(True)

    def toggle_var(self, var):
        _slice = {}
        # need to remove dimensions that don't exist in the new var
        # if that's the case (Ex, toggling from 3D to 2D var)
        for dim, val in self._state.current_slice.items():
            if dim in self._state.ds[var].dims:
                _slice[dim] = val

        self._state.da = last2d(self._state.ds[var].sel(**_slice))
        self.hslice(var_changed=True)

    def toggle_time(self, timestamp):
        _slice = self._state.current_slice.copy()
        for key in _slice.keys():
            if "time" in key:
                _slice[key] = self.timeSelector.currentText()

                self._state.da = last2d(
                    self._state.ds[self._state.var].sel(**_slice))
                self.hslice()
                break

    def toggle_lev(self, lev):
        _slice = self._state.current_slice.copy()
        for key in _slice.keys():
            if "s_rho" in key:
                _slice[key] = self.levSelector.currentText()

                self._state.da = last2d(
                    self._state.ds[self._state.var].sel(**_slice))
                self.hslice()
                break

    def set_range(self):
        try:
            vmin = float(self.vmin.text())
            self._state.vmin = vmin
        except:
            vmin = self._state.vmin

        try:
            vmax = float(self.vmax.text())
            self._state.vmax = vmax
        except:
            vmax = self._state.vmax

        self._plot.set_norm(mpl.colors.Normalize(
            self._state.vmin, self._state.vmax))
        self.mplcanvas.draw()

    def set_colorbar(self, cbar):
        if hasattr(self._plot, "set_cmap"):
            self._plot.set_cmap(getattr(plt.cm, cbar))
            self.mplcanvas.draw()
        else:
            not_found_dialog("Colorbar does not apply to this plot")

    def set_alpha(self, val):
        self._plot.set_alpha(val / 100)
        self.mplcanvas.draw()


class VsliceDialog(Ui, QDialog):
    def __init__(self, title='ROMSView dialog', *args, **kwargs):
        # super().__init__(*args, **kwargs)
        QDialog.__init__(self, *args, **kwargs)
        self.setWindowTitle(title)
        self.generalLayout = QHBoxLayout()
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.centralWidget.setLayout(self.generalLayout)

        # self._createMenu()
        # self._createToolBar()
        self._createSideBar()
        self._createMplCanvas()
        self._createStatusBar()

    def vslice(self):
        # leaving some hints on how to expand lon/lat dims for vslices
        # pcolormesh(
        #     self._state.ds.lat_rho.isel(xi_rho=0).expand_dims(
        #         {"s_rho": range(self._state.ds.dims["s_rho"])}, 0
        #     ),
        #     self._state.ds.z_rho.isel(ocean_time=0, xi_rho=0),
        #     self._state.dstemp.isel(ocean_time=0, xi_rho=0),
        # )
        # coords can't be masked or have nans (z_rho), so need to work that out

        pass


class TseriesDialog(Ui, QDialog):
    def __init__(self, title='ROMSView dialog', *args, **kwargs):
        # super().__init__(*args, **kwargs)
        QDialog.__init__(self, *args, **kwargs)
        self.setWindowTitle(title)
        self.generalLayout = QHBoxLayout()
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.centralWidget.setLayout(self.generalLayout)

        # self._createMenu()
        # self._createToolBar()
        self._createSideBar()
        self._createMplCanvas()
        self._createStatusBar()

    def tseries(self):
        pass


def detect_roms_file(filepath):
    for key in RomsNCFiles.__dataclass_fields__.keys():
        if key in basename(filepath):
            return key


def numpydatetime2str(numpydatetime):
    return str(numpydatetime).split(".")[0].replace("T", " ")


def last2d(da):
    if da.ndim <= 2:
        return da

    slc = [0] * (da.ndim - 2)
    slc += [slice(None), slice(None)]
    slc = {d: s for d, s in zip(da.dims, slc)}

    return da.isel(**slc)


def not_found_dialog(message="Coming soon..."):
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setInformativeText(message)
    msg.setWindowTitle("Not found")
    msg.exec_()


# Client code
def main():
    # Create an instance of QApplication
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    # Now use a palette to switch to dark colors:
    # palette = QPalette()
    # palette.setColor(QPalette.Window, QColor(53, 53, 53))
    # palette.setColor(QPalette.WindowText, Qt.white)
    # palette.setColor(QPalette.Base, QColor(25, 25, 25))
    # palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    # palette.setColor(QPalette.ToolTipBase, Qt.white)
    # palette.setColor(QPalette.ToolTipText, Qt.white)
    # palette.setColor(QPalette.Text, Qt.white)
    # palette.setColor(QPalette.Button, QColor(53, 53, 53))
    # palette.setColor(QPalette.ButtonText, Qt.white)
    # palette.setColor(QPalette.BrightText, Qt.red)
    # palette.setColor(QPalette.Link, QColor(42, 130, 218))
    # palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    # palette.setColor(QPalette.HighlightedText, Qt.black)
    # app.setPalette(palette)

    # Show the UI
    view = Ui()
    view.setGeometry(2500, 60, 1000, 800)
    view.show()

    if len(sys.argv) > 1:
        view.onOpenFile(sys.argv[1])

    # Create instances of the model and the controller
    # model = evaluateExpression
    # Controller(model=model, view=view)
    # Execute the calculator's main loop
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
