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
    dialogs = []

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
            "What to do when clicking horizontal slice point (s)")
        self.plotSelector.addItems(
            ["MPL Toolbar", "Tseries on click", "Vslice on click"])
        self.plotSelector.setDisabled(True)
        self.plotSelector.activated[str].connect(self.toggle_plot)
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
        self.cid = None  # initialize variable to reprensent mpl_connect
        self.init_plot()

        # Create toolbar, passing canvas as first parameter, parent (self, the MainWindow) as second.
        self.mpltoolbar = NavigationToolbar(self.mplcanvas, self)

        layout = QVBoxLayout()
        layout.addWidget(self.mpltoolbar)
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
        self._state.current_file = filename
        self._load_dataset(filename)
        self._state.filetype = detect_roms_file(filename)
        # getting a representative var based on settings.rep_var
        rep_var = getattr(REP_VAR, self._state.filetype)
        self._state.da = last2d(self._state.ds[rep_var])
        self.plot(var_changed=True)

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

    def plot(self, var_changed=False):
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

        self._state.clicked_points.append([evt.xdata, evt.ydata])
        self._state.vslice_ref.append(self.mplcanvas.axes.plot(evt.xdata, evt.ydata,
                                                               'wo', markeredgecolor='k', zorder=10))
        self.mplcanvas.draw()

        if 'Vslice' in self.plotSelector.currentText():
            if len(self._state.clicked_points) < 2:
                return
            else:
                self._state.vslice_ref.append(self.mplcanvas.axes.plot(
                    *pairs2lists(self._state.clicked_points), 'k', zorder=9))
                self.mplcanvas.draw()
                dialog = VsliceDialog(
                    parent=self, title="Vertical Slice")
                dialog.setGeometry(2000, 60, 900, 500)
                dialog.show()
                dialog.plot()
                self.dialogs.append(dialog)

        if 'Tseries' in self.plotSelector.currentText():
            dialog = TseriesDialog(parent=self, title="Time Series")
            dialog.setGeometry(2000, 60, 900, 500)
            dialog.show()
            dialog.plot()
            self.dialogs.append(dialog)

        return

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

    def toggle_plot(self, plot_type):
        if plot_type in ["Tseries on click", "Vslice on click"]:
            if not self.cid:
                self.cid = self.mplcanvas.mpl_connect(
                    'button_press_event', self.timeseries_or_vslice)
            self.mpltoolbar.setDisabled(True)
        else:
            self.mplcanvas.mpl_disconnect(self.cid)
            self.cid = None
            self.mpltoolbar.setEnabled(True)

    def toggle_var(self, var):
        _slice = {}
        # need to remove dimensions that don't exist in the new var
        # if that's the case (Ex, toggling from 3D to 2D var)
        for dim, val in self._state.current_slice.items():
            if dim in self._state.ds[var].dims:
                _slice[dim] = val

        self._state.da = last2d(self._state.ds[var].sel(**_slice))
        self.plot(var_changed=True)

    def toggle_time(self, timestamp):
        _slice = self._state.current_slice.copy()
        for key in _slice.keys():
            if "time" in key:
                _slice[key] = self.timeSelector.currentText()

                self._state.da = last2d(
                    self._state.ds[self._state.var].sel(**_slice))
                self.plot()
                break

    def toggle_lev(self, lev):
        _slice = self._state.current_slice.copy()
        for key in _slice.keys():
            if "s_rho" in key:
                _slice[key] = self.levSelector.currentText()

                self._state.da = last2d(
                    self._state.ds[self._state.var].sel(**_slice))
                self.plot()
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
    def __init__(self, parent=None, title='ROMSView dialog', *args, **kwargs):
        # super().__init__(*args, **kwargs)
        QDialog.__init__(self, *args, **kwargs)
        self.parent = parent
        self._state = AppState()
        self._state.ds = self.parent._state.ds
        self._state.filetype = self.parent._state.filetype
        self.setWindowTitle(title)
        self.generalLayout = QHBoxLayout()
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.centralWidget.setLayout(self.generalLayout)
        self._createSideBar()
        self.plotSelector.setDisabled(True)
        self.levSelector.setDisabled(True)
        self._createMplCanvas()
        self._createStatusBar()
        self.status.showMessage(
            f"Current file: {self.parent._state.current_file}")

    def plot(self, var_changed=False):
        self._reset_mpl_axes()
        if not var_changed:
            var = self.parent._state.var
        # getting time index from parent plot
        for key, val in self.parent._state.current_slice.items():
            if "time" in key:
                sel = {key: val}
                self._state.da = self._state.ds[var].sel(**sel)
                self._state.za = self._state.ds['z_rho'].sel(**sel)
                # TODO open a dialog to load grid if z_rho not available, and
                # compute zlevels (clm files or outputs that did not save z_rho)

        xsec, zsec, vsec, xaxis = self._extract_slice()

        self._plot = self.mplcanvas.axes.contourf(xsec, zsec, vsec, 20)
        self.mplcanvas.axes.set_ylabel('z_rho')
        self.mplcanvas.axes.set_xlabel(xaxis)
        self.mplcanvas.axes.figure.colorbar(self._plot)

        self.cbarSelector.setEnabled(True)
        self.rangeBox.setEnabled(True)
        self.plotSelector.setDisabled(True)

        if (
            hasattr(self._state, "vmin") and hasattr(self._state, "vmax")
        ) and not var_changed:
            self._plot.set_norm(
                mpl.colors.Normalize(self._state.vmin, self._state.vmax)
            )
        else:
            self._reset_range(np.nanmin(self._state.da),
                              np.nanmax(self._state.da))

        self.mplcanvas.draw()
        self._update_vars()
        self._update_times()

    def closeEvent(self, event):
        self.parent._state.clicked_points.clear()

        for plot in self.parent._state.vslice_ref:
            for pl in plot:
                pl.remove()

        self.parent.mplcanvas.draw()

        self.parent._state.vslice_ref.clear()

    def _extract_slice(self):
        p1 = self.parent._state.clicked_points[0]
        p2 = self.parent._state.clicked_points[1]
        xi_rho, eta_rho = np.meshgrid(self.parent._state.ds.xi_rho.values,
                                      self.parent._state.ds.eta_rho.values)
        xaxis = eta_or_xi(p1, p2)

        xs, ys = get_segment(xi_rho, eta_rho, p1, p2)

        # finding nearest depths and data values to segment points
        for idx in range(xs.size):
            line, col = near2d(xi_rho, eta_rho, xs[idx], ys[idx])
            try:
                vsec = np.hstack((vsec, self._state.da.isel(
                    eta_rho=line, xi_rho=col).values[:, None]))
                zsec = np.hstack((zsec, self._state.za.isel(
                    eta_rho=line, xi_rho=col).values[:, None]))
            except:
                vsec = self._state.da.isel(
                    eta_rho=line, xi_rho=col).values[:, None]
                zsec = self._state.za.isel(
                    eta_rho=line, xi_rho=col).values[:, None]

        if xaxis == 'xi_rho':
            xsec = np.atleast_2d(xs).repeat(self._state.da.s_rho.size, axis=0)
        else:
            xsec = np.atleast_2d(ys).repeat(self._state.da.s_rho.size, axis=0)

        return xsec, zsec, vsec, xaxis


class TseriesDialog(VsliceDialog):
    def __init__(self, parent=None, title='ROMSView dialog', *args, **kwargs):
        # super().__init__(*args, **kwargs)
        QDialog.__init__(self, *args, **kwargs)
        self.parent = parent
        self._set_state()
        self._set_da()
        self.setWindowTitle(title)
        self.generalLayout = QHBoxLayout()
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.centralWidget.setLayout(self.generalLayout)
        self._createSideBar()
        self.plotSelector.setDisabled(True)
        self.levSelector.setDisabled(True)
        self._createMplCanvas()
        self._createStatusBar()
        self.status.showMessage(
            f"Current file: {self.parent._state.current_file}")

    def _set_state(self):
        """Grabing some state variables from parent"""
        self._state = AppState()
        self._state.ds = self.parent._state.ds
        self._state.filetype = self.parent._state.filetype
        self._state.clicked_points = self.parent._state.clicked_points

    def _set_da(self):

        for key, val in self.parent._state.current_slice.items():
            if key == 's_rho':
                sel = {key: val}
                break
            else:
                sel = {}

        xi_rho, eta_rho = np.meshgrid(self._state.ds.xi_rho.values,
                                      self._state.ds.eta_rho.values)

        sel['eta_rho'], sel['xi_rho'] = near2d(
            xi_rho,
            eta_rho,
            self._state.clicked_points[0][0],
            self._state.clicked_points[0][1])

        self._state.da = self._state.ds[self.parent._state.var].sel(**sel)

    def plot(self, var_changed=False):
        self._reset_mpl_axes()
        if not var_changed:
            try:
                var = self._state.var
            except:
                var = self.parent._state.var

        self._plot = self._state.da.plot(ax=self.mplcanvas.axes)

        self.mplcanvas.draw()
        self._update_vars()
        self._update_times()

        if 's_rho' in self._state.da.coords.keys():
            self.levSelector.setEnabled(True)
            self._update_levels()
        else:
            self.levSelector.setDisabled(True)

    def closeEvent(self, event):
        self.parent._state.clicked_points.clear()

        for plot in self.parent._state.vslice_ref:
            for pl in plot:
                pl.remove()

        self.parent.mplcanvas.draw()

        self.parent._state.vslice_ref.clear()


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


def pairs2lists(pairs):
    """Transforms list of pairs of [(x1, y1), (x2, y2), ...]
       in lists of coords [x1, x2, ...], [y1, y2, ...]
    """
    x, y = [], []
    for pair in pairs:
        x.append(pair[0])
        y.append(pair[1])

    return x, y


def eta_or_xi(p1, p2):
    """Chooses preferred X axis based on alignment of two points that 
       define a transect
    """
    if np.abs(p1[0] - p2[0]) >= np.abs(p1[1] - p2[1]):
        return 'xi_rho'
    else:
        return 'eta_rho'


def near2d(x, y, x0, y0):
    """
    Find the indexes of the grid point that is
    nearest a chosen (x0, y0).
    Usage: line, col = near2d(x, y, x0, y0)
    """
    dx = np.abs(x - x0)
    dx = dx / dx.max()
    dy = np.abs(y - y0)
    dy = dy / dy.max()
    dn = dx + dy
    fn = np.where(dn == dn.min())
    line = int(fn[0])
    col = int(fn[1])
    return line, col


def get_segment(xg, yg, p1, p2):
    """Creates a segment from p1 to p2 based on average resolution of 
       xg, yg 2D generic grid

    Args:
        xg, yg [numpy 2D arrays]: x, y grid coordinates
        p1, p2 [list or tuple]: x, y pair of start and end of the segment 

    Returns:
        xs, ys [numpy 1D arrays]: x, y coordinates of the segment
    """
    dl = (np.gradient(xg)[1].mean() +
          np.gradient(yg)[0].mean()) / 2
    size = int(np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2) / dl)
    xs = np.linspace(p1[0], p2[0], size)
    ys = np.linspace(p1[1], p2[1], size)

    return xs, ys


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
