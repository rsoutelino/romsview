#!/usr/bin/env python
######################################################
# GUI to vizualize ROMS input/output files
# Sep 2021
# rsoutelino@gmail.com
######################################################
import os
import wx
import datetime as dt

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Navbar
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import scipy.io as sp
import netCDF4 as nc

from lib import *

# TO-DO LIST: ====================================================
#   - correct bug with date selection: somehow the times re-start
#       every 00z
#   - need to decide which x-axis to use, lon or lat
# ================================================================

# NICE TIP TO DEBUG THIS PROGRAM: ================================
#   - comment out app.MainLoop at the last line of this script
#   - ipython --gui=wx
#   - run pyromsgui.py
#   - trigger the events and check out the objects in the shell
# ================================================================


global currentDirectory
currentDirectory = os.getcwd()

PROJECT_DIR = os.path.abspath(os.path.dirname(__file__))
DEFAULT_VMIN = 0
DEFAULT_VMAX = 1.5
DEFAULT_CMAP = plt.cm.BrBG
DEFAULT_DEPTH_FOR_LAND = -50


class App(wx.App):
    def OnInit(self):
        self.frame = Interface("PyRomsGUI 0.1.0", size=(1024, 800))
        self.frame.Show()
        return True


class Interface(wx.Frame):
    def __init__(self, title=wx.EmptyString, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE,
                 *args, **kwargs):
        wx.Frame.__init__(self, None, -1, "PyRomsGUI 0.1.0", pos=pos,
                          size=size, style=style, *args, **kwargs)

        # Initializing toolbar
        self.toolbar = MainToolBar(self)

        # BASIC LAYOUT OF THE NESTED SIZERS ======================
        panel1 = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel.SetBackgroundColour("WHITE")

        # BOX 1 is the main sizer
        box1 = wx.BoxSizer(wx.HORIZONTAL)
        box1.Add(panel1, 1, wx.EXPAND)
        box1.Add(mplpanel, 4, wx.EXPAND)

        # BOX 2 is the inner sizer of the left big control panel
        box2 = wx.BoxSizer(wx.VERTICAL)

        # BOX 3 is the sizer of the right big parent panel(panel1), the one that will
        #    serve as base for two child panels which will hold
        #    the two matplotlib canvas's
        box3 = wx.BoxSizer(wx.VERTICAL)

        # panel 1 content ========================================
        variable = wx.StaticText(panel1, label="Variable")
        box2.Add(variable, proportion=0, flag=wx.CENTER)
        self.var_select = wx.ComboBox(panel1, value='Choose variable')
        box2.Add(self.var_select, proportion=0, flag=wx.CENTER)
        self.var_select.Bind(wx.EVT_COMBOBOX, self.toolbar.OnUpdateHslice)

        time = wx.StaticText(panel1, label="Time record")
        box2.Add(time, proportion=0, flag=wx.CENTER)
        self.time_select = wx.ComboBox(panel1, value='Choose time step')
        box2.Add(self.time_select, proportion=0, flag=wx.CENTER)
        self.time_select.Bind(wx.EVT_COMBOBOX, self.toolbar.OnUpdateHslice)

        # mplpanel content ========================================
        self.mplpanel = SimpleMPLCanvas(mplpanel)
        box3.Add(self.mplpanel.canvas, 1, flag=wx.CENTER)

        # FINAL LAYOUT CONFIGURATIONS ============================
        self.SetAutoLayout(True)
        panel1.SetSizer(box2)
        mplpanel.SetSizer(box3)

        self.SetSizer(box1)

        self.InitMenu()
        self.Layout()
        self.Centre()

    def InitMenu(self):
        menubar = wx.MenuBar()
        fileMenu = wx.Menu()
        fileMenu.Append(wx.ID_OPEN, u'&Open ROMS grid file')
        fileMenu.Append(wx.ID_OPEN, u'&Open coastline file')
        fileMenu.Append(wx.ID_SAVE, '&Save grid')
        fileMenu.AppendSeparator()

        qmi = wx.MenuItem(fileMenu, wx.ID_EXIT, '&Quit\tCtrl+W')
        opf = wx.MenuItem(fileMenu, wx.ID_OPEN, '&Open\tCtrl+O')
        opc = wx.MenuItem(fileMenu, wx.ID_OPEN, '&Open\tCtrl+O+C')
        svf = wx.MenuItem(fileMenu, wx.ID_SAVE, '&Save\tCtrl+S')
        fileMenu.AppendItem(qmi)
        # fileMenu.AppendItem(svf)

        self.Bind(wx.EVT_MENU, self.OnQuit, qmi)
        self.Bind(wx.EVT_MENU, self.toolbar.OnLoadFile, opf)
        self.Bind(wx.EVT_MENU, self.toolbar.OnLoadCoastline, opc)
        self.Bind(wx.EVT_MENU, self.toolbar.OnPlotVslice, svf)

        menubar.Append(fileMenu, u'&PyRomsGUI')
        self.SetMenuBar(menubar)

    def OnQuit(self, e):
        """Fecha o programa"""
        self.Close()
        self.Destroy()

    def OnCloseWindow(self, e):
        self.Destroy()


class SimpleMPLCanvas(object):
    """docstring for SimpleMPLCanvas"""

    def __init__(self, parent):
        super(SimpleMPLCanvas, self).__init__()
        self.parent = parent
        self.plot_properties()
        self.make_navbar()

    def make_navbar(self):
        self.navbar = Navbar(self.canvas)
        self.navbar.SetPosition(wx.Point(0, 0))  # this is not working !!

    def plot_properties(self):
        # Create matplotlib figure
        self.fig = Figure(facecolor='w', figsize=(12, 8))
        self.canvas = FigureCanvas(self.parent, -1, self.fig)

        self.ax = self.fig.add_subplot(111)
        # tit = self.ax1.set_title("ROMS mask_rho", fontsize=12, fontweight='bold')
        # tit.set_position([0.9, 1.05])


class MainToolBar(object):
    def __init__(self, parent):
        self.currentDirectory = os.getcwd()
        self.parent = parent
        self.toolbar = parent.CreateToolBar(style=1, id=1,
                                            name="Toolbar")
        self.tools_params = {
            'load_file': (load_bitmap('grid.png'), u"Load ROMS netcdf file",
                          "Load ocean_???.nc ROMS netcdf file"),
            'load_coastline': (load_bitmap('coast.png'), u"Load coastline",
                               "Load *.mat coastline file [lon / lat poligons]"),
            'plot_vslice': (load_bitmap('save.png'), u"Plot vertical slice",
                            "Plot vertical slice of some variable"),
            'settings': (load_bitmap('settings.png'), u"PyRomsGUI settings",
                         "PyRomsGUI configurations"),
            'quit': (load_bitmap('exit.png'), u"Quit",
                     "Quit PyRomsGUI"),
        }

        self.createTool(self.toolbar, self.tools_params['load_file'],
                        self.OnLoadFile)
        self.createTool(self.toolbar, self.tools_params['load_coastline'],
                        self.OnLoadCoastline)

        self.toolbar.AddSeparator()
        # from IPython import embed; embed()
        self.plot_vslice = self.createTool(self.toolbar,
                                           self.tools_params['plot_vslice'],
                                           self.OnPlotVslice)

        self.toolbar.AddSeparator()

        self.createTool(self.toolbar, self.tools_params['settings'],
                        self.OnSettings)
        self.createTool(self.toolbar, self.tools_params['quit'],
                        self.parent.OnQuit)

        self.toolbar.Realize()

    def createTool(self, parent, params, evt, isToggle=False):
        tool = parent.AddTool(wx.NewId(), 'a', params[0], shortHelp=params[1])
        self.parent.Bind(wx.EVT_TOOL, evt, id=tool.GetId())
        return tool

    def OnLoadFile(self, evt):
        openFileDialog = wx.FileDialog(self.parent, "Open roms netcdf file [*.nc]",
                                       "/ops/hindcast/roms/", " ",
                                       "netcdf files (*.nc)|*.nc",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...

        filename = openFileDialog.GetPath()
        self.ncfile = nc.Dataset(filename)

        # this function is intended to return relevant information on the file
        varlist, axeslist, time = taste_ncfile(self.ncfile)

        timelist = romsTime2string(time)
        app.frame.var_select.SetItems(varlist)
        app.frame.time_select.SetItems(timelist)
        app.frame.time_select.SetValue(timelist[0])

        # opening ROMS grid
        openFileDialog = wx.FileDialog(self.parent, "Open roms GRID netcdf file [*_grd.nc]",
                                       "/ops/hindcast/roms/", " ",
                                       "netcdf files (*.nc)|*.nc",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...

        grdname = openFileDialog.GetPath()
        self.grd = nc.Dataset(grdname)

        lon = self.grd.variables['lon_rho'][:]
        lat = self.grd.variables['lat_rho'][:]
        h = self.grd.variables['h'][:]

        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        self.pcolor = ax.pcolormesh(lon, lat, h, cmap=plt.cm.terrain_r)
        ax.set_xlim([lon.min(), lon.max()])
        ax.set_ylim([lat.min(), lat.max()])
        ax.set_aspect('equal')

        mplpanel.canvas.draw()

    def OnUpdateHslice(self, evt):
        # from IPython import embed; embed()
        varname = app.frame.var_select.GetValue()
        var = self.ncfile.variables[varname]
        dimensions = var.dimensions
        grid = dimensions[-1].split('_')[-1]
        lon = self.grd.variables['lon_'+grid][:]
        lat = self.grd.variables['lat_'+grid][:]

        # time index
        varlist, axeslist, time = taste_ncfile(self.ncfile)
        timestr = app.frame.time_select.GetValue()
        selected_time = string2romsTime(timestr, self.ncfile)
        # from IPython import embed; embed()
        tindex = np.where(time[:] == selected_time)[0][0]

        if len(dimensions) == 3:
            arr = var[tindex, ...]
        if len(dimensions) == 4:
            arr = var[tindex, -1, ...]

        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        ax.clear()
        ax.pcolormesh(lon, lat, arr, cmap=plt.cm.jet)
        ax.set_xlim([lon.min(), lon.max()])
        ax.set_ylim([lat.min(), lat.max()])
        ax.set_title("%s   %s" % (varname, timestr))
        ax.set_aspect('equal')

        mplpanel.canvas.draw()

    def OnLoadCoastline(self, evt):
        openFileDialog = wx.FileDialog(self.parent, "Open coastline file - MATLAB Seagrid-like format",
                                       "/home/rsoutelino/metocean/projects/mermaid", " ",
                                       "MAT files (*.mat)|*.mat",
                                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if openFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...

        filename = openFileDialog.GetPath()
        coast = sp.loadmat(filename)
        lon, lat = coast['lon'], coast['lat']

        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        ax.plot(lon, lat, 'k')

        try:
            ax.set_xlim([self.grd.lonr.min(), self.grd.lonr.max()])
            ax.set_ylim([self.grd.latr.min(), self.grd.latr.max()])
        except AttributeError:  # just in case a grid was not loaded before
            ax.set_xlim([np.nanmin(lon), np.nanmax(lon)])
            ax.set_ylim([np.nanmin(lat), np.nanmax(lat)])

        ax.set_aspect('equal')
        mplpanel.canvas.draw()

    def OnPlotVslice(self, evt):
        mplpanel = app.frame.mplpanel
        self.cid = mplpanel.canvas.mpl_connect(
            'button_press_event', self.vslice)

    def OnSettings(self, evt):
        pass

    def vslice(self, evt):
        if evt.inaxes != app.frame.mplpanel.ax:
            return
        mplpanel = app.frame.mplpanel
        ax = mplpanel.ax
        x, y = evt.xdata, evt.ydata
        button = evt.button
        p = ax.plot(x, y, 'wo', markeredgecolor='k')
        try:
            self.points.append(p)
            self.area.append((x, y))
        except AttributeError:
            self.points = [p]
            self.area = [(x, y)]

        if len(self.points) == 2:
            ax.plot([self.area[0][0], self.area[1][0]],
                    [self.area[0][1], self.area[1][1]], 'k')

            p1, p2 = self.area[0], self.area[1]

        mplpanel.canvas.draw()

        if len(self.points) == 2:

            # assigning relevant variables
            varname = app.frame.var_select.GetValue()
            var = self.ncfile.variables[varname]
            dimensions = var.dimensions
            grid = dimensions[-1].split('_')[-1]
            lon = self.grd.variables['lon_'+grid][:]
            lat = self.grd.variables['lat_'+grid][:]

            ts = self.ncfile.variables['theta_s'][:]
            tb = self.ncfile.variables['theta_b'][:]
            hc = self.ncfile.variables['hc'][:]
            nlev = var.shape[1]
            sc = (np.arange(1, nlev + 1) - nlev - 0.5) / nlev
            sigma = self.ncfile.variables['Cs_r'][:]

            dl = (np.gradient(lon)[1].mean() + np.gradient(lat)[0].mean()) / 2
            siz = int(np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2) / dl)
            xs = np.linspace(p1[0], p2[0], siz)
            ys = np.linspace(p1[1], p2[1], siz)

            # time index
            varlist, axeslist, time = taste_ncfile(self.ncfile)
            timestr = app.frame.time_select.GetValue()
            selected_time = string2romsTime(timestr, self.ncfile)
            tindex = np.where(time[:] == selected_time)[0][0]

            # getting nearest values
            hsec, zeta, vsec = [], [], []
            for ind in range(xs.size):
                line, col = near2d(lon, lat, xs[ind], ys[ind])
                vsec.append(var[tindex, :, line, col])
                hsec.append(self.grd.variables['h'][line, col])
                zeta.append(self.ncfile.variables['zeta'][tindex, line, col])

            vsec = np.array(vsec).transpose()
            hsec, zeta = np.array(hsec), np.array(zeta)
            xs = xs.reshape(1, xs.size).repeat(nlev, axis=0)
            ys = ys.reshape(1, ys.size).repeat(nlev, axis=0)
            zsec = get_zlev(hsec, sigma,  5, sc, ssh=zeta, Vtransform=2)

            xs = np.ma.masked_where(vsec > 1e20, xs)
            ys = np.ma.masked_where(vsec > 1e20, ys)
            zsec = np.ma.masked_where(vsec > 1e20, zsec)
            vsec = np.ma.masked_where(vsec > 1e20, vsec)

            self.vslice_dialog = VsliceDialog(app.frame, xs, ys, zsec, vsec)
            del self.points, self.area

        mplpanel.canvas.draw()


class VsliceDialog(wx.Dialog):
    def __init__(self, parent, xs, ys, zsec, vsec, *args, **kwargs):
        wx.Dialog.__init__(self, parent, -1, "VARIABLE Vertical Slice, TIMERECORD", pos=(0, 0),
                           size=(1200, 600), style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        self.xs, self.ys, self.zsec, self.vsec = xs, ys, zsec, vsec

        # BASIC LAYOUT OF THE NESTED SIZERS ======================
        panel1 = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel = wx.Panel(self, wx.ID_ANY, style=wx.SUNKEN_BORDER)
        mplpanel.SetBackgroundColour("WHITE")

        # BOX 1 is the main sizer
        box1 = wx.BoxSizer(wx.HORIZONTAL)
        box1.Add(panel1, 1, wx.EXPAND)
        box1.Add(mplpanel, 4, wx.EXPAND)

        # BOX 2 is the inner sizer of the left control panel
        box2 = wx.BoxSizer(wx.VERTICAL)
        # BOX 3 is the sizer of the panel1
        box3 = wx.BoxSizer(wx.VERTICAL)

        # panel 1 content ========================================
        plot_type = wx.StaticText(panel1, label="Plot type")
        box2.Add(plot_type, proportion=0, flag=wx.CENTER)
        self.plot_select = wx.ComboBox(panel1, value='scatter')
        box2.Add(self.plot_select, proportion=0, flag=wx.CENTER)
        self.plot_select.Bind(wx.EVT_COMBOBOX, self.OnUpdatePlot)
        self.plot_select.SetItems(['scatter',  'pcolormesh',
                                   'contourf', 'contour'])

        minmax = wx.StaticText(panel1, label="Range")
        box2.Add(minmax, proportion=0, flag=wx.CENTER)
        self.max = wx.TextCtrl(panel1, value=str(vsec.max()))
        self.min = wx.TextCtrl(panel1, value=str(vsec.min()))
        box2.Add(self.max, proportion=0, flag=wx.CENTER)
        box2.Add(self.min, proportion=0, flag=wx.CENTER)

        scale = wx.StaticText(panel1, label="Scatter scale")
        box2.Add(scale, proportion=0, flag=wx.CENTER)
        self.scatter_scale = wx.SpinCtrl(panel1, value='50')
        box2.Add(self.scatter_scale, proportion=0, flag=wx.CENTER)

        # mplpanel content ========================================
        self.mplpanel = SimpleMPLCanvas(mplpanel)
        box3.Add(self.mplpanel.canvas, 1, flag=wx.CENTER)

        ax = self.mplpanel.ax
        pl = ax.scatter(xs.ravel(), zsec.ravel(), s=50, c=vsec.ravel(),
                        edgecolors='none', cmap=plt.cm.jet)
        self.mplpanel.ax2 = self.mplpanel.fig.add_axes(
            [0.93, 0.15, 0.015, 0.7])
        ax2 = self.mplpanel.ax2
        cbar = self.mplpanel.fig.colorbar(pl, cax=ax2)

        ax.set_xlim([xs.min(), xs.max()])
        ax.set_ylim([zsec.min(), zsec.max()])

        self.mplpanel.canvas.draw()

        # FINAL LAYOUT CONFIGURATIONS ============================
        self.SetAutoLayout(True)
        panel1.SetSizer(box2)
        mplpanel.SetSizer(box3)

        self.SetSizer(box1)
        self.Show()

    def OnUpdatePlot(self, evt):
        xs, ys, zsec, vsec = self.xs, self.ys, self.zsec, self.vsec
        ax, ax2 = self.mplpanel.ax, self.mplpanel.ax2
        ax.clear()
        ax2.clear()

        vmin, vmax = float(self.min.GetValue()), float(self.max.GetValue())
        plot_type = self.plot_select.GetValue()
        sc = self.scatter_scale.GetValue()

        if plot_type == 'scatter':
            pl = ax.scatter(xs.ravel(), zsec.ravel(), s=sc, c=vsec.ravel(),
                            vmin=vmin, vmax=vmax, edgecolors='none', cmap=plt.cm.jet)
        elif plot_type == 'pcolormesh':
            pl = ax.pcolormesh(xs, zsec, vsec, vmin=vmin,
                               vmax=vmax, cmap=plt.cm.jet)
        elif plot_type == 'contourf':
            zsec = np.array(zsec)
            f = np.where(np.isnan(zsec) == True)
            zsec[f] = 0
            levs = np.linspace(vmin, vmax, 50)
            pl = ax.contourf(xs, zsec, vsec, levs, cmap=plt.cm.jet)
        elif plot_type == 'contour':
            zsec = np.array(zsec)
            f = np.where(np.isnan(zsec) == True)
            zsec[f] = 0
            levs = np.linspace(vmin, vmax, 50)
            pl = ax.contour(xs, zsec, vsec, levs, cmap=plt.cm.jet)

        ax.set_xlim([xs.min(), xs.max()])
        ax.set_ylim([zsec.min(), zsec.max()])
        cbar = self.mplpanel.fig.colorbar(pl, cax=ax2)

        self.mplpanel.canvas.draw()


def taste_ncfile(ncfile):
    try:
        if "history" in ncfile.type:
            filetype = 'his'
        elif 'restart' in ncfile.type:
            filetype = 'rst'
    except AttributeError:
        print "Not a standard ROMS file !"
        filetype = 'clim'  # old wrapper

    varlist = ROMSVARS[filetype]['variables']
    axeslist = ROMSVARS[filetype]['axes']

    for axes in axeslist:
        if 'time' in axes:
            try:
                time = ncfile.variables[axes]
            except KeyError:
                time = ncfile.variables['time']  # for non-default axes name
        else:
            pass

    return varlist, axeslist, time


def romsTime2string(nctime):
    """
    nctime  :  netCDF4 variable
    """
    timeunits = nctime.units
    units = timeunits.split(' ')[0]
    tstart = dt.datetime.strptime(timeunits.split(' ')[-2], "%Y-%m-%d")
    timelist = []
    for t in nctime[:]:
        if units == 'seconds':
            current = tstart + dt.timedelta(seconds=t)
        if units == 'days':
            current = tstart + dt.timedelta(seconds=t*86400)

        timelist.append(current.strftime("%Y-%m-%d  %H h"))

    return timelist


def string2romsTime(timelist, ncfile):
    if not isinstance(timelist, list):
        timelist = [timelist]

    varlist, axeslist, time = taste_ncfile(ncfile)
    timeunits = time.units
    units = timeunits.split(' ')[0]
    tstart = dt.datetime.strptime(timeunits.split(' ')[-2], "%Y-%m-%d")

    romstime = []
    for timestr in timelist:
        dttime = dt.datetime.strptime(timestr, "%Y-%m-%d  %H h")
        delta = dttime - tstart
        if units == 'seconds':
            current = delta.seconds
        if units == 'days':
            current = delta.days

        romstime.append(current)

    if len(romstime) == 1:
        return romstime[0]
    else:
        return romstime


def load_bitmap(filename, direc=None):
    """
    Load a bitmap file from the ./icons subdirectory. 
    The filename parameter should not
    contain any path information as this is determined automatically.

    Returns a wx.Bitmap object
    copied from matoplotlib resources
    """

    if not direc:
        basedir = os.path.join(PROJECT_DIR, 'icons')
    else:
        basedir = os.path.join(PROJECT_DIR, direc)

    bmpFilename = os.path.normpath(os.path.join(basedir, filename))
    if not os.path.exists(bmpFilename):
        raise IOError('Could not find bitmap file "%s"; dying' % bmpFilename)

    bmp = wx.Bitmap(bmpFilename)
    return bmp


if __name__ == "__main__":
    app = App(False)
    app.MainLoop()
