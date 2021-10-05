from dataclasses import dataclass, field

import numpy as np


@dataclass
class RomsNCFiles:
    grd: str
    nud: str
    ini: str
    bry: str
    clm: str
    his: str
    rst: str
    avg: str
    qks: str
    dia: str
    frc: str
    rivers: str
    tides: str


@dataclass
class AppState:
    ds: object = None
    da: object = None
    cbar: str = "viridis"
    filetype: str = "grd"
    clicked_points: list = field(default_factory=list)
    current_file: str = ""

    @property
    def current_slice(self):
        sliced = {}
        for dim, coord in self.da.coords.items():
            if coord.size == 1:
                sliced[dim] = coord.values

        return sliced

    @property
    def var(self):
        return self.da.name


# representative variable for each ROMS file
REP_VAR = RomsNCFiles(
    grd="h",
    nud="temp_NudgeCoef",
    ini="temp",
    bry="temp_west",
    clm="temp",
    his="temp",
    rst="temp",
    avg="temp",
    qks="temp",
    dia="undefined",
    frc="Pair",
    rivers="river_salt",
    tides="tide_Eamp",
)
