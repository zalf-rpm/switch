#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from netCDF4 import Dataset
import monica_run_lib
import numpy as np


def update_config(config, argv, print_config=False, allow_new_keys=False):
    if len(argv) > 1:
        for arg in argv[1:]:
            k, v = arg.split("=", maxsplit=1)
            if allow_new_keys or k in config:
                config[k] = v.lower() == "true" if v.lower() in ["true", "false"] else v
        if print_config:
            print(config)

def get_lat_0_lon_0_resolution_from_grid_metadata(metadata):
    lat_0 = float(metadata["yllcorner"]) \
                + (float(metadata["cellsize"]) * float(metadata["nrows"])) \
                - (float(metadata["cellsize"]) / 2.0)
    lon_0 = float(metadata["xllcorner"]) + (float(metadata["cellsize"]) / 2.0)
    resolution = float(metadata["cellsize"])
    return {"lat_0": lat_0, "lon_0": lon_0, "res": resolution}


def check_for_nill_dates(mgmt):
    for key, value in mgmt.items():
        if "date" in key and value == "Nill":
            return False
    return True


def mgmt_date_to_rel_date(mgmt_date):
    if mgmt_date[:5] == "0000-":
        return mgmt_date

    day_str, month_short_name = mgmt_date.split("-")
    month_str = "00"
    if month_short_name == "Jan":
        month_str = "01"
    elif month_short_name == "Feb":
        month_str = "02"
    elif month_short_name == "Mar":
        month_str = "03"
    elif month_short_name == "Apr":
        month_str = "04"
    elif month_short_name == "May":
        month_str = "05"
    elif month_short_name == "Jun":
        month_str = "06"
    elif month_short_name == "Jul":
        month_str = "07"
    elif month_short_name == "Aug":
        month_str = "08"
    elif month_short_name == "Sep":
        month_str = "09"
    elif month_short_name == "Oct":
        month_str = "10"
    elif month_short_name == "Nov":
        month_str = "11"
    elif month_short_name == "Dec":
        month_str = "12"

    return f"0000-{month_str}-{int(day_str):02}"


class GlobalSoilDataSet:
    """Global Soil Dataset for Earth System Modeling"""
    def __init__(self, path_to_soil_dir, resolution):
        # open netcdfs
        path_to_soil_netcdfs = path_to_soil_dir + "/" + resolution + "/"
        if resolution == "5min":
            self.soil_data = {
                "sand": {"var": "SAND", "file": "SAND5min.nc", "conv_factor": 0.01},  # % -> fraction
                "clay": {"var": "CLAY", "file": "CLAY5min.nc", "conv_factor": 0.01},  # % -> fraction
                "corg": {"var": "OC", "file": "OC5min.nc", "conv_factor": 0.01},  # scale factor
                "bd": {"var": "BD", "file": "BD5min.nc", "conv_factor": 0.01 * 1000.0},  # scale factor * 1 g/cm3 = 1000 kg/m3
            }
        else:
            self.soil_data = None  # ["Sand5min.nc", "Clay5min.nc", "OC5min.nc", "BD5min.nc"]
        self.soil_datasets = {}
        self.soil_vars = {}
        for elem, data in self.soil_data.items():
            ds = Dataset(path_to_soil_netcdfs + data["file"], "r", format="NETCDF4")
            self.soil_datasets[elem] = ds
            self.soil_vars[elem] = ds.variables[data["var"]]

    def create_soil_profile(self, row, col):
        # skip first 4.5cm layer and just use 7 layers
        layers = []

        layer_depth = 8
        # find the fill value for the soil data
        for elem2 in self.soil_data.keys():
            for i in range(8):
                if np.ma.is_masked(self.soil_vars[elem2][i, row, col]):
                    if i < layer_depth:
                        layer_depth = i
                    break
                    # return None
        layer_depth -= 1

        if layer_depth < 4:
            return None

        for i, real_depth_cm, monica_depth_m in [(0, 4.5, 0), (1, 9.1, 0.1), (2, 16.6, 0.1), (3, 28.9, 0.1),
                                                 (4, 49.3, 0.2), (5, 82.9, 0.3), (6, 138.3, 0.6), (7, 229.6, 0.7)][1:]:
            if i <= layer_depth:
                layers.append({
                    "Thickness": [monica_depth_m, "m"],
                    "SoilOrganicCarbon": [self.soil_vars["corg"][i, row, col] * self.soil_data["corg"]["conv_factor"], "%"],
                    "SoilBulkDensity": [self.soil_vars["bd"][i, row, col] * self.soil_data["bd"]["conv_factor"], "kg m-3"],
                    "Sand": [self.soil_vars["sand"][i, row, col] * self.soil_data["sand"]["conv_factor"], "fraction"],
                    "Clay": [self.soil_vars["clay"][i, row, col] * self.soil_data["clay"]["conv_factor"], "fraction"]
                })
        return layers


def load_grid_cached(path_to_grid, val_type, print_path=False):
    if not hasattr(load_grid_cached, "cache"):
        load_grid_cached.cache = {}

    if path_to_grid in load_grid_cached.cache:
        return load_grid_cached.cache[path_to_grid]

    md, _ = monica_run_lib.read_header(path_to_grid)
    grid = np.loadtxt(path_to_grid, dtype=type, skiprows=len(md))
    print("read: ", path_to_grid)
    ll0r = get_lat_0_lon_0_resolution_from_grid_metadata(md)

    def col(lon):
        return int((lon - ll0r["lon_0"]) / ll0r["res"])

    def row(lat):
        return int((ll0r["lat_0"] - lat) / ll0r["res"])

    def value(lat, lon, return_no_data=False):
        c = col(lon)
        r = row(lat)
        if 0 <= r < md["nrows"] and 0 <= c < md["ncols"]:
            val = val_type(grid[r, c])
            if val != md["nodata_value"] or return_no_data:
                return val
        return None

    cache_entry = {
        "metadata": md, "grid": grid, "ll0r": ll0r,
        "col": lambda lon: col(lon),
        "row": lambda lat: row(lat),
        "value": lambda lat, lon, ret_no_data: value(lat, lon, ret_no_data)
    }
    load_grid_cached.cache[path_to_grid] = cache_entry
    return cache_entry

