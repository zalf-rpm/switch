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

import capnp
from datetime import date, timedelta, datetime
import json
from netCDF4 import Dataset
import numpy as np
import os
from pathlib import Path
import sys
import time
import zmq

import monica_run_lib
import shared

PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO / "../mas-infrastructure"
PATH_TO_PYTHON_CODE = PATH_TO_MAS_INFRASTRUCTURE_REPO / "src/python"
if str(PATH_TO_PYTHON_CODE) not in sys.path:
    sys.path.insert(1, str(PATH_TO_PYTHON_CODE))

from pkgs.common import common
from pkgs.model import monica_io3

PATH_TO_CAPNP_SCHEMAS = (PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas").resolve()
abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
fbp_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "fbp.capnp"), imports=abs_imports)

PATHS = {
    # adjust the local path to your environment
    "mbm-local-local": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        # "path-to-soil-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        "path-to-soil-dir": "/home/berg/Desktop/soil/",
        "monica-path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "mbm-local-remote": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        # "path-to-soil-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        "path-to-soil-dir": "/home/berg/Desktop/soil/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
    "hpc-local-remote": {
        #"path-to-climate-dir": "/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        # mounted path to archive or hard drive with climate data
        "path-to-soil-dir": "/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
        "path-debug-write-folder": "./debug-out/",
    },
}


def run_producer(server=None, port=None):
    context = zmq.Context()
    socket = context.socket(zmq.PUSH)  # pylint: disable=no-member

    config = {
        "mode": "mbm-local-remote",
        "server-port": port if port else "6666",  # local: 6667, remote 6666
        "server": server if server else "login01.cluster.zalf.de",
        "start_lat": "83.95833588",
        "end_lat": "-55.95833206",
        "start_lon": "-179.95832825",
        "end_lon": "179.50000000",
        "region": "africa",
        "resolution": "5min",  # 30sec,
        "sim.json": "sim.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups_africa_calibration.csv",
        "run-setups": "[1]",
        "reader_sr": None,
        "test_mode": "false",
        "path_to_out": "out/",
        "only_country_ids": "[]",  # "[10]",
    }

    common.update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    path_to_out_file = config["path_to_out"] + "/producer.out"
    if not os.path.exists(config["path_to_out"]):
        try:
            os.makedirs(config["path_to_out"])
        except OSError:
            print("run-calibration-producer.py: Couldn't create dir:", config["path_to_out"], "!")
    with open(path_to_out_file, "a") as _:
        _.write(f"config: {config}\n")

    only_country_ids = json.loads(config["only_country_ids"])

    s_resolution = {"5min": 5 / 60., "30sec": 30 / 3600.}[config["resolution"]]
    s_res_scale_factor = {"5min": 60., "30sec": 3600.}[config["resolution"]]

    region_to_lat_lon_bounds = {
        "nigeria": {"tl": {"lat": 14.0, "lon": 2.7}, "br": {"lat": 4.25, "lon": 14.7}},
        "africa": {"tl": {"lat": 37.4, "lon": -17.55}, "br": {"lat": -34.9, "lon": 51.5}},
        "earth": {
            "5min": {"tl": {"lat": 83.95833588, "lon": -179.95832825},
                     "br": {"lat": -55.95833206, "lon": 179.50000000}},
            "30sec": {"tl": {"lat": 83.99578094, "lon": -179.99583435},
                      "br": {"lat": -55.99583435, "lon": 179.99568176}}
        }
    }

    # select paths
    paths = PATHS[config["mode"]]
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["server-port"]))

    # read setup from csv file
    setups = monica_run_lib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"])

    # transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    # soil_crs_to_x_transformers = {}
    # wgs84_crs = CRS.from_epsg(4326)
    # utm32_crs = CRS.from_epsg(25832)
    # transformers[wgs84] = Transformer.from_crs(wgs84_crs, gk5_crs, always_xy=True)

    # open netcdfs
    path_to_soil_netcdfs = paths["path-to-soil-dir"] + "/" + config["resolution"] + "/"
    if config["resolution"] == "5min":
        soil_data = {
            "sand": {"var": "SAND", "file": "SAND5min.nc", "conv_factor": 0.01},  # % -> fraction
            "clay": {"var": "CLAY", "file": "CLAY5min.nc", "conv_factor": 0.01},  # % -> fraction
            "corg": {"var": "OC", "file": "OC5min.nc", "conv_factor": 0.01},  # scale factor
            "bd": {"var": "BD", "file": "BD5min.nc", "conv_factor": 0.01 * 1000.0},  # scale factor * 1 g/cm3 = 1000 kg/m3
        }
    else:
        soil_data = None  # ["Sand5min.nc", "Clay5min.nc", "OC5min.nc", "BD5min.nc"]
    soil_datasets = {}
    soil_vars = {}
    for elem, data in soil_data.items():
        ds = Dataset(path_to_soil_netcdfs + data["file"], "r", format="NETCDF4")
        soil_datasets[elem] = ds
        soil_vars[elem] = ds.variables[data["var"]]

    def create_soil_profile(row, col):
        # skip first 4.5cm layer and just use 7 layers
        layers = []

        layer_depth = 8
        # find the fill value for the soil data
        for elem2 in soil_data.keys():
            for i in range(8):
                if np.ma.is_masked(soil_vars[elem2][i, row, col]):
                    if i < layer_depth:
                        layer_depth = i
                    break
        layer_depth -= 1

        if layer_depth < 4:
            return None
        
        for i, real_depth_cm, monica_depth_m in [(0, 4.5, 0), (1, 9.1, 0.1), (2, 16.6, 0.1), (3, 28.9, 0.1),
                                                 (4, 49.3, 0.2), (5, 82.9, 0.3), (6, 138.3, 0.6), (7, 229.6, 0.7)][1:]:
            if i <= layer_depth:
                layers.append({
                    "Thickness": [monica_depth_m, "m"],
                    "SoilOrganicCarbon": [soil_vars["corg"][i, row, col] * soil_data["corg"]["conv_factor"], "%"],
                    "SoilBulkDensity": [soil_vars["bd"][i, row, col] * soil_data["bd"]["conv_factor"], "kg m-3"],
                    "Sand": [soil_vars["sand"][i, row, col] * soil_data["sand"]["conv_factor"], "fraction"],
                    "Clay": [soil_vars["clay"][i, row, col] * soil_data["clay"]["conv_factor"], "fraction"]
                })
        return layers

    sent_env_count = 0
    start_time = time.perf_counter()

    if len(run_setups) > 1 and run_setups[0] not in setups:
        return
    else:
        setup_id = run_setups[0]

    conman = common.ConnectionManager()
    reader = conman.try_connect(config["reader_sr"], cast_as=fbp_capnp.Channel.Reader, retry_secs=1)
    if reader:
        while True:
            msg = reader.read().wait()
            # check for end of data from in port
            if msg.which() == "done":
                break

            env_template = None
            start_setup_time = None
            try:
                in_ip = msg.value.as_struct(fbp_capnp.IP)
                s: str = in_ip.content.as_text()
                params = json.loads(s)  # keys: MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
                if "only_country_ids" in params:
                    only_country_ids = params["only_country_ids"]
                    del params["only_country_ids"]

                start_setup_time = time.perf_counter()

                setup = setups[setup_id]
                gcm = setup["gcm"]
                scenario = setup["scenario"]
                ensmem = setup["ensmem"]
                crop = setup["crop"]

                region = setup["region"] if "region" in setup else config["region"]
                lat_lon_bounds = region_to_lat_lon_bounds.get(region, {
                    "tl": {"lat": float(config["start_lat"]), "lon": float(config["start_lon"])},
                    "br": {"lat": float(config["end_lat"]), "lon": float(config["end_lon"])}
                })

                if setup["region"] == "nigeria":
                    planting = setup["planting"].lower()
                    nitrogen = setup["nitrogen"].lower()
                    management_file = f"{planting}_planting_{nitrogen}_nitrogen.csv"
                    # load management data
                    management = monica_run_lib.read_csv(paths["path-to-data-dir"] +
                                                  "/agro_ecological_regions_nigeria/" + management_file, key="id")
                else:
                    planting = nitrogen = management = None

                eco_data = shared.load_grid_cached(
                    paths["path-to-data-dir"] +
                    "/agro_ecological_regions_nigeria/agro-eco-regions_0.038deg_4326_wgs84_nigeria.asc", int)
                country_id_data = shared.load_grid_cached(
                    paths["path-to-data-dir"] + "country-id_0.083deg_4326_wgs84_africa.asc", int)
                crop_mask_data = shared.load_grid_cached(
                    paths["path-to-data-dir"] + f"{setup['crop']}-mask_0.083deg_4326_wgs84_africa.asc.gz", int)
                planting_data = shared.load_grid_cached(
                    paths["path-to-data-dir"] + f"{setup['crop']}-planting-doy_0.5deg_4326_wgs84_africa.asc", int)
                harvest_data = shared.load_grid_cached(
                    paths["path-to-data-dir"] + f"{setup['crop']}-harvest-doy_0.5deg_4326_wgs84_africa.asc", int)
                height_data = shared.load_grid_cached(setup["path_to_dem_asc_grid"], float)
                slope_data = shared.load_grid_cached(setup["path_to_slope_asc_grid"], float)

                # read template sim.json
                with open(setup.get("sim.json", config["sim.json"])) as _:
                    sim_json = json.load(_)
                # change start and end date acording to setup
                if setup["start_date"]:
                    sim_json["climate.csv-options"]["start-date"] = str(setup["start_date"])
                if setup["end_date"]:
                    end_year = int(setup["end_date"].split("-")[0])
                    sim_json["climate.csv-options"]["end-date"] = str(setup["end_date"])

                    # read template site.json
                with open(setup.get("site.json", config["site.json"])) as _:
                    site_json = json.load(_)

                if len(scenario) > 0 and scenario[:3].lower() == "ssp":
                    site_json["EnvironmentParameters"]["rcp"] = f"rcp{scenario[-2:]}"

                # read template crop.json
                with open(setup.get("crop.json", config["crop.json"])) as _:
                    crop_json = json.load(_)
                    # set current crop
                    for ws in crop_json["cropRotation"][0]["worksteps"]:
                        if "Sowing" in ws["type"]:
                            ws["crop"][2] = crop
                    # set value of calibration params
                    ps = crop_json["crops"][crop]["cropParams"]
                    for pname, pval in params.items():
                        if pname in ps["species"]:
                            ps["species"][pname] = pval
                        elif pname in ps["cultivar"]:
                            ps["cultivar"][pname] = pval

                crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = setup[
                    "use_vernalisation_fix"] if "use_vernalisation_fix" in setup else False

                # create environment template from json templates
                env_template = monica_io3.create_env_json_from_json_config({
                    "crop": crop_json,
                    "site": site_json,
                    "sim": sim_json,
                    "climate": ""
                })

                c_lon_0 = -179.75
                c_lat_0 = +89.25
                c_resolution = 0.5

                s_lat_0 = region_to_lat_lon_bounds["earth"][config["resolution"]]["tl"]["lat"]
                s_lon_0 = region_to_lat_lon_bounds["earth"][config["resolution"]]["tl"]["lon"]
                b_lat_0 = lat_lon_bounds["tl"]["lat"]
                b_lon_0 = lat_lon_bounds["tl"]["lon"]

                lats_scaled = range(int(lat_lon_bounds["tl"]["lat"] * s_res_scale_factor),
                                    int(lat_lon_bounds["br"]["lat"] * s_res_scale_factor) - 1,
                                    -int(s_resolution * s_res_scale_factor))
                no_of_lats = len(lats_scaled)
                s_row_0 = int((s_lat_0 - (lats_scaled[0] / s_res_scale_factor)) / s_resolution)
                for lat_scaled in lats_scaled:
                    lat = lat_scaled / s_res_scale_factor
                    #print("lat:"+str(round(lat,2)))
                    print(str(round(lat, 2)), end=" ", flush=True)

                    lons_scaled = range(int(lat_lon_bounds["tl"]["lon"] * s_res_scale_factor),
                                        int(lat_lon_bounds["br"]["lon"] * s_res_scale_factor) + 1,
                                        int(s_resolution * s_res_scale_factor))
                    no_of_lons = len(lons_scaled)
                    s_col_0 = int(((lons_scaled[0] / s_res_scale_factor) - s_lon_0) / s_resolution)
                    for lon_scaled in lons_scaled:
                        lon = lon_scaled / s_res_scale_factor
                        #print("lon:"+str(round(lon,2)), end=" ")
                        #print(".", end="", flush=True)

                        c_col = int((lon - c_lon_0) / c_resolution)
                        c_row = int((c_lat_0 - lat) / c_resolution)

                        s_col = int((lon - s_lon_0) / s_resolution)
                        s_row = int((s_lat_0 - lat) / s_resolution)

                        # set management
                        mgmt = None
                        aer = None
                        if setup["region"] == "nigeria":
                            aer = eco_data["value"](lat, lon, False)
                            if aer and aer > 0 and aer in management:
                                mgmt = management[aer]
                        else:
                            mgmt = {}
                            planting_doy = planting_data["value"](lat, lon, False)
                            if planting_doy:
                                d = date(2023, 1, 1) + timedelta(days=planting_doy-1)
                                mgmt["Sowing date"] = f"0000-{d.month:02}-{d.day:02}"
                            harvest_doy = harvest_data["value"](lat, lon, False)
                            if harvest_doy:
                                d = date(2023, 1, 1) + timedelta(days=harvest_doy - 1)
                                mgmt["Harvest date"] = f"0000-{d.month:02}-{d.day:02}"

                        valid_mgmt = False
                        if mgmt and shared.check_for_nill_dates(mgmt) and len(mgmt) > 1:
                            valid_mgmt = True
                            for ws in env_template["cropRotation"][0]["worksteps"]:
                                if ws["type"] == "Sowing" and "Sowing date" in mgmt:
                                    ws["date"] = shared.mgmt_date_to_rel_date(mgmt["Sowing date"])
                                    if "Planting density" in mgmt:
                                        ws["PlantDensity"] = [float(mgmt["Planting density"]), "plants/m2"]
                                elif ws["type"] == "Harvest" and "Harvest date" in mgmt:
                                    ws["date"] = shared.mgmt_date_to_rel_date(mgmt["Harvest date"])
                                elif ws["type"] == "AutomaticHarvest" and "Harvest date" in mgmt:
                                    ws["latest-date"] = shared.mgmt_date_to_rel_date(mgmt["Harvest date"])
                                elif ws["type"] == "Tillage" and "Tillage date" in mgmt:
                                    ws["date"] = shared.mgmt_date_to_rel_date(mgmt["Tillage date"])
                                elif ws["type"] == "MineralFertilization" and mgmt[:2] == "N " and mgmt[-5:] == " date":
                                    app_no = int(ws["application"])
                                    app_str = str(app_no) + ["st", "nd", "rd", "th"][app_no - 1]
                                    ws["date"] = shared.mgmt_date_to_rel_date(mgmt[f"N {app_str} date"])
                                    ws["amount"] = [float(mgmt[f"N {app_str} application (kg/ha)"]), "kg"]
                        else:
                            mgmt = None
                        if not mgmt or not valid_mgmt:
                            continue

                        crop_mask_value = crop_mask_data["value"](lat, lon, False)
                        if not crop_mask_value or crop_mask_value == 0:
                            continue

                        country_id = country_id_data["value"](lat, lon, False)
                        if not country_id or (len(only_country_ids) > 0 and country_id not in only_country_ids):
                            continue

                        height_nn = height_data["value"](lat, lon, False)
                        if not height_nn:
                            continue

                        slope = slope_data["value"](lat, lon, False)
                        if not slope:
                            slope = 0

                        soil_profile = create_soil_profile(s_row, s_col)
                        if not soil_profile or len(soil_profile) == 0:
                            continue

                        env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup[
                            "LeafExtensionModifier"]

                        env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                        if setup["elevation"]:
                            env_template["params"]["siteParameters"]["heightNN"] = height_nn

                        if setup["slope"]:
                            if setup["slope_unit"] == "degree":
                                s = slope / 90.0
                            else:
                                s = slope
                            env_template["params"]["siteParameters"]["slope"] = s

                        if setup["latitude"]:
                            env_template["params"]["siteParameters"]["Latitude"] = lat

                        if setup["FieldConditionModifier"]:
                            for ws in env_template["cropRotation"][0]["worksteps"]:
                                if "Sowing" in ws["type"]:
                                    if "|" in setup["FieldConditionModifier"] and aer and aer > 0:
                                        fcms = setup["FieldConditionModifier"].split("|")
                                        fcm = float(fcms[aer-1])
                                        if fcm > 0:
                                            ws["crop"]["cropParams"]["species"]["FieldConditionModifier"] = fcm
                                    else:
                                        ws["crop"]["cropParams"]["species"]["FieldConditionModifier"] = \
                                            setup["FieldConditionModifier"]

                        env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup[
                            "fertilization"]
                        env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]
                        env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
                        env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup[
                            "WaterDeficitResponseOn"]
                        env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup[
                            "EmergenceMoistureControlOn"]
                        env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup[
                            "EmergenceFloodingControlOn"]

                        env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
                        hist_sub_path = "isimip/3b_v1.1_CMIP6/csvs/{gcm}/historical/{ensmem}/row-{crow}/col-{ccol}.csv.gz".format(
                            gcm=gcm, ensmem=ensmem, crow=c_row, ccol=c_col)
                        sub_path = "isimip/3b_v1.1_CMIP6/csvs/{gcm}/{scenario}/{ensmem}/row-{crow}/col-{ccol}.csv.gz".format(
                            gcm=gcm, scenario=scenario, ensmem=ensmem, crow=c_row, ccol=c_col
                        )
                        if setup["incl_historical"] and scenario != "historical":
                            climate_data_paths = [
                                paths["monica-path-to-climate-dir"] + hist_sub_path,
                                paths["monica-path-to-climate-dir"] + sub_path
                            ]
                        else:
                            climate_data_paths = [paths["monica-path-to-climate-dir"] + sub_path]
                        env_template["pathToClimateCSV"] = climate_data_paths
                        #print("pathToClimateCSV:", env_template["pathToClimateCSV"])

                        env_template["customId"] = {
                            "setup_id": setup_id,
                            "lat": lat, "lon": lon,
                            "no_of_s_cols": no_of_lons, "no_of_s_rows": no_of_lats,
                            "env_id": sent_env_count+1,
                            "nodata": False,
                            "country_id": country_id,
                        }

                        socket.send_json(env_template)
                        #with open(path_to_out_file, "a") as _:
                        #    _.write(f"sent env {sent_env_count} customId: {env_template['customId']}\n")
                        #print("sent env ", sent_env_count, " customId: ", env_template["customId"])

                        sent_env_count += 1

                        if config["test_mode"] == "true" and sent_env_count == 100:
                            raise Exception("leave early for test")
            except Exception as e:
                with open(path_to_out_file, "a") as _:
                    _.write(f"raised exception: {e}\n")
                #print("Exception raised:", e)
                raise e
            #    pass

            # send a last message will be just forwarded by monica to signify last
            if env_template:
                env_template["pathToClimateCSV"] = ""
                env_template["customId"] = {
                    "no_of_sent_envs": sent_env_count,
                    "nodata": True
                }
                socket.send_json(env_template)

            stop_setup_time = time.perf_counter()
            print("Setup ", sent_env_count, " envs took ", (stop_setup_time - start_setup_time), " seconds")
            with open(path_to_out_file, "a") as _:
                _.write(f"{datetime.now()} Setup {sent_env_count} envs took {stop_setup_time - start_setup_time} seconds\n")
            sent_env_count = 0

    stop_time = time.perf_counter()

    # write summary of used json files
    try:
        print("sending ", (sent_env_count - 1), " envs took ", (stop_time - start_time), " seconds")
        print("exiting run_producer()")
    except Exception:
        raise


if __name__ == "__main__":
    run_producer()
