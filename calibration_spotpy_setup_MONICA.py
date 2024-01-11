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
from datetime import datetime
import json
import os
from pathlib import Path

import numpy as np
import spotpy
import re

PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO / "../mas-infrastructure"
PATH_TO_CAPNP_SCHEMAS = (PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas").resolve()
abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
fbp_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "fbp.capnp"), imports=abs_imports)

class spot_setup(object):
    def __init__(self, user_params, observations, prod_writer, cons_reader, path_to_out, only_nuts3_region_ids):
        self.user_params = user_params
        self.params = []
        self.observations = observations
        self.obs_flat_list = list(map(lambda d: d["value"], observations))
        self.prod_writer = prod_writer
        self.cons_reader = cons_reader
        self.path_to_out_file = path_to_out + "/spot_setup.out"
        self.only_nuts3_region_ids = only_nuts3_region_ids

        if not os.path.exists(path_to_out):
            try:
                os.makedirs(path_to_out)
            except OSError:
                print("spot_setup.__init__: Couldn't create dir:", path_to_out, "!")

        with open(self.path_to_out_file, "a") as _:
            _.write(f"observations: {self.observations}\n")
            _.write(f"obs_flat_list: {self.obs_flat_list}\n")

        for par in user_params:
            par_name = par["name"]
            if "array" in par:
                par["name"] = f"{par_name}_{par['array']}"  # spotpy does not allow two parameters to have the same name
                del par["array"]
            if "derive_function" not in par:  # spotpy does not care about derived params
                self.params.append(spotpy.parameter.Uniform(**par))

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # vector = MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
        msg_content = dict(zip(vector.name, vector))
        msg_content["only_nuts3_region_ids"] = self.only_nuts3_region_ids
        out_ip = fbp_capnp.IP.new_message(content=json.dumps(msg_content))
        self.prod_writer.write(value=out_ip).wait()
        with open(self.path_to_out_file, "a") as _:
            _.write(f"{datetime.now()} sent params to monica setup: {vector}\n")
        print("sent params to monica setup:", vector, flush=True)

        msg = self.cons_reader.read().wait()
        # check for end of data from in port
        if msg.which() == "done":
            return

        in_ip = msg.value.as_struct(fbp_capnp.IP)
        s: str = in_ip.content.as_text()
        nuts3_region_id_and_year_to_avg_yield = json.loads(s)
        # print("received monica results:", country_id_and_year_to_avg_yield, flush=True)

        # remove all simulation results which are not in the observed list
        sim_list = []
        for d in self.observations:
            key = f"{d['id']}|{d['year']}"
            if key in nuts3_region_id_and_year_to_avg_yield:
                if np.isnan(d["value"]):
                    sim_list.append(np.nan)
                else:
                    sim_list.append(nuts3_region_id_and_year_to_avg_yield[key])
            else:
                sim_list.append(np.nan)

        print("len(sim_list):", len(sim_list), "== len(self.obs_list):", len(self.obs_flat_list), flush=True)
        with open(self.path_to_out_file, "a") as _:
            #_.write(f"received monica results: {country_id_and_year_to_avg_yield}\n")
            _.write(f"{datetime.now()}  len(sim_list): {len(sim_list)} == len(self.obs_list): {len(self.obs_flat_list)}\n")
            #_.write(f"sim_list: {sim_list}\n")
            #_.write(f"obs_list: {self.obs_flat_list}\n")
        # besides the order the length of observation results and simulation results should be the same
        assert len(sim_list) == len(self.obs_flat_list)
        return sim_list if len(sim_list) > 0 else None

    def evaluation(self):
        return self.obs_flat_list

    def objectivefunction(self, simulation, evaluation):
        return spotpy.objectivefunctions.rmse(evaluation, simulation)
