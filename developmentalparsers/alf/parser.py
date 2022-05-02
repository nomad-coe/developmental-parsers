#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
import logging
import h5py

from nomad.datamodel.metainfo.simulation.run import Run
from nomad.datamodel.metainfo.simulation.calculation import Calculation

from .metainfo.alf import x_alf_observable, x_alf_observable_data


class ALFParser:
    def __init__(self):
        pass

    def parse(self, filepath, archive, logger):
        self.filepath = os.path.abspath(filepath)
        self.archive = archive
        self.logger = logging.getLogger(__name__) if logger is None else logger

        try:
            data = h5py.File(self.filepath)
        except Exception:
            self.logger.error('Error opening h5 file.')
            data = None

        if data is None:
            return

        sec_run = archive.m_create(Run)
        sec_calc = sec_run.m_create(Calculation)

        observables = [
            'Den_eq', 'Den_tau', 'Ener_scal', 'Green_eq', 'Green_tau', 'Kin_scal',
            'Part_scal', 'Pot_scal', 'SpinT_eq', 'SpinXY_eq', 'SpinXY_tau', 'SpinZ_eq',
            'SpinZ_tau'
        ]

        def parse_obser(source, target):
            # observables have arbitrary shapes, need to iterate until we until 1-d
            if len(source.shape) == 1:
                target.x_alf_obser.append(x_alf_observable_data(x_afl_data=source))
                return
            for data_n in source:
                parse_obser(data_n, target)

        for observable in observables:
            sec_observable = x_alf_observable()
            data_observable = data[observable]
            if 'sign' in data_observable.keys():
                sec_observable.x_alf_sign = data_observable['sign']
            if 'obser' in data_observable.keys():
                sec_observable.x_alf_obser = data_observable['obser'][:]
                # TODO determine if we need to flatten data
                # parse_obser(data_observable['obser'][:], sec_observable)
            # TODO parse other properties

            setattr(sec_calc, f'x_alf_{observable.lower()}', sec_observable)
