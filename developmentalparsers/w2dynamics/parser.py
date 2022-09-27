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
import re
import numpy as np

from nomad.parsing.file_parser import TextParser, Quantity
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.calculation import Calculation
from nomad.datamodel.metainfo.simulation.method import Method
from nomad.datamodel.metainfo.simulation.system import System
from .metainfo.w2dynamics import (
    x_w2dynamics_axes, x_w2dynamics_quantities, x_w2dynamics_config_parameters,
    x_w2dynamics_config_atoms_parameters, x_w2dynamics_config_general_parameters,
    x_w2dynamics_config_qmc_parameters, DMFT
)

re_n = r'[\n\r]'


class LogParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = [
            Quantity(
                'program_version', r'Version\s*([0-9.]+)',
                dtype=str, flatten=False
            )
        ]


class W2DynamicsParser:
    def __init__(self):
        self._re_namesafe = re.compile(r'[^\w]')
        self.log_parser = LogParser()

        self._dmft_keys_mapping = {
            'hamiltonian': 'local_hamiltonian',
            'nat': 'number_of_atoms_per_unit_cell',
            'nd': 'number_of_correlated_bands',
            'totdens': 'number_of_correlated_electrons',
            'udd': 'U',
            'vdd': 'Up',
            'jdd': 'JH',
            'beta': 'inverse_temperature',
            'magnetism': 'magnetic_state',
            'ntau': 'number_of_tau',
            'niw': 'number_of_matsubara_freq',
            'dc': 'double_counting_mixing'
        }

        self._dataset_run_mapping = {
            '.axes': x_w2dynamics_axes,
            '.quantities': x_w2dynamics_quantities
        }

    def parse_program_version(self):
        # read program version from .log file if present
        log_files = [f for f in os.listdir(self.maindir) if f.endswith('.log')]
        if log_files:
            if len(log_files) > 1:
                self.logger.warn('Multiple logging files found.')

            self.log_parser.mainfile = os.path.join(self.maindir, log_files[0])

            return self.log_parser.get('program_version', None)

    def parse_dataset(self, source, target, include=[]):
        for key in source.keys():
            if include and key not in include:
                continue
            # resolve value from 'value'
            value = source[key]
            if isinstance(value, h5py.Group) and 'value' in value.keys():
                value = value['value']
            if not isinstance(value, h5py.Dataset):
                continue
            name = self._re_namesafe.sub('_', key)
            if value.shape:
                setattr(target, f'x_w2dynamics_{name}', value[:])
            # mu is a single value
            if key == 'mu':
                setattr(target, f'x_w2dynamics_{name}', value)

    def parse_method(self, data):
        sec_run = self.archive.run[-1]

        sec_method = sec_run.m_create(Method)

        # Parse Method.x_w2dynamics_config quantities
        sec_config = sec_method.m_create(x_w2dynamics_config_parameters)
        config_sections = [
            x_w2dynamics_config_atoms_parameters, x_w2dynamics_config_general_parameters,
            x_w2dynamics_config_qmc_parameters]
        for subsection in config_sections:
            sec_config_subsection = sec_config.m_create(subsection)
            for key in data.get('.config').attrs.keys():
                parameters = data.get('.config').attrs.get(key)
                keys_mod = (key.replace('-', '_')).split('.')
                # TODO parse correctly sec_config_atoms as a repeat subsection
                setattr(sec_config_subsection, f'x_w2dynamics_{keys_mod[-1]}', parameters)

        # DMFT section
        sec_dmft = sec_method.m_create(DMFT)
        for key in data.get('.config').attrs.keys():
            parameters = data.get('.config').attrs.get(key)
            keys_mod = (key.replace('-', '_')).split('.')[-1]

            # asign DMFT specific metadata
            if keys_mod in self._dmft_keys_mapping:
                setattr(sec_dmft, self._dmft_keys_mapping[keys_mod], parameters)

            # check if umatrix.dat file exists
            if keys_mod == 'umatrix':
                sec_dmft.is_U_matrix_file = True
                umat_files = [f for f in os.listdir(self.maindir) if f.startswith(keys_mod)]
                if not len(umat_files):
                    sec_dmft.is_U_matrix_file = False
        sec_dmft.impurity_solver = 'CT-HYB'

    def parse_system(self):
        sec_run = self.archive.run[-1]

        sec_system = sec_run.m_create(System)

        # TODO wait for data from Wurzburg

    def parse_scc(self, data):
        sec_run = self.archive.run[-1]

        # order calculations
        calc_keys = [key for key in data.keys() if key.startswith('dmft-') or key.startswith('stat-')]
        calc_keys.sort()

        # NOT converged results; these are stored in the dmft-last (or stat-last) group
        # calc_keys.insert(0, 'start')
        # calc_keys.insert(-1, 'finish')

        calc_quantities = [
            'dc-latt', 'gdensnew', 'gdensold', 'glocnew-lattice', 'glocold-lattice', 'mu']
        for key in calc_keys:
            if key not in data.keys():
                continue
            sec_calc = sec_run.m_create(Calculation)
            self.parse_dataset(data[key], sec_calc, include=calc_quantities)
            for sub_key in data[key].keys():
                if sub_key.startswith('ineq-'):
                    sec_ineq = sec_calc.m_create(x_w2dynamics_quantities)
                    self.parse_dataset(data[key][sub_key], sec_ineq)

    def parse(self, filepath, archive, logger):
        self.filepath = filepath
        self.archive = archive
        self.maindir = os.path.dirname(self.filepath)
        self.logger = logging.getLogger(__name__) if logger is None else logger

        self._method_type = 'DMFT'
        # TODO move mainfile as hdf5 to MatchingParserInterface list and
        # create init_parser() for the mainfile?
        try:
            data = h5py.File(self.filepath)
        except Exception:
            self.logger.error('Error opening hdf5 file.')
            data = None

        if data is None:
            return

        sec_run = archive.m_create(Run)

        # Program section
        sec_run.program = Program(
            name='w2dynamics', version=self.parse_program_version())

        # run.x_w2dynamics_axes section
        sec_axes = sec_run.m_create(x_w2dynamics_axes)
        self.parse_dataset(data.get('.axes'), sec_axes)

        # Method section
        self.parse_method(data)

        # Calculation section
        self.parse_scc(data)

        # System section
        self.parse_system()
