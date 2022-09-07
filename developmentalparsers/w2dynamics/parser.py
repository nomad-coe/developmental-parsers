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

from nomad.parsing.file_parser import TextParser, Quantity
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.calculation import Calculation
from nomad.datamodel.metainfo.simulation.method import Method
from nomad.datamodel.metainfo.simulation.system import System
from .metainfo.w2dynamics import (
    x_w2dynamics_quantities, x_w2dynamics_axes, x_w2dynamics_input_parameters,
    x_w2dynamics_atom_parameters, DMFT
)

re_n = r'[\n\r]'
k_B_eV = 8.617333262e-5 # Boltzmann constant in eV K^-1


class ParameterParser(TextParser):
    def init_quantities(self):
        section_quantities = [
            Quantity(
                'parameter', f'{re_n} *([\w-]+ *= *.+)',
                #repeats=True, str_operation=lambda x: x.split(' = ', 1))
                repeats=True, str_operation=lambda x: [v.strip() for v in x.split('=', 1)])
        ]

        self._quantities = [
            Quantity(
                'general', r'(General\][\s\S]+?)(?:\]|\Z)',
                sub_parser=TextParser(quantities=section_quantities)
            ),
            Quantity(
                'atom', r'(\[\d+\]\][\s\S]+?)(?:\]|\Z)', repeats=True,
                sub_parser=TextParser(quantities=section_quantities)
            ),
            Quantity(
                'qmc', r'(QMC\][\s\S]+?)(?:\]|\Z)',
                sub_parser=TextParser(quantities=section_quantities)
            )
        ]


class InfoParser(TextParser):
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
        self.parameter_parser = ParameterParser()
        self.info_parser = InfoParser()

        self._method_keys_mapping = {
            'NAt' : 'number_of_atoms_per_unit_cell',
            'totdens' : 'number_of_correlated_electrons',
            'beta' : 'temperature',
            'magnetism' : 'magnetic_state',
            'dc' : 'self_energy_mixing',
            'Hamiltonian' : 'local_hamiltonian',
            'Nd' : 'number_of_correlated_bands',
            'Udd' : 'U',
            'Vdd' : 'Up',
            'Jdd' : 'JH',
            'Ntau' : 'number_of_tau',
            'Niw' : 'number_of_matsubara_freq'
        }

        self._method_default = {
            'general' : {
                'NAt' : 1,
                'totdens' : 0.0,
                'beta' : 100.0,
                'magnetism' : 'para',
                'dc' : 'anisimov'
            },
            'atom': {
                'Hamiltonian' : 'Density',
                'Nd' : None,
                'Udd' : 0.0,
                'Vdd' : 0.0,
                'Jdd' : 0.0
            },
            'qmc': {
                'Ntau' : 1000,
                'Niw' : 2000
            }
        }

    def parse(self, filepath, archive, logger):
        self.filepath = os.path.abspath(filepath)
        self.archive = archive
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.maindir = os.path.dirname(self.filepath)

        try:
            data = h5py.File(self.filepath)
        except Exception:
            self.logger.error('Error opening hdf5 file.')
            data = None

        if data is None:
            return

        sec_run = archive.m_create(Run)

        def parse_quantities(source, target, include=[]):
            for key in source.keys():
                if include and key not in include:
                    continue
                # resolve value from 'value'
                value = source[key]
                if isinstance(value, h5py.Group) and 'value' in value.keys():
                    value = value['value']
                if not isinstance(value, h5py.Dataset):
                    continue
                if value.shape:
                    name = self._re_namesafe.sub('_', key)
                    setattr(target, f'x_w2dynamics_{name}', value[:])

        # order calculations
        calc_keys = [key for key in data.keys() if key.startswith('dmft-') or key.startswith('stat-')]
        calc_keys.sort()
        calc_keys.insert(0, 'start')
        calc_keys.insert(-1, 'finish')

        calc_quantities = ['dc-latt', 'gdensnew', 'gdensold', 'glocnew-lattice', 'glocold-lattice', 'mu']
        for key in calc_keys:
            if key not in data.keys():
                continue
            sec_calc = sec_run.m_create(Calculation)
            parse_quantities(data[key], sec_calc, include=calc_quantities)
            for sub_key in data[key].keys():
                if sub_key.startswith('ineq-'):
                    sec_ineq = sec_calc.m_create(x_w2dynamics_quantities)
                    parse_quantities(data[key][sub_key], sec_ineq)

        if '.quantities' in data.keys():
            sec_quantities = sec_run.m_create(x_w2dynamics_quantities)
            parse_quantities(data['.quantities'], sec_quantities)

        if '.axes' in data.keys():
            sec_axes = sec_run.m_create(x_w2dynamics_axes)
            parse_quantities(data['.axes'], sec_axes)

        # TODO determine if parameters/program version can be read from hdf5 file
        # read parameters from parameter file if present
        in_files = [f for f in os.listdir(self.maindir) if f.endswith('.in')]
        if in_files:
            if len(in_files) > 1:
                self.logger.warn('Multiple parameter files found.')

            self.parameter_parser.mainfile = os.path.join(self.maindir, in_files[0])

            sec_system = sec_run.m_create(System)
            sec_method = sec_run.m_create(Method)
            sec_dmft = sec_method.m_create(DMFT)

            sec_input_parameters = x_w2dynamics_input_parameters()

            # [QMC]
            for key in ['qmc']:
                parameters = {key: val for key, val in self.parameter_parser.get(key, {}).get('parameter', [])}
                if parameters:
                    setattr(sec_input_parameters, f'x_w2dynamics_{key}', parameters)

                    # Method(DMFT) quantities
                    for key_mapping in self._method_keys_mapping.keys():
                        if key_mapping in self._method_default[key]:
                            setattr(sec_dmft, self._method_keys_mapping[key_mapping], self._method_default[key][key_mapping])
                            if key_mapping in sec_input_parameters.x_w2dynamics_qmc.keys():
                                setattr(sec_dmft, self._method_keys_mapping[key_mapping], sec_input_parameters.x_w2dynamics_qmc[key_mapping])

            # [General]
            for key in ['general']:
                parameters = {key: val for key, val in self.parameter_parser.get(key, {}).get('parameter', [])}
                if parameters:
                    setattr(sec_input_parameters, f'x_w2dynamics_{key}', parameters)

                    if parameters['DOS']:
                        sec_system.dos = parameters['DOS']

                    # Method(DMFT) quantities
                    for key_mapping in self._method_keys_mapping.keys():
                        if key_mapping in self._method_default[key]:
                            setattr(sec_dmft, self._method_keys_mapping[key_mapping], self._method_default[key][key_mapping])
                            if key_mapping in sec_input_parameters.x_w2dynamics_general.keys():
                                setattr(sec_dmft, self._method_keys_mapping[key_mapping], sec_input_parameters.x_w2dynamics_general[key_mapping])

                    # TODO decide on storing temperature or inverse_temperature
                    # Storing the temperature in kelvin from beta = 1/(kB*T)
                    sec_dmft.temperature = 1.0/(k_B_eV*parameters['beta'])

            sec_dmft.x_w2dynamics_input_parameters = sec_input_parameters

            # [Atoms]
            for atom in self.parameter_parser.get('atom', []):
                parameters = {key: val for key, val in atom.get('parameter', [])}
                if parameters:
                    sec_atom = sec_dmft.m_create(x_w2dynamics_atom_parameters)
                    sec_atom.x_w2dynamics_atom = parameters

                    # Method(DMFT) quantities
                    for key_mapping in self._method_keys_mapping.keys():
                        if key_mapping in self._method_default['atom']:
                            setattr(sec_dmft, self._method_keys_mapping[key_mapping], self._method_default['atom'][key_mapping])
                            if key_mapping in sec_atom.x_w2dynamics_atom.keys():
                                setattr(sec_dmft, self._method_keys_mapping[key_mapping], sec_atom.x_w2dynamics_atom[key_mapping])

            # TODO check if this is always true (EDcheck?)
            if sec_system.dos == 'EDcheck':
                sec_dmft.impurity_solver = 'ED'
            else:
                sec_dmft.impurity_solver = 'CT-HYB'

        sec_program = sec_run.m_create(Program)
        sec_program.name = 'w2dynamics'
        # read program version from .log file if present
        log_files = [f for f in os.listdir(self.maindir) if f.endswith('.log')]
        if log_files:
            if len(log_files) > 1:
                self.logger.warn('Multiple logging files found.')

            self.info_parser.mainfile = os.path.join(self.maindir, log_files[0])

            version = self.info_parser.get('program_version', '')
            if version:
                sec_program.version = version
