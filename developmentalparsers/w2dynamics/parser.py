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
from nomad.datamodel.metainfo.simulation.run import Run
from nomad.datamodel.metainfo.simulation.calculation import Calculation
from nomad.datamodel.metainfo.simulation.method import Method
from .metainfo.w2dynamics import (
    x_w2dynamics_quantities, x_w2dynamics_axes, x_w2dynamics_input_parameters, 
    x_w2dynamics_atom_parameters
)

re_n = r'[\n\r]'


class ParameterParser(TextParser):
    def init_quantities(self):
        section_quantities = [
            Quantity(
                'parameter', f'{re_n} *(\w+ *= *.+)',
                repeats=True, str_operation=lambda x: x.split(' = ', 1))
        ]

        self._quantities = [
            Quantity(
                'general', r'(General\][\s\S]+?)(\]|\Z)',
                sub_parser=TextParser(quantities=section_quantities)
            ),
            Quantity(
                'atom', r'(\[\d+\]\][\s\S]+?)(\[|\Z)', repeats=True,
                sub_parser=TextParser(quantities=section_quantities)
            ),
            Quantity(
                'qmc', r'(QMC\][\s\S]+?)(\]|\Z)',
                sub_parser=TextParser(quantities=section_quantities)
            )
        ]


class W2DynamicsParser:
    def __init__(self):
        self._re_namesafe = re.compile(r'[^\w]')
        self.parameter_parser = ParameterParser()

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

        # read parameters from parameter file
        # TODO determine if parameters can be read from hdf5 file
        in_files = [f for f in os.listdir(self.maindir) if f.endswith('.in')]
        if in_files:
            if len(in_files) > 1:
                self.logger.warn('Multiple parameter files found.')

            self.parameter_parser.mainfile = os.path.join(self.maindir, in_files[0])

            sec_input_parameters = x_w2dynamics_input_parameters()
            for key in ['general', 'qmc']:
                parameters = {key: val for key, val in self.parameter_parser.get(key, {}).get('parameter', [])}
                if parameters:
                    setattr(sec_input_parameters, f'x_w2dynamics_{key}', parameters)

            sec_method = sec_run.m_create(Method)
            sec_method.x_w2dynamics_input_parameters = sec_input_parameters

            for atom in self.parameter_parser.get('atom', []):
                parameters = {key: val for key, val in atom.get('parameter', [])}
                if parameters:
                    sec_atom = sec_method.m_create(x_w2dynamics_atom_parameters)
                    sec_atom.x_w2dynamics_atom = parameters
