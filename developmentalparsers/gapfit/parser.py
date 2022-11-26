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

import logging
import numpy as np

from nomad.units import ureg
from nomad.parsing.file_parser import TextParser, Quantity
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.system import System, Atoms
from nomad.datamodel.metainfo.simulation.calculation import Calculation, Energy, EnergyEntry
from developmentalparsers.gapfit.metainfo.gapfit import (
    x_gapfit_gap_params, x_gapfit_gap_data, x_gapfit_e0, x_gapfit_gpsparse, x_gapfit_gp_coordinates,
    x_gapfit_sparsex)

import xml.etree.ElementTree as ET

re_f = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'

class XYZParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity(
                'frame',
                rf'(Lattice[\s\S]+?)(?:\s+ *\d+ *[\n\r]|\Z)', repeats=True,
                sub_parser=TextParser(quantities=[
                    Quantity(
                        'lattice',
                        r'Lattice=\"(.+?)\"', dtype=np.dtype(np.float64), shape=[3, 3]
                    ),
                    Quantity(
                        'labels_positions',
                        rf'([A-Z][a-z]*\s+{re_f}\s+{re_f}\s+{re_f})',
                        repeats=True
                    ),
                    Quantity('pbc', r'pbc=\"(.+?)\"'),
                    Quantity('energy', rf'energy=({re_f})', dtype=float)
                ])
            )
        ]


class GAPFitParser:
    def __init__(self) -> None:
        self.xyz_parser = XYZParser()

    def parse(self, filepath, archive, logger=None):
        sec_run = archive.m_create(Run)
        sec_run.program = Program(name='gapfit')
        sec_current = None

        logger = logger if logger is not None else logging.getLogger(__name__)

        for event, elem in ET.iterparse(filepath, events=('start', 'end')):
            if event == 'start':
                if elem.tag == 'GAP_params':
                    sec_params = sec_run.m_create(x_gapfit_gap_params)
                    sec_current = sec_params
                elif elem.tag == 'GAP_data':
                    sec_current = sec_current.m_create(x_gapfit_gap_data)
                elif elem.tag == 'gpSparse':
                    sec_current = sec_params.m_create(x_gapfit_gpsparse)
                elif elem.tag == 'gpCoordinates':
                    sec_current = sec_current.m_create(x_gapfit_gp_coordinates)
            if event == 'end':
                attrib = elem.attrib
                if elem.tag == 'Potential':
                    sec_run.x_gapfit_potential_label = attrib.get('label')
                    sec_run.x_gapfit_potential_init_args = attrib.get('init_args')
                elif elem.tag == 'GAP_params':
                    sec_current.x_gapfit_label = attrib.get('label')
                    sec_current.x_gapfit_gap_version = attrib.get('gap_version')
                elif elem.tag == 'GAP_data':
                    sec_current.do_core = attrib.get('do_core', '').lower().startswith('t')
                    sec_current = sec_run.x_gapfit_gap_params
                elif elem.tag == 'e0':
                    sec_e0 = sec_current.m_create(x_gapfit_e0)
                    sec_e0.x_gapfit_Z = int(attrib.get('Z', 0))
                    sec_e0.x_gapfit_value = float(attrib.get('value', 0))
                elif elem.tag == 'gpSparse':
                    sec_current.x_gapfit_label = attrib.get('label')
                    sec_current.x_gapfit_n_coordinate = int(attrib.get('n_coordinate', 0))
                    sec_current.x_gapfit_fitted = attrib.get('fitted', '').lower().startswith('t')
                    sec_current = sec_run.x_gapfit_gap_params
                elif elem.tag == 'gpCoordinates':
                    sec_current.x_gapfit_label = attrib.get('label')
                    sec_current.x_gapfit_dimensions = int(attrib.get('dimensions', 0))
                    sec_current.x_gapfit_signal_variance = float(attrib.get('signal_variance', 0))
                    sec_current.x_gapfit_signal_mean = float(attrib.get('signal_mean', 0))
                    sec_current.x_gapfit_sparsified = attrib.get('sparsified', '').lower().startswith('t')
                    sec_current.x_gapfit_n_permutations = int(attrib.get('n_permutations', 0))
                    sec_current.x_gapfit_covariance_type = int(attrib.get('covariance_type', 0))
                    sec_current.x_gapfit_zeta = float(attrib.get('zeta', 0))
                    sec_current.x_gapfit_n_sparseX = int(attrib.get('n_sparseX', 0))
                    sec_current.x_gapfit_sparseX_filename = attrib.get('sparseX_filename')
                    sec_current.x_gapfit_sparseX_md5sum = attrib.get('sparseX_md5sum')
                    sec_current = sec_run.x_gapfit_gap_params.x_gapfit_gpsparse[-1]
                elif elem.tag == 'descriptor':
                    sec_current.x_gapfit_descriptor = elem.text.strip()
                elif elem.tag == 'theta':
                    sec_current.x_gapfit_theta = float(elem.text)
                elif elem.tag == 'permutation':
                    sec_current.x_gapfit_permutation = int(elem.text)
                elif elem.tag == 'sparseX':
                    sec_sparsex = sec_current.m_create(x_gapfit_sparsex)
                    sec_sparsex.x_gapfit_i = int(attrib.get('i'))
                    sec_sparsex.x_gapfit_alpha = float(attrib.get('alpha'))
                    sec_sparsex.x_gapfit_sparseCutoff == float(attrib.get('sparseCutoff'))
                elif elem.tag == 'XYZ_data':
                    # TODO this is a hack to the text parser since, provide support for StringiO
                    self.xyz_parser._file_handler = elem.text.encode()
                    self.xyz_parser.parse()
                    for frame in self.xyz_parser.results.get('frame', []):
                        sec_system = sec_run.m_create(System)
                        atoms = Atoms()
                        lattice = frame.results.get('lattice')
                        if lattice is not None:
                            atoms.lattice_vectors = lattice * ureg.angstrom
                        labels_positions = frame.results.get('labels_positions')
                        if labels_positions is not None:
                            atoms.labels = [d[0] for d in labels_positions]
                            atoms.positions = np.array([d[1:4] for d in labels_positions], np.float64) * ureg.angstrom
                        pbc = frame.results.get('pbc')
                        if pbc is not None:
                            atoms.periodic = [d.lower().startswith('t') for d in pbc]
                        sec_system.atoms = atoms
                        energy = frame.results.get('energy')
                        if energy is not None:
                            sec_calc = sec_run.m_create(Calculation)
                            sec_calc.energy = Energy(total=EnergyEntry(value=energy * ureg.eV))
