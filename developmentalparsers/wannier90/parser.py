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
from nomad.datamodel.metainfo.simulation.method import Method, KMesh
from nomad.datamodel.metainfo.simulation.system import System
from .metainfo.wannier90 import Projections

re_n = r'[\n\r]'


class WOutParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        kmesh_quantities = [
            Quantity(
                'n_points', r'Total points[\s=]*(\d+)', dtype=int,
                repeats=False),
            Quantity(
                'k_points', r'\|[\s\d]*(-*\d.[^\|]+)', repeats=True)]

        disentangle_quantities = [
            Quantity(
                'outer', r'\|\s*Outer:\s*([-\d.]+)\s*\w*\s*([-\d.]+)\s*\((?P<__unit>\w+)\)',
                dtype=float, repeats=False),
            Quantity(
                'inner', r'\|\s*Inner:\s*([-\d.]+)\s*\w*\s*([-\d.]+)\s*\((?P<__unit>\w+)\)',
                dtype=float, repeats=False)]

        self._quantities = [
            Quantity(
                Program.version, r'\s*\|\s*Release\:\s*([\d\.]+)\s*', repeats=False),
            Quantity(
                'k_mesh', rf'\s*(K-POINT GRID[\s\S]+?)(?:-\s*MAIN)', repeats=False,
                sub_parser=TextParser(quantities=kmesh_quantities)),
            Quantity(
                'Nwannier', r'\|\s*Number of Wannier Functions\s*\:\s*(\d+)',
                repeats=False),
            Quantity(
                'Nband', r'\|\s*Number of input Bloch states\s*\:\s*(\d+)',
                repeats=False),
            Quantity(
                'Niter', r'\|\s*Total number of iterations\s*\:\s*(\d+)',
                repeats=False),
            Quantity(
                'conv_tol', r'\|\s*Convergence tolerence\s*\:\s*([\d.eE-]+)',
                repeats=False),
            Quantity(
                'energy_windows', rf'(\|\s*Energy\s*Windows\s*\|[\s\S]+?)(?:Number of target bands to extract:)',
                repeats=False, sub_parser=TextParser(quantities=disentangle_quantities))]


class BandParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = [
            Quantity(
                'kpoint_eigenvalue', rf'{re_n}([\s\S]+)(?:\Z)', repeats=False)]
                #sub_parser=TextParser(quantities=[Quantity(
                #    'value', r'  \d.\d*E[\-\+]\d*  ', dtype=float, repeats=True)]))]


class HrParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        pass


class Wannier90Parser:
    def __init__(self):
        self.wout_parser = WOutParser()
        self.band_parser = BandParser()
        #self.hr_parser = HrParser()

        self._input_projection_mapping = {
            'Nwannier': 'number_of_projected_orbitals',
            'Nband': 'number_of_bands',
            'conv_tol': 'convergence_tolerance_ml'
        }

    def parse_method(self):
        sec_run = self.archive.run[-1]
        sec_method = sec_run.m_create(Method)
        sec_proj = sec_method.m_create(Projections)

        # k_mesh section
        sec_k_mesh = sec_proj.m_create(KMesh)
        sec_k_mesh.n_points = self.wout_parser.get('k_mesh', None).n_points
        sec_k_mesh.points = self.wout_parser.get('k_mesh', None).k_points[::2]

        # Wannier90 section
        for key in self._input_projection_mapping.keys():
            setattr(sec_proj, self._input_projection_mapping[key], self.wout_parser.get(key))
        if self.wout_parser.get('Niter') is not None:
            sec_proj.is_maximally_localise = False
            if self.wout_parser.get('Niter')>1:
                sec_proj.is_maximally_localise = True
        sec_proj.outer_energy_window = self.wout_parser.get('energy_windows').outer
        sec_proj.inner_energy_window = self.wout_parser.get('energy_windows').inner

    def parse_bandstructure(self):
        band_files = [f for f in os.listdir(self.maindir) if f.endswith('_band.dat')]
        if not band_files:
            return

        print(self.band_parser.get('kpoint_eigenvalue'))

    def parse_scc(self):
        sec_run = self.archive.run[-1]
        sec_calc = sec_run.m_create(Calculation)

        # Hopping matrices
        #sec_calc.hopping_matrix = ...

        # Wannier band structure
        #self.parse_bandstructure()

    def parse(self, filepath, archive, logger):
        self.filepath = filepath
        self.archive = archive
        self.maindir = os.path.dirname(self.filepath)
        self.logger = logging.getLogger(__name__) if logger is None else logger

        sec_run = archive.m_create(Run)

        # TODO move mainfile as hdf5 to MatchingParserInterface list and
        # create init_parser() for the mainfile?
        wout_files = [f for f in os.listdir(self.maindir) if f.endswith('.wout')]
        if not wout_files:
            self.logger.error('Error finding the woutput file.')

        self.wout_parser.mainfile = os.path.join(self.maindir, wout_files[0])

        # Program section
        sec_run.program = Program(
            name='Wannier90', version=self.wout_parser.get('version', ''))
        # TODO TimeRun section

        # Method section
        self.parse_method()

        self.parse_bandstructure()