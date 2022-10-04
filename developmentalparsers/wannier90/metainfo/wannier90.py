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
import numpy as np

from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference, MEnum, JSON
)
from nomad.datamodel.metainfo import simulation
from nomad.datamodel.metainfo.simulation.method import KMesh


m_package = Package()


class Projections(MSection):
    '''
    Section containing the various parameters that define a Wannier90 projection
    '''

    m_def = Section(validate=False)

    k_mesh = SubSection(sub_section=KMesh.m_def, repeats=False)

    number_of_projected_orbitals = Quantity(
        type=np.dtype(np.int32),
        description='''
        Number of Wannier orbitals used to fit the DFT band structure
        ''')

    number_of_bands = Quantity(
        type=np.dtype(np.int32),
        description='''
        Number of input Bloch bands.
        ''')

    is_maximally_localise = Quantity(
        type=bool,
        description='''
        Are the projected orbitals maximally localized or just a single-shot projection?
        ''')

    convergence_tolerance_ml = Quantity(
        type=np.dtype(np.float64),
        description='''
        Convergence tolerance for maximal localization of the projected orbitals.
        ''')

    outer_energy_window = Quantity(
        type=np.dtype(np.float64),
        unit='electron_volt',
        shape=[2],
        description='''
        Bottom and top of the outer energy window used for the projection.
        ''')

    inner_energy_window = Quantity(
        type=np.dtype(np.float64),
        unit='electron_volt',
        shape=[2],
        description='''
        Bottom and top of the inner energy window used for the projection.
        ''')


class Method(simulation.method.Method):
    '''
    Section containing the various parameters that define the theory and the
    approximations (convergence, thresholds, etc.) behind the calculation.
    '''

    m_def = Section(validate=False, extends_base_section=True)

    projection = SubSection(sub_section=Projections.m_def, repeats=False)


class x_wannier90_hopping_parameters(MSection):
    '''
    Section containing the Wannier90 hopping parameters
    '''

    m_def = Section(validate=False)

    nrpts = Quantity(
        type=np.dtype(np.int32),
        description='''
        Number of Wigner-Seitz real points.
        ''')

    degeneracy_factors = Quantity(
        type=np.dtype(np.int32),
        shape=['nrpts'],
        description='''
        Degeneracy of each Wigner-Seitz grid point.
        ''')

    hopping_matrix = Quantity(
        type=np.dtype(np.float64),
        shape=['nrpts', 7],
        description='''
        Real space hopping matrix for each Wigner-Seitz grid point. The elements are
        defined as follows:
            n_x   n_y   n_z   orb_1   orb_2   real_part   imag_part
        where (n_x, n_y, n_z) define the Wigner-Seitz cell vector in fractional coordinates,
        (orb_1, orb_2) indicates the hopping amplitude between orb_1 and orb_2, and the
        real and imaginary parts of the hopping in electron_volt.
        ''')


class Run(simulation.run.Run):

    m_def = Section(validate=False, extends_base_section=True)

    x_wannier90_hoppings = SubSection(sub_section=x_wannier90_hopping_parameters, repeats=False)
