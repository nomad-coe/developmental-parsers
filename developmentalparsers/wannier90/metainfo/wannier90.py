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


m_package = Package()


class Projections(MSection):
    '''
    Section containing the various parameters that define a Wannier90 projection
    '''

    m_def = Section(validate=False)


class Method(simulation.method.Method):
    '''
    Section containing the various parameters that define the theory and the
    approximations (convergence, thresholds, etc.) behind the calculation.
    '''

    m_def = Section(validate=False, extends_base_section=True)

    projection = SubSection(sub_section=Projections.m_def, repeats=False)


class Calculation(simulation.calculation.Calculation):
    '''
    Every calculation section contains the values computed
    during a *single configuration calculation*, i.e. a calculation performed on a given
    configuration of the system (as defined in section_system) and a given computational
    method (e.g., exchange-correlation method, basis sets, as defined in section_method).

    The link between the current section calculation and the related
    system and method sections is established by the values stored in system_ref and
    method_ref, respectively.

    The reason why information on the system configuration and computational method is
    stored separately is that several *single configuration calculations* can be performed
    on the same system configuration, viz. several system configurations can be evaluated
    with the same computational method. This storage strategy avoids redundancies.
    '''

    m_def = Section(validate=False, extends_base_section=True)

    x_w2dynamics_dc_latt = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_nd', 2, 'x_w2dynamics_nd', 2],
        description='''
        Double counting correction on the lattice.
        ''')

    x_w2dynamics_gdensnew = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_nd', 2, 'x_w2dynamics_nd', 2],
        description='''
        Densities with new self-energy after adjustment of mu.
        ''')

    x_w2dynamics_gdensold = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_nd', 2, 'x_w2dynamics_nd', 2],
        description='''
        Densities with old self-energy before adjustment of mu.
        ''')

    x_w2dynamics_glocnew_lattice = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_nd', 2, 'Niw'],
        description='''
        Local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands.
        ''')

    x_w2dynamics_glocold_lattice = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_nd', 2, 'Niw'],
        description='''
        Local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands.
        ''')

    x_w2dynamics_mu = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        units='electron_volt',
        description='''
        Chemical potential (lattice model).
        ''')

    x_w2dynamics_ineq = SubSection(sub_section=x_w2dynamics_quantities.m_def, repeats=True)


class Run(simulation.run.Run):

    m_def = Section(validate=False, extends_base_section=True)

    x_w2dynamics_axes = SubSection(sub_section=x_w2dynamics_axes.m_def, repeats=False)

    # TODO add config, environment variables
