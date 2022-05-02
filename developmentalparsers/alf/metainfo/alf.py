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
import numpy as np            # pylint: disable=unused-import

from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.datamodel.metainfo import simulation


m_package = Package()


class x_alf_lattice(MSection):

    m_def = Section(validate=False)

    x_alf_n = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of lattice points.
        ''')

    x_alf_cell = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        ''')


class System(simulation.system.System):

    m_def = Section(validate=False, extends_base_section=True)

    x_alf_lattice = SubSection(sub_section=x_alf_lattice)


class x_alf_observable_data(MSection):

    m_def = Section(validate=False)

    x_alf_n = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_alf_data = Quantity(
        type=np.dtype(np.float64),
        shape=['x_alf_n'],
        description='''
        ''')


class x_alf_observable(MSection):

    m_def = Section(validate=False)

    x_alf_n = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_alf_sign = Quantity(
        type=np.dtype(np.float64),
        shape=['x_alf_n'],
        description='''
        ''')

    x_alf_shape = Quantity(
        type=np.dtype(np.int32),
        shape=['0..*'],
        description='''
        ''')

    x_alf_obser = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')
    # x_alf_obser = SubSection(sub_section=x_alf_observable_data, repeats=True)


class Calculation(simulation.calculation.Calculation):

    m_def = Section(validate=False, extends_base_section=True)

    x_alf_den_eq = SubSection(sub_section=x_alf_observable)

    x_alf_den_tau = SubSection(sub_section=x_alf_observable)

    x_alf_ener_scal = SubSection(sub_section=x_alf_observable)

    x_alf_green_eq = SubSection(sub_section=x_alf_observable)

    x_alf_green_tau = SubSection(sub_section=x_alf_observable)

    x_alf_kin_scal = SubSection(sub_section=x_alf_observable)

    x_alf_part_scal = SubSection(sub_section=x_alf_observable)

    x_alf_pot_scal = SubSection(sub_section=x_alf_observable)

    x_alf_spint_eq = SubSection(sub_section=x_alf_observable)

    x_alf_spint_tau = SubSection(sub_section=x_alf_observable)

    x_alf_spinxy_eq = SubSection(sub_section=x_alf_observable)

    x_alf_spinxy_tau = SubSection(sub_section=x_alf_observable)

    x_alf_spinz_eq = SubSection(sub_section=x_alf_observable)

    x_alf_spinq_tau = SubSection(sub_section=x_alf_observable)
