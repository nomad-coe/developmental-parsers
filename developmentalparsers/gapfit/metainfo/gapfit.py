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
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.datamodel.metainfo import simulation


m_package = Package()


class x_gapfit_e0(MSection):

    m_def = Section(validate=False)

    x_gapfit_Z = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')


class x_gapfit_gap_data(MSection):

    m_def = Section(validate=False)

    x_gapfit_do_core = Quantity(
        type=bool,
        shape=[],
        description='''
        ''')

    x_gapfit_e0 = SubSection(sub_section=x_gapfit_e0, repeats=True)


class x_gapfit_sparsex(MSection):

    m_def = Section(validate=False)

    x_gapfit_i = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_alpha = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')

    x_gapfit_sparseCutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')


class x_gapfit_gp_coordinates(MSection):

    m_def = Section(validate=False)

    x_gapfit_label = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_dimensions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_signal_variance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')

    x_gapfit_signal_mean = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')

    x_gapfit_sparsified = Quantity(
        type=bool,
        shape=[],
        description='''
    ''')

    x_gapfit_n_permutations = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_covariance_type = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_zeta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')

    x_gapfit_n_sparseX = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_sparseX_filename = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_sparseX_md5sum = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_theta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ''')

    x_gapfit_descriptor = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_permutation = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_sparsex = SubSection(sub_section=x_gapfit_sparsex, repeats=True)


class x_gapfit_gpsparse(MSection):

    m_def = Section(validate=False)

    x_gapfit_label = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_n_coordinate = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        ''')

    x_gapfit_fitted = Quantity(
        type=bool,
        shape=[],
        description='''
        ''')

    x_gapfit_gp_coordinates = SubSection(sub_section=x_gapfit_gp_coordinates, repeats=True)


class x_gapfit_gap_params(MSection):

    m_def = Section(validate=False)

    x_gapfit_label = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_gap_version = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_gap_data = SubSection(sub_section=x_gapfit_gap_data, repeats=True)

    x_gapfit_gpsparse = SubSection(sub_section=x_gapfit_gpsparse, repeats=True)


class Run(simulation.run.Run):

    m_def = Section(validate=False, extends_base_section=True)

    x_gapfit_potential_label = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_potential_init_args = Quantity(
        type=str,
        shape=[],
        description='''
        ''')

    x_gapfit_gap_params = SubSection(sub_section=x_gapfit_gap_params)
