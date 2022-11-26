#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest

from nomad.datamodel import EntryArchive
from developmentalparsers.gapfit import GAPFitParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return GAPFitParser()


def test_gapfit(parser):
    archive = EntryArchive()
    parser.parse('tests/data/gapfit/gap.xml', archive, None)

    assert archive.run[0].x_gapfit_potential_init_args == 'IP GAP label=GAP_2022_11_2_60_12_30_17_127'
    sec_params = archive.run[0].x_gapfit_gap_params
    assert sec_params.x_gapfit_label == 'GAP_2022_11_2_60_12_30_17_127'
    assert sec_params.x_gapfit_gap_version == '1663247896'
    sec_data =  sec_params.x_gapfit_gap_data[0]
    assert not sec_data.x_gapfit_do_core
    assert len(sec_data.x_gapfit_e0) == 118
    assert sec_data.x_gapfit_e0[13].x_gapfit_value == approx(-1481.2905682064156)
    assert sec_data.x_gapfit_e0[17].x_gapfit_Z == 18
    sec_gpsparse = sec_params.x_gapfit_gpsparse[0]
    assert sec_gpsparse.x_gapfit_n_coordinate == 2
    assert sec_gpsparse.x_gapfit_fitted
    assert len(sec_gpsparse.x_gapfit_gp_coordinates) == 2
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_label == 'GAP_2022_11_2_60_12_30_17_1271'
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_dimensions == 681
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_signal_variance == 10
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_signal_mean == 0
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_sparsified
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_n_permutations == 1
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_covariance_type == 1
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_n_sparseX == 1000
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_sparseX_filename == 'gap.xml.sparseX.GAP_2022_11_2_60_12_30_17_1271'
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_zeta == 2
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_sparseX_md5sum == 'fc0ae82220eb499fe58919abdd346328'
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_descriptor.startswith('soap add_species=F')
    assert sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_theta == 0.75
    assert sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_permutation == 1
    assert len(sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_sparsex) == 100
    assert len(sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_sparsex) == 1000
    sec_gpsparse.x_gapfit_gp_coordinates[0].x_gapfit_sparsex[3].x_gapfit_alpha == approx(-18028445.496617578)
    sec_gpsparse.x_gapfit_gp_coordinates[1].x_gapfit_sparsex[500].x_gapfit_sparseCutoff == 1.0

    sec_system = archive.run[0].system
    assert len(sec_system) == 2203
    assert len(sec_system[1].atoms.positions) == 2
    assert sec_system[10].atoms.positions[4][1].to('angstrom').magnitude == approx(1.34725000)
    assert sec_system[17].atoms.lattice_vectors[1][1].to('angstrom').magnitude == approx(-0.00583853)
    assert sec_system[25].atoms.periodic[0]

    sec_calc = archive.run[0].calculation
    assert len(sec_calc) == 2203
    assert sec_calc[2194].energy.total.value.to('eV').magnitude == approx(-15.216196150000002)