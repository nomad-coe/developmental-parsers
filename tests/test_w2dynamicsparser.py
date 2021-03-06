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
from developmentalparsers.w2dynamics import W2DynamicsParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return W2DynamicsParser()


def test_dmft(parser):
    archive = EntryArchive()
    parser.parse('tests/data/w2dynamics/srvo3/SrVO3_beta60-2021-12-03-Fri-13-38-46.hdf5', archive, None)

    assert len(archive.run[0].calculation) == 6
    # TODO assert value of quantities
