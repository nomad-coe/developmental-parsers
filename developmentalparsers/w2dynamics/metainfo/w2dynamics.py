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


class x_w2dynamics_axes(MSection):

    m_def = Section(validate=False)

    x_w2dynamics_iw = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwb_g4 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwb_p2 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwb_p3 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwf_g2 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwf_g4 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwf_p3 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwb_g4 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_iwb_p3 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_pos_iw = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_tau = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_tau_g4 = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_taubin = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_tauf = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_tausus = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')

    x_w2dynamics_w_dos = Quantity(
        type=np.dtype(np.float64),
        description='''
        ''')


class x_w2dynamics_quantities(MSection):

    m_def = Section(validate=False)

    x_w2dynamics_accept_flavourchange = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for flavourchange moves
        ''')

    x_w2dynamics_accept_glob = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for global symmetry moves
        ''')

    x_w2dynamics_accept_ins = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for insertions into local trace
        ''')

    x_w2dynamics_accept_ins4 = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for insertions of 4 operators into local trace
        ''')

    x_w2dynamics_accept_pair_tau = Quantity(
        type=np.dtype(np.float64),
        description='''
        Total number of accepted pair insertions and removals in tau bins
        ''')

    x_w2dynamics_accept_rem = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for removals of local trace operators
        ''')

    x_w2dynamics_accept_rem4 = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for removal of 4 operators in the local trace
        ''')

    x_w2dynamics_accept_shift = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for tau-shift moves
        ''')

    x_w2dynamics_accept_worm_ins = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for worm insertions into local trace
        ''')

    x_w2dynamics_accept_worm_rem = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for worm removals from local trace
        ''')

    x_w2dynamics_accept_worm_rep = Quantity(
        type=np.dtype(np.float64),
        description='''
        Acceptance rate for worm replacement moves
        ''')

    x_w2dynamics_at2neq = Quantity(
        type=np.dtype(np.float64),
        description='''
        corresponding inequivalent index to atom index
        ''')

    x_w2dynamics_contrib_sst = Quantity(
        type=np.dtype(np.float64),
        description='''
        Trace contribution per outer superstate
        ''')

    x_w2dynamics_contrib_state = Quantity(
        type=np.dtype(np.float64),
        description='''
        Trace contribution per outer state
        ''')

    x_w2dynamics_dc = Quantity(
        type=np.dtype(np.float64),
        description='''
        double counting correction
        ''')

    x_w2dynamics_dc_latt = Quantity(
        type=np.dtype(np.float64),
        description='''
        double counting correction on the lattice
        ''')

    x_w2dynamics_densitymatrix = Quantity(
        type=np.dtype(np.float64),
        description='''
        Density matrix in occupation number basis measured at beta/2
        ''')

    x_w2dynamics_double_occ = Quantity(
        type=np.dtype(np.float64),
        description='''
        Double occupancies <n_i n_j>
        ''')

    x_w2dynamics_energies_eigenstates = Quantity(
        type=np.dtype(np.float64),
        description='''
        Energies of eigenstates by state index
        ''')

    x_w2dynamics_expresdensitymatrix = Quantity(
        type=np.dtype(np.float64),
        description='''
        density matrix resolved by expansion order
        ''')

    x_w2dynamics_fiw = Quantity(
        type=np.dtype(np.float64),
        description='''
        hybridisation function in Matsubara frequencies
        ''')

    x_w2dynamics_fmom = Quantity(
        type=np.dtype(np.float64),
        description='''
        moments of the Hybridisation function in iw
        ''')

    x_w2dynamics_ftau = Quantity(
        type=np.dtype(np.float64),
        description='''
        hybridisation function in imaginary time
        ''')

    x_w2dynamics_ftau_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        hybridisation function in imaginary time
        ''')

    x_w2dynamics_g0iw = Quantity(
        type=np.dtype(np.float64),
        description='''
        Weiss field
        ''')

    x_w2dynamics_g0iw_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        Weiss field including spin/orbital off-diagonal terms
        ''')

    x_w2dynamics_g2iw = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in Matsubara basis
        ''')

    x_w2dynamics_g4iw = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in Matsubara basis (particle-hole channel)
        ''')

    x_w2dynamics_g4iw_pp = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in Matsubara basis (particle-particle channel)
        ''')

    x_w2dynamics_g4iw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in particle-hole Matsubara frequencies from worm sampling
        There are two conventions:
        0: (v+w) tau_1 - v tau_2 + v' tau_3 - (v'+w) tau_4
        1: v tau_1 - (v-w) tau_2 + (v'-w) tau_3 - v' tau_4
        ''')

    x_w2dynamics_g4iwpp_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in particle-particle Matsubara frequencies from worm sampling
        Convention: v tau_1 - (w-v') tau_2 + (w-v) tau_3 - v' tau_4
        ''')

    x_w2dynamics_g4leg = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in Legendre/Matsubara basis
        ''')

    x_w2dynamics_g4tau = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle Green's function in tau12, tau34, tau14
        ''')

    x_w2dynamics_gdensnew = Quantity(
        type=np.dtype(np.float64),
        description='''
        densities with new self-energy after adjustment of mu
        ''')

    x_w2dynamics_gdensold = Quantity(
        type=np.dtype(np.float64),
        description='''
        densities with old self-energy before adjustment of mu
        ''')

    x_w2dynamics_giw = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function used as input for Dyson equation
        ''')

    x_w2dynamics_giw_cov = Quantity(
        type=np.dtype(np.float64),
        description='''
        covariance of diagonal impurity Green's function
        ''')

    x_w2dynamics_giw_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function used as input for Dyson equation
        ''')

    x_w2dynamics_giw_meas = Quantity(
        type=np.dtype(np.float64),
        description='''
        Impurity Green's function in Matsubara
        ''')

    x_w2dynamics_giw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Impurity Green's function in Matsubara from worm sampling
        ''')

    x_w2dynamics_gleg = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function in Legendre expansion
        ''')

    x_w2dynamics_gleg_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        full impurity Green's function in Legendre expansion
        ''')

    x_w2dynamics_glocnew = Quantity(
        type=np.dtype(np.float64),
        description='''
        local Green's function in Matsubara (new self-energy)
        ''')

    x_w2dynamics_glocnew_lattice = Quantity(
        type=np.dtype(np.float64),
        description='''
        local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands
        ''')

    x_w2dynamics_glocold = Quantity(
        type=np.dtype(np.float64),
        description='''
        local Green's function in Matsubara (old self-energy)
        ''')

    x_w2dynamics_glocold_lattice = Quantity(
        type=np.dtype(np.float64),
        description='''
        local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands
        ''')

    x_w2dynamics_gsigmaiw = Quantity(
        type=np.dtype(np.float64),
        description='''
        Improved estimators in Matsubara basis
        ''')

    x_w2dynamics_gsigmaiw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Improved estimators in Matsubara basis
        ''')

    x_w2dynamics_gtau = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function on the imaginary time axis
        ''')

    x_w2dynamics_gtau_blocks = Quantity(
        type=np.dtype(np.float64),
        description='''
        blocking analysis for the impurity Green's function
        ''')

    x_w2dynamics_gtau_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function on the imaginary time axis, full offdiagonal
        ''')

    x_w2dynamics_gtau_mean_step = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function averaged over tau, resolved in time (= Monte carlo steps)
        ''')

    x_w2dynamics_gtau_mid_step = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity Green's function averaged over tau between 0.4 * beta and 0.6 * beta, resolved in time (= Monte carlo steps)
        ''')

    x_w2dynamics_gtau_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Impurity Green's function in imaginary time from worm sampling
        ''')

    x_w2dynamics_h_mean = Quantity(
        type=np.dtype(np.float64),
        description='''
        Mean over k-points of Hamiltonian H(k)
        ''')

    x_w2dynamics_h_mom2 = Quantity(
        type=np.dtype(np.float64),
        description='''
        Second moment of the density of states
        ''')

    x_w2dynamics_h4iw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two-particle improved estimator in Matsubara basis from worm sampling
        ''')

    x_w2dynamics_hist = Quantity(
        type=np.dtype(np.float64),
        description='''
        Histogram for expansion orders
        ''')

    x_w2dynamics_hist_seg = Quantity(
        type=np.dtype(np.float64),
        description='''
        number of empty flavours in trace
        ''')

    x_w2dynamics_hist_sst = Quantity(
        type=np.dtype(np.float64),
        description='''
        Histogram for outer superstates
        ''')

    x_w2dynamics_hist_state = Quantity(
        type=np.dtype(np.float64),
        description='''
        Histogram for outer states
        ''')

    x_w2dynamics_hk = Quantity(
        type=np.dtype(np.float64),
        description='''
        Full Hamiltonian H(k)
        ''')

    x_w2dynamics_hk_mean = Quantity(
        type=np.dtype(np.float64),
        description='''
        Mean of diagonals of Hamiltonian H(k)
        ''')

    x_w2dynamics_jbar = Quantity(
        type=np.dtype(np.float64),
        description='''
        shell-averaged J values
        ''')

    x_w2dynamics_lda_dens = Quantity(
        type=np.dtype(np.float64),
        description='''
        LDA orbital-resolved densities
        ''')

    x_w2dynamics_lda_dos = Quantity(
        type=np.dtype(np.float64),
        description='''
        LDA density of states
        ''')

    x_w2dynamics_lda_mu = Quantity(
        type=np.dtype(np.float64),
        description='''
        chemical potential (non-interacting)
        ''')

    x_w2dynamics_leadsiw_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        Weiss field including spin/orbital off-diagonal terms
        ''')

    x_w2dynamics_lhist = Quantity(
        type=np.dtype(np.float64),
        description='''
        Histogram for Lanczos steps in the time evolution
        ''')

    x_w2dynamics_mu = Quantity(
        type=np.dtype(np.float64),
        description='''
        chemical potential (lattice model)
        ''')

    x_w2dynamics_muimp = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity chemical potential
        ''')

    x_w2dynamics_muimp_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        impurity chemical potential
        ''')

    x_w2dynamics_neq2at = Quantity(
        type=np.dtype(np.float64),
        description='''
        atom index first occurrence of inequivalent atom
        ''')

    x_w2dynamics_nneq = Quantity(
        type=np.dtype(np.float64),
        description='''
        number of non-equivalent atoms in unit cell
        ''')

    x_w2dynamics_nqqdag_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Part of Symmetric Improved 2P estimator in Matsubara frequencies
        ''')

    x_w2dynamics_ntau_n0 = Quantity(
        type=np.dtype(np.float64),
        description='''
        density-density correlation function in tau bins
        ''')

    x_w2dynamics_occ = Quantity(
        type=np.dtype(np.float64),
        description='''
        Occupancies
        ''')

    x_w2dynamics_occbasis_mapping = Quantity(
        type=np.dtype(np.float64),
        description='''
        Occupation number eigenvalues by state index as used in densitymatrix
        ''')

    x_w2dynamics_p2iw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two legged two-particle Green's function ph-convention (1 bosonic)
        ''')

    x_w2dynamics_p2iwpp_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two legged two-particle Green's function pp-convention (1 bosonic)
        ''')

    x_w2dynamics_p2tau_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two legged two-particle Green's function in imaginary time
        ''')

    x_w2dynamics_p2taupp_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Two legged two-particle Green's function pp in imaginary time
        ''')

    x_w2dynamics_p3iw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Three legged two-particle Green's function ph-convention (1 fermionic, 1 bosonic)
        ''')

    x_w2dynamics_p3iwpp_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Three legged two-particle Green's function pp-convention (1 fermionic, 1 bosonic)
        ''')

    x_w2dynamics_qqdd_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Part of Symmetric Improved 2P estimator in Matsubara frequencies
        ''')

    x_w2dynamics_qqiw_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 1P estimator in Matsubara basis
        ''')

    x_w2dynamics_qqqq_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 2P estimators in Matsubara basis
        ''')

    x_w2dynamics_qqtau_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 1P estimator in imaginary time
        ''')

    x_w2dynamics_quddag_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Improved estimators in Matsubara basis
        ''')

    x_w2dynamics_rhist = Quantity(
        type=np.dtype(np.float64),
        description='''
        Histogram for additional Lanczos steps needed for reversible trace
        ''')

    x_w2dynamics_rho1 = Quantity(
        type=np.dtype(np.float64),
        description='''
        < c^+ c >
        ''')

    x_w2dynamics_rho2 = Quantity(
        type=np.dtype(np.float64),
        description='''
        < c^+ c^+ c c >
        ''')

    x_w2dynamics_sigma_hartree = Quantity(
        type=np.dtype(np.float64),
        description='''
        Hartree self-energy
        ''')

    x_w2dynamics_sign = Quantity(
        type=np.dtype(np.float64),
        description='''
        Mean sign
        ''')

    x_w2dynamics_sign_sst = Quantity(
        type=np.dtype(np.float64),
        description='''
        Summed total sign per outer superstate
        ''')

    x_w2dynamics_sign_state = Quantity(
        type=np.dtype(np.float64),
        description='''
        Summed total sign per outer state
        ''')

    x_w2dynamics_sign_step = Quantity(
        type=np.dtype(np.float64),
        description='''
        Sign of a configuration's total weight, resolved in time (= Monte carlo steps)
        ''')

    x_w2dynamics_single_occ = Quantity(
        type=np.dtype(np.float64),
        description='''
        Single occupancies <n_i>
        ''')

    x_w2dynamics_siw = Quantity(
        type=np.dtype(np.float64),
        description='''
        Band-spin-diagonal self-energy in matsubara expansion
        ''')

    x_w2dynamics_siw_cov = Quantity(
        type=np.dtype(np.float64),
        description='''
        covariance of diagonal self-energy in matsubara expansion
        ''')

    x_w2dynamics_siw_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        Full self-energy in matsubara expansion (with jackknife error)
        ''')

    x_w2dynamics_siw_trial = Quantity(
        type=np.dtype(np.float64),
        description='''
        Full trial self-energy in matsubara expansion
        ''')

    x_w2dynamics_smom = Quantity(
        type=np.dtype(np.float64),
        description='''
        Moments of the self-energy (in Matsubara)
        ''')

    x_w2dynamics_smom_full = Quantity(
        type=np.dtype(np.float64),
        description='''
        Moments of the self-energy including all terms
        ''')

    x_w2dynamics_ssts_states = Quantity(
        type=np.dtype(np.float64),
        description='''
        Assignment of states to superstates by state index
        ''')

    x_w2dynamics_steps_worm_partition = Quantity(
        type=np.dtype(np.float64),
        description='''
        Time spent in worm space and partition space
        ''')

    x_w2dynamics_time = Quantity(
        type=np.dtype(np.float64),
        description='''
        Total elapsed walltime of the run
        ''')

    x_w2dynamics_time_g4iw_add = Quantity(
        type=np.dtype(np.float64),
        description='''
        Time spent on constructing and adding G4(iw,iw,iW)
        ''')

    x_w2dynamics_time_g4iw_ft = Quantity(
        type=np.dtype(np.float64),
        description='''
        Time spent on fourier-transforming G4(iw,iw,iW)
        ''')

    x_w2dynamics_time_giw = Quantity(
        type=np.dtype(np.float64),
        description='''
        Time spent on measuring G(iw)
        ''')

    x_w2dynamics_time_qmc = Quantity(
        type=np.dtype(np.float64),
        description='''
        Mean CPU seconds spent in CT-QMC per process
        ''')

    x_w2dynamics_time_sampling = Quantity(
        type=np.dtype(np.float64),
        description='''
        CPU time used for the QMC simulation excluding measurements
        ''')

    x_w2dynamics_time_simulation = Quantity(
        type=np.dtype(np.float64),
        description='''
        CPU time used for the QMC simulation
        ''')

    x_w2dynamics_time_warmup = Quantity(
        type=np.dtype(np.float64),
        description='''
        CPU time used for the QMC warmup
        ''')

    x_w2dynamics_ubar = Quantity(
        type=np.dtype(np.float64),
        description='''
        shell-averaged U values
        ''')

    x_w2dynamics_ucaca_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 2P estimators in Matsubara basis
        ''')

    x_w2dynamics_ucacatau_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 2P estimators in Matsubara basis
        ''')

    x_w2dynamics_uccaa_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 2P estimators in Matsubara basis
        ''')

    x_w2dynamics_uccaatau_worm = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm Symmetric Improved 2P estimators in Matsubara basis
        ''')

    x_w2dynamics_worm_eta = Quantity(
        type=np.dtype(np.float64),
        description='''
        Worm balancing factor
        ''')


class x_w2dynamics_atom_parameters(MSection):

    m_def = Section(validate=False)

    x_w2dynamics_atom = Quantity(
        type=JSON,
        description='''
        ''')


class x_w2dynamics_input_parameters(MSection):

    m_def = Section(validate=False)

    x_w2dynamics_general = Quantity(
        type=JSON,
        description='''
        ''')

    x_w2dynamics_qmc = Quantity(
        type=JSON,
        description='''
        ''')


class System(simulation.system.System):
    '''
    Contains parameters describing a system of atomic configuration. These inclue the
    compound name, atomic positions, lattice vectors, contraints on the atoms, etc.
    '''

    m_def = Section(validate=False, extends_base_section=True)

    # TODO extend this quantity
    dos = Quantity(
        type = str,
        shape = [],
        description = '''
        DOS choice
        ''')


class DMFT(MSection):
    '''
    Section containing the various parameters that define a DMFT calculation
    '''

    m_def = Section(validate=False)

    x_w2dynamics_input_parameters = SubSection(sub_section=x_w2dynamics_input_parameters.m_def, repeats=False)

    x_w2dynamics_atom_parameters = SubSection(sub_section=x_w2dynamics_atom_parameters.m_def, repeats=True)

    local_hamiltonian = Quantity(
        type = MEnum(["Density", "Kanamori", "Coulomb", "ReadUmatrix", "ReadNormalUmatrix"]),
        shape = [],
        description = '''
        Form of the local Coulomb interacting Hamiltonian:

        | Name                  | Description                               |
        | --------------------- | ----------------------------------------- |
        | `"Density"`           | Density-density interaction               |
        | `"Kanamori"`          | Slater-Kanamori interaction               |
        | `"Coulomb"`           | Slater integrals F0, F2, F4, F6 specified |
        | `"ReadUmatrix"`       | Read umatrix file (spin-dependent)        |
        | `"ReadNormalUmatrix"` | Read umatrix file (spin-independent)      |
        ''')

    number_of_atoms_per_unit_cell = Quantity(
        type = int,
        shape = [],
        description = '''
        Number of atoms per unit cell
        ''')

    number_of_correlated_bands = Quantity(
        type = int,
        shape = [],
        description = '''
        Number of correlated bands for the local interactions.
        ''')

    number_of_correlated_electrons = Quantity(
        type = np.dtype(np.float64),
        shape = [],
        description = '''
        Number of valence electrons per atom in the correlated bands
        ''')

    U = Quantity(
        type = np.dtype(np.float64),
        shape = [],
        unit = 'electron_volt',
        description = '''
        Value of the intra-orbital local interaction for d-orbitals
        ''')

    Up = Quantity(
        type = np.dtype(np.float64),
        shape = [],
        unit = 'electron_volt',
        description = '''
        Value of the inter-orbital local interaction between d-orbitals
        ''')

    JH = Quantity(
        type = np.dtype(np.float64),
        shape = [],
        unit = 'electron_volt',
        description = '''
        Value of the Hund's coupling between d-orbitals
        ''')

    temperature = Quantity(
        type = np.dtype(np.float64),
        shape = [],
        unit = 'kelvin',
        description = '''
        Temperature in Kelvin calculated from the inverse temperature `beta` as:
            T = 1/(kB*beta)
        ''')

    magnetic_state = Quantity(
        type = MEnum(["para", "ferro", "antiferro"]),
        shape = [],
        description = '''
        Magnetic state in which the DMFT calculation is done:

        | Name          | State                   |
        | ------------- | ----------------------- |
        | `"para"`      | paramagnetic state      |
        | `"ferro"`     | ferromagnetic state     |
        | `"antiferro"` | antiferromagnetic state |
        ''')

    number_of_tau = Quantity(
        type = int,
        shape = [],
        description = '''
        Number of tau (imaginary times)
        ''')

    number_of_matsubara_freq = Quantity(
        type = int,
        shape = [],
        description = '''
        Number of Matsubara frequencies (imaginary frequencies)
        ''')

    self_energy_mixing = Quantity(
        type = MEnum(["anisimov", "fll", "amf", "trace", "siginfbar", "sigzerobar"]),
        shape = [],
        description = '''
        LDA+DMFT method for mixing the self-energy:

        | Name          | Description                   |
        | ------------- | ----------------------- |
        | `"anisimov"`      | ???????      |
        | `"fll"`      | ???????      |
        | `"amf"`      | ???????      |
        | `"trace"`      | ???????      |
        | `"siginfbar"`      | ???????      |
        | `"sigzerobar"`      | ???????      |
        ''')

    impurity_solver = Quantity(
        type = MEnum(["CT-INT", "CT-HYB", "CT-AUX", "ED", "NRG", "MPS", "IPT", "NCA", "OCA", "slave_bosons"]),
        shape = [],
        description = '''
        Impurity solver for the single impurity Anderson model:

        | Name          | Description                   |
        | ------------- | ----------------------- |
        | `"CT-INT"`      | ???????      |
        | `"CT-HYB"`      | ???????      |
        | `"CT-AUX"`      | ???????      |
        | `"ED"`      | ???????      |
        | `"NRG"`      | ???????      |
        | `"MPS"`      | ???????      |
        | `"IPT"`      | ???????      |
        | `"NCA"`      | ???????      |
        | `"OCA"`      | ???????      |
        | `"slave_bosons"`      | ???????      |
        ''')


class Method(simulation.method.Method):
    '''
    Section containing the various parameters that define the theory and the
    approximations (convergence, thresholds, etc.) behind the calculation.
    '''

    m_def = Section(validate=False, extends_base_section=True)

    dmft = SubSection(sub_section=DMFT.m_def, repeats=False)


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
        shape=['x_w2dynamics_band1', 'x_w2dynamics_spin1', 'x_w2dynamics_band2', 'x_w2dynamics_spin2'],
        description='''
        double counting correction on the lattice
        ''')

    x_w2dynamics_gdensnew = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_lda_band1', 'x_w2dynamics_spin1', 'x_w2dynamics_lda_band2', 'x_w2dynamics_spin2'],
        description='''
        densities with new self-energy after adjustment of mu
        ''')

    x_w2dynamics_gdensold = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_lda_band1', 'x_w2dynamics_spin1', 'x_w2dynamics_lda_band2', 'x_w2dynamics_spin2'],
        description='''
        densities with old self-energy before adjustment of mu
        ''')

    x_w2dynamics_glocnew_lattice = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_lda_band', 'x_w2dynamics_spin', 'x_w2dynamics_iw'],
        description='''
        local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands
        ''')

    x_w2dynamics_glocold_lattice = Quantity(
        type=np.dtype(np.float64),
        shape=['x_w2dynamics_lda_band', 'x_w2dynamics_spin', 'x_w2dynamics_iw'],
        description='''
        local Green's function in Matsubara (old self-energy), diagonal part for all lda-bands
        ''')

    x_w2dynamics_mu = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        chemical potential (lattice model)
        ''')

    x_w2dynamics_ineq = SubSection(sub_section=x_w2dynamics_quantities.m_def, repeats=True)


class Run(simulation.run.Run):

    m_def = Section(validate=False, extends_base_section=True)

    x_w2dynamics_axes = SubSection(sub_section=x_w2dynamics_axes.m_def, repeats=False)

    x_w2dynamics_quantities = SubSection(sub_section=x_w2dynamics_quantities.m_def, repeats=True)

    # TODO add config, environment variables
