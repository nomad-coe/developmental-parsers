#!/bin/bash

module load slurm_setup

source /dss/dsshome1/0E/di36fed2/.w2dynamics_profile
slurm_export_end_time

set -e
set -x

export BETA=60
export NIW=1200
export NFTAU=1200
export NCORR=290
export NWARMUPS=200000
export NWARMUPS2PLUS=
export REUSEMCCONFIG=
export NLOOKUP_NFFT=1000000
export NMEAS_LONGEST=200001

export FNPREFIX="SrVO3_beta${BETA}"
export OLDFILE=$(ls -t SrVO3_beta*hdf* | head -n 1)
export READOLD=-1
export TOTDENS=1.0

mpiexec -n ${SLURM_NTASKS} /dss/dsshome1/0E/di36fed2/clusterversion/DMFT.py General.SLURMOnlyContinueIfTime=50 General.DMFTsteps=3 General.beta="${BETA}" General.FileNamePrefix="${FNPREFIX}" General.mixing_strategy=modbroyden General.mixing_diis_history=100 General.mixing=0.95 General.dump_mixer=mixer.bin General.read_mixer=mixer.bin General.fileold=${OLDFILE} General.readold=${READOLD} General.mu_search=mixer General.totdens=${TOTDENS} General.EPSN=0.000001 QMC.Nmeas=${NMEAS_LONGEST} QMC.Nwarmups=${NWARMUPS} QMC.Niw=${NIW} QMC.Nftau=${NFTAU} QMC.NLookup_nfft=${NLOOKUP_NFFT} QMC.NCorr=${NCORR} Parameters.in
