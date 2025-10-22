# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N BIOS_act9test
#PBS -P rp23
# express / normal / copyq (2x24, cascadelake)
#PBS -q normal
# Typical for global or Aust continent at 0.25, 192 GB memory and 48 cpus,
# maybe 12 hours walltime
# Typical for small runs, fewer cpus than pixels
#PBS -l walltime=24:00:00
#PBS -l mem=20GB
#PBS -l ncpus=1
# #PBS -l jobfs=1GB
#PBS -l storage=gdata/rp23+scratch/rp23
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -S /bin/bash
#PBS -M lachlan.whyborn@anu.edu.au
#PBS -m ae

source run_cable-pop_lib.sh

# MPI run or single processor run
# nproc should fit with job tasks
dompi=0   # 0: normal run: ./cable
          # 1: MPI run: mpiexec -n ${nproc} ./cable_mpi
nproc=1   # Number of cores for MPI runs
          # must be same as above: SBATCH -n nproc or PBS -l ncpus=nproc

# --------------------------------------------------------------------
#
# Full Cable run with biomass spinup, POP, land-use change, etc.
#
# This script uses CRU-JRA forcing.
#
# Global meteo and land-use change data can be used with a mask giving land points.
# In step 0, the land mask can be extracted for a single point, an area
# or a number of random points can be chosen.
# Alternatively, single site met, LUH2, forcing and mask can be extracted from the global data sets.
#
# The run sequence is as follows:
#   1. Create a climate restart file using Cable's default vegetation distribution.
#   2. First phase of spinup with static land use, fixed atmospheric CO2 and
#      N deposition from 1700, and 30 years of repeated meteorology.
#      Zero initial biomass stocks, POP and climate.
#   3. Bring biomass stocks into equilibrium, restricting labile P and mineral N pools.
#      Repeat a and b several times.
#      a) Start from restart files.
#      b) Use dump files for the biophysics. Quasi-equilibrium of soil and litter pools
#         using an analytic steady-state solution.
#   4. Same as 3 but without any restriction on labile P and mineral N pools.
#      Repeat a and b several times.
#   5. Second phase of spinup with dynamic land use, atmospheric CO2 and N deposition.
#      a) Dynamic land use from 1580 to 1699, using still fixed atmospheric CO2 and
#         N deposition from 1700, and 30 years of repeated meteorology.
#      b) Run from 1700 to 1899 with dynmic land use, varying atmospheric CO2 and N deposition,
#         but still with 30 years of repeated meteorology.
#   6. Final historical run, everything dynamic from 1900 to 2017.
#   7. Future run, everything dynamic (not all met types).
#
# Written,  Matthias Cuntz, Aug 2019, following the run scripts and namelists provided by Vanessa Haverd
# Modified, Jurgen Knauer, 2020      - gm_explicit, coordination, acclimation
#                                    - bios, plume, future runs
#           Matthias Cuntz, Mar 2021 - functions into run_cable-pop_lib.sh
#
# --------------------------------------------------------------------

#ASKJK - changes in comparison to gm_acclim_coord
# 1. Removed switch Test_parallel=1 -> used [[ dompi -eq 1 ]]
# 2. ScriptsPath="${workpath}/scripts"
#    but
#    namelistpath="$(dirname ${workpath})/namelists"
#    with
#      workpath="/home/599/jk8585/CABLE_run/gm_acclim_coord/global_runs"
#      cablehome="/home/599/jk8585/CABLE_code"
#    -> changed to similar of ScriptsPath="$(dirname ${workpath})/scripts"
# 3. Need plume.nml, bios.nml, gm_LUT_*.nc
# 4. output%grid = "mask" (in gm_acclim_coord) or "land" (default before) [PB:Always has to be "mask"]
# 5. What should be for Run in bios.nml after 1. Climate restart?
# 6. Why is YearEnd different for plume (1849) compared to cru (1699) in 5a. First dynamic land use?
# 7. Do we need chunking in 5a, 5b, 6, and 7? [PB: unnecessary at this point, but possibly for 0.05degs]
# 8. Do we need cropping output to latlon region at the end: is this not in step 0 with ${doextractsite} -eq 1? [PB: Not needed for # BIOS
#
#ASKJK - changes in comparison to gm_acclim_coord

# --------------------------------------------------------------------
# Sequence switches
#

# Step 1
doclimate=1     # 1/0: Do/Do not create climate restart file
# Step 2
dofromzero=1    # 1/0: Do/Do not first spinup phase from zero biomass stocks
# Step 3
doequi1=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
nequi1=4       #      number of times to repeat steps in doequi1
# Step 4
doequi2=1       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
nequi2=4        #      number of times to repeat steps in doequi2
# Step 5a
# [TODO] change this to run LUC spin-up when LUC is active
doiniluc=1      # 1/0: Do/Do not spinup with dynamic land use (5a)
# Step 5b
doinidyn=1      # 1/0: Do/Do not full dynamic spinup from 1700 to 1899 (5b)
# Step 6
dofinal=1       # 1/0: Do/Do not final run from 1900 to 2017
# Step 7
dofuture=0      # 1/0: Do/Do not future runs (plume only)

# --------------------------------------------------------------------
# Other switches
#

# MetType
mettype='bios'     # 'cru', 'plume', 'bios'

# Cable
explicit_gm=0       # 1/0: explicit (finite) or implicit mesophyll conductance
use_LUTgm=1         # 1/0: Do/Do not use lookup table for parameter conversion accounting for gm (only used if explicit_gm=1)
Rubisco_params="Bernacchi_2002"   # "Bernacchi_2002" or "Walker_2013"
coordinate_photosyn=1 # 1/0: Do/Do not coordinate photosynthesis
coord=F               # T/F: version of photosyn. optimisation (optimised(F) or forced (T))
acclimate_photosyn=1  # 1/0: Do/Do not acclimate photosynthesis
call_pop=1          # 1/0: Do/Do not use POP population dynamics model, coupled to CASA
doc13o2=0           # 1/0: Do/Do not calculate 13C
c13o2_simple_disc=0 # 1/0: simple or full 13C leaf discrimination

# --------------------------------------------------------------------
# Setup
#

trap cleanup 1 2 3 6

pid=$$
isdir="${PWD}"
prog=$0
pprog=$(basename ${prog})
pdir=$(dirname ${prog})
tmp=${TMPDIR:-"/tmp"}
system=$(echo ${system} | tr A-Z a-z)
sys=${system#*@}
user=${system%@*}
sys="gadi"
# Special things on specific computer system such as loading modules

pdir=${isdir}
. /etc/bashrc
module purge
# module load intel-compiler/2019.5.281
# module load intel-mpi/2019.5.281
# module load netcdf/4.6.3
# module load intel-compiler/2021.5.0
# module load intel-mpi/2021.5.1
# module load netcdf/4.8.0
# # module load hdf5/1.10.5
module load intel-compiler-llvm/2023.0.0
module load intel-mpi/2021.8.0
module load netcdf/4.9.2
export mpiexecdir=/apps/intel-mpi/2019.5.281/intel64/bin
if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi

# Directories of things
# Relative directories must be relative to the directory of this script,
#   not relative to the directory from which this script is launched (if different)
#   nor relative to the run path.
#

# Run directory: runpath="${sitepath}/run"
experiment=ext_newblazempi_off
sitepath="/scratch/rp23/lw5085/BIOStests/${experiment}" # [TODO] Location to write results
workpath="/scratch/rp23/lw5085/BIOStests/runs" # [TODO] # Location of namelists+parameter files
cablehome="${HOME}/CABLE-POP-BIOS/" # [TODO] Location of code
# Cable executable
if [[ ${dompi} -eq 1 ]] ; then
  exe="${cablehome}/offline/cable-mpi"
else
  exe="${cablehome}/offline/cable"
fi

# CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
aux=""
BlazeDataPath="/g/data/rp23/experiments/2024-04-17_BIOS3-merge/Data_BLAZE"
# Global Mask
SurfaceFile="/g/data/rp23/experiments/2024-04-17_BIOS3-merge/ag9761/data_005/gridinfo_CSIRO_CRU005x005_4tiles_v2.nc"

# Global Met
# For specific tests [TODO]
GlobalLandMaskFile="/g/data/rp23/experiments/2024-04-17_BIOS3-merge/BIOS3_forcing/acttest9/acttest9" # no file extension
GlobalMetPath="/g/data/rp23/experiments/2024-04-17_BIOS3-merge/BIOS3_forcing/acttest9/met/"          # last slash is needed - updated 29/3/2024
ParamPath="/g/data/rp23/experiments/2024-04-17_BIOS3-merge/BIOS3_forcing/acttest9/params/"           # only in bios.nml
GlobalTransitionFilePath="/g/data/rp23/experiments/2024-04-17_BIOS3-merge/LUH2/v3h/0.05deg_aust/EXTRACT"
 
# Other scripts
ScriptsPath="${cablehome}/scripts"

runpath="${sitepath}/run"

# Cable parameters
namelistpath="${workpath}/namelists_bios"
filename_veg="${workpath}/params_bios/def_veg_params.txt"
filename_soil="${workpath}/params_bios/def_soil_params.txt"
casafile_cnpbiome="${workpath}/params_bios/pftlookup.csv"
gm_lut_bernacchi_2002=/g/data/rp23/data/no_provenance/parameter/gm_LUT_351x3601x7_1pt8245_Bernacchi2002.nc
gm_lut_walker_2013=/g/data/rp23/data/no_provenance/parameter/gm_LUT_351x3601x7_1pt8245_Walker2013.nc

ClimateFile="${sitepath}/mask/bios_climate_rst.nc"
MetPath=$(abspath ${GlobalMetPath})
TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
LandMaskFile=$(absfile ${GlobalLandMaskFile})

# --------------------------------------------------------------------
# Start Script
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Helper functions, most functions are in plumber_cable-pop_lib.sh
#

# usage of script
function usage()
{
    printf "${pprog} [-h]\n"
    printf "Runs Cable on a single grid cell with spinup, POP, land-use change, etc.\n"
    printf "Behaviour of the script is controlled by switches at the top of the script (ca. line 101ff).\n"
    printf "\n"
    printf "Options\n"
    printf "    -h    Prints this help screen.\n"
}

# cleanup at end or at trap
function cleanup()
{
    \rm -f ${tmp}/*.${pid}*
    exit 1
}

# returns argument to extract lat and lon with ncks
function nckslatlon()
{
    vars=$(ncvarlist ${1})
    if [[ -z $(isin latitude ${vars}) ]] ; then ilat="lat" ; else ilat="latitude" ; fi
    if [[ -z $(isin longitude ${vars}) ]] ; then ilon="lon" ; else ilon="longitude" ; fi
    if [[ -z $(echo ${2} | cut -f 3 -d ",") || -z $(echo ${2} | cut -f 4 -d ",") ]] ; then
        iilat=$(echo ${2} | cut -f 1 -d ",")
        iilon=$(echo ${2} | cut -f 2 -d ",")
        echo "-d ${ilat},${iilat} -d ${ilon},${iilon}"
    else
        iilat1=$(echo ${2} | cut -f 1 -d ",")
        iilat2=$(echo ${2} | cut -f 2 -d ",")
        iilon1=$(echo ${2} | cut -f 3 -d ",")
        iilon2=$(echo ${2} | cut -f 4 -d ",")
        echo "-d ${ilat},${iilat1},${iilat2} -d ${ilon},${iilon1},${iilon2}"
    fi
}

# --------------------------------------------------------------------------------------------------
# Preparation
#
# Get options
while getopts "h" option ; do
    case ${option} in
        h) usage; exit;;
        *) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2;  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

#
# get directories
mkdir -p ${sitepath}/mask
pdir=$(abspath ${pdir})
cd ${pdir}
adir=$(abspath ${aux})
exe=$(absfile ${exe})
mkdir -p ${runpath}
rdir=$(abspath ${runpath})
ndir=$(abspath ${namelistpath})
sdir=$(abspath ${ScriptsPath})

#
# prepare run directory
cd ${rdir}
mkdir -p logs
mkdir -p outputs
mkdir -p restart
ln -sf ${adir}
# ln -sf ${exe}
cp ${exe} ./
iexe=$(basename ${exe})
cd ${pdir}

#
# set stacksize to unlimited if permitted, otherwise to 15 bit if possible
set +e
ulimit -s unlimited 2> /dev/null || ulimit -s 32768
set -e

# --------------------------------------------------------------------
# Info
#
t1=$(date +%s)
printf "Started at %s\n" "$(date)"

printf "\nSetup\n"
printf "    Serial / Parallel\n"
printf "        dompi=${dompi}\n"
printf "            nproc=${nproc}\n"
printf "\n"
printf "    Sequence\n"
printf "        imeteo=${imeteo}\n"
printf "        doextractsite=${doextractsite}\n"
printf "            experiment=${experiment}\n"
printf "            randompoints=${randompoints}\n"
printf "            latlon=${latlon}\n"
printf "        doclimate=${doclimate}\n"
printf "        dofromzero=${dofromzero}\n"
printf "        doequi1=${doequi1}\n"
printf "            nequi1=${nequi1}\n"
printf "        doequi2=${doequi2}\n"
printf "            nequi2=${nequi2}\n"
printf "        doiniLUC=${doiniLUC}\n"
printf "        doinidyn=${doinidyn}\n"
printf "        dofinal=${dofinal}\n"
printf "        dofuture=${dofuture}\n"
printf "\n"
printf "    Options\n"
printf "        mettype=${mettype}\n"
printf "        metmodel=${metmodel}\n"
printf "        RCP=${RCP}\n"
printf "        explicit_gm=${explicit_gm}\n"
printf "        use_LUTgm=${use_LUTgm}\n"
printf "        Rubisco_params=${Rubisco_params}\n"
printf "        coordinate_photosyn=${coordinate_photosyn}\n"
printf "        coord=${coord}\n"
printf "        acclimate_photosyn=${acclimate_photosyn}\n"
printf "        call_pop=${call_pop}\n"
printf "        doc13o2=${doc13o2}\n"
printf "        c13o2_simple_disc=${c13o2_simple_disc}\n"
printf "\n"
printf "    Directories\n"
printf "        sitepath=${sitepath}\n"
printf "        cablehome=${cablehome}\n"
printf "        exe=${exe}\n"
printf "        aux=${aux}\n"
printf "        GlobalLandMaskFile=${GlobalLandMaskFile}\n"
printf "        SurfaceFile=${SurfaceFile}\n"
printf "        GlobalMetPath=${GlobalMetPath}\n"
printf "        GlobalTransitionFilePath=${GlobalTransitionFilePath}\n"
printf "        runpath=${runpath}\n"
printf "        namelistpath=${namelistpath}\n"
printf "        filename_veg=${filename_veg}\n"
printf "        filename_soil=${filename_soil}\n"
printf "        casafile_cnpbiome=${casafile_cnpbiome}\n"
printf "        LandMaskFile=${LandMaskFile}\n"
printf "        MetPath=${MetPath}\n"
printf "        ClimateFile=${ClimateFile}\n"
printf "        TransitionFilePath=${TransitionFilePath}\n"
printf "        gm_lut_bernacchi_2002=${gm_lut_bernacchi_2002}\n"
printf "        gm_lut_walker_2013=${gm_lut_walker_2013}\n"
#printf "        filename_d13c_atm=${filename_d13c_atm}\n"
printf "\n"

# --------------------------------------------------------------------

# absolute pathes of other parameter files
ClimateFile=$(absfile ${ClimateFile})
filename_veg=$(absfile ${filename_veg})
filename_soil=$(absfile ${filename_soil})
casafile_cnpbiome=$(absfile ${casafile_cnpbiome})
gm_lut_bernacchi_2002=$(absfile ${gm_lut_bernacchi_2002})
gm_lut_walker_2013=$(absfile ${gm_lut_walker_2013})
# This is not necessary while c13o2 is .false.
#filename_d13c_atm=$(absfile ${filename_d13c_atm})
if [[ "${Rubisco_params}" == "Bernacchi_2002" ]] ; then
    filename_gm_lut=${gm_lut_bernacchi_2002}
elif [[ "${Rubisco_params}" == "Walker_2013" ]] ; then
    filename_gm_lut=${gm_lut_walker_2013}
else
    filename_gm_lut=""
fi


# delete all restart files if required
if [[ ${purge_restart} -eq 1 ]] ; then
    rm -f ${rdir}/restart/*
fi

cat > ${tmp}/sedtmp.${pid} << EOF
  met_path         = "${MetPath}/"
  param_path       = "${ParamPath}"
  landmaskflt_file = "${GlobalLandMaskFile}.flt"
  landmaskhdr_file = "${GlobalLandMaskFile}.hdr"
  rain_file        = "1900010120231231_rain_recal_b2311.bin"
  swdown_file      = "1900010120231231_rad_b2311.bin"
  tairmax_file     = "1900010120231231_tmax_noclim_b2311.bin"
  tairmin_file     = "1900010120231231_tmin_noclim_b2311.bin"
  wind_file        = "1900010120231231_windspeed_ms_b2311.bin"
  vp0900_file      = "1900010120231231_vph09_b2311.bin"
  vp1500_file      = "1900010120231231_vph15_b2311.bin"
  co2_file         = "1700_2023_trendy_global_co2_ann.bin"
EOF
applysed ${tmp}/sedtmp.${pid} ${ndir}/bios.nml ${rdir}/bios_${experiment}.nml

# global landuse change namelist
cat > ${tmp}/sedtmp.${pid} << EOF
    TransitionFilePath = "${TransitionFilePath}"
    ClimateFile        = "${ClimateFile}"
    YearStart          = 1700
    YearEnd            = 2017
EOF
applysed ${tmp}/sedtmp.${pid} ${ndir}/LUC.nml ${rdir}/LUC_${experiment}.nml

# global Cable namelist
cat > ${tmp}/sedtmp.${pid} << EOF
    filename%met                       = "${mettype}"
    filename%veg                       = "${filename_veg}"
    filename%soil                      = "${filename_soil}"
    filename%type                      = "${SurfaceFile}"
    filename%out                       = "outputs/${mettype}_out_cable.nc"
    filename%restart_in                = "restart/${mettype}_cable_rst.nc"
    filename%restart_out               = "restart/${mettype}_cable_rst.nc"
    casafile%cnpbiome                  = "${casafile_cnpbiome}"
    casafile%out                       = "outputs/${mettype}_out_casa.nc"
    casafile%cnpipool                  = "restart/${mettype}_casa_rst.nc"
    casafile%cnpepool                  = "restart/${mettype}_casa_rst.nc"
    cable_user%CASA_OUT_FREQ           = "monthly"
    cable_user%POP_restart_in          = "restart/pop_${mettype}_ini.nc"
    cable_user%POP_restart_out         = "restart/pop_${mettype}_ini.nc"
    cable_user%LUC_restart_in          = "restart/${mettype}_LUC_rst.nc"
    cable_user%LUC_restart_out         = "restart/${mettype}_LUC_rst.nc"
    cable_user%LUC_outfile             = "outputs/${mettype}_out_LUC.nc"
    cable_user%climate_restart_in      = "restart/${mettype}_climate_rst.nc"
    cable_user%climate_restart_out     = "restart/${mettype}_climate_rst.nc"
    cable_user%RunIden                 = "${mettype}"
    cable_user%MetType                 = "${mettype}"
    output%averaging                   = "monthly"
    output%grid                        = "land"
    leaps                              = .false.
    cable_user%SOIL_STRUC              = "sli"
    cable_user%Rubisco_parameters      = "${Rubisco_params}"
    cable_user%CALL_POP                = .false.
    cable_user%coordinate_photosyn     = .false.
    cable_user%acclimate_photosyn      = .false.
    cable_user%explicit_gm             = .false.
    cable_user%gm_LUT_file             = "${filename_gm_lut}"
    cable_user%c13o2                   = .false.
    cable_user%c13o2_simple_disc       = .false.
    cable_user%c13o2_delta_atm_file    = "${filename_d13c_atm}"
    cable_user%c13o2_outfile           = "outputs/${mettype}_out_casa_c13o2.nc"
    cable_user%c13o2_restart_in_flux   = "restart/${mettype}_c13o2_flux_rst.nc"
    cable_user%c13o2_restart_out_flux  = "restart/${mettype}_c13o2_flux_rst.nc"
    cable_user%c13o2_restart_in_pools  = "restart/${mettype}_c13o2_pools_rst.nc"
    cable_user%c13o2_restart_out_pools = "restart/${mettype}_c13o2_pools_rst.nc"
    cable_user%c13o2_restart_in_LUC    = "restart/${mettype}_c13o2_LUC_rst.nc"
    cable_user%c13o2_restart_out_LUC   = "restart/${mettype}_c13o2_LUC_rst.nc"
    cable_user%CALL_BLAZE              = .false.
EOF
if [[ ${call_pop} -eq 1 ]] ; then
    sed -i -e "/cable_user%CALL_POP/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${coordinate_photosyn} -eq 1 ]] ; then
    sed -i -e "/cable_user%coordinate_photosyn/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${acclimate_photosyn} -eq 1 ]] ; then
    sed -i -e "/cable_user%acclimate_photosyn/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${explicit_gm} -eq 1 ]] ; then
    sed -i -e "/cable_user%explicit_gm/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
fi
if [[ ${doc13o2} -eq 1 ]] ; then
    sed -i -e "/cable_user%c13o2/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
    if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
        sed -i -e "/cable_user%c13o2_simple_disc/s/=.*/= .true./" ${tmp}/sedtmp.${pid}
    fi
fi
applysed ${tmp}/sedtmp.${pid} ${ndir}/cable.nml ${rdir}/cable_${experiment}.nml


# --------------------------------------------------------------------
# Sequence
#

# --------------------------------------------------------------------
# 1. Create climate restart file
if [[ ${doclimate} -eq 1 ]] ; then
    echo "1. Create climate restart file"
    rid="climate_restart"
    # Met forcing
    cat > ${tmp}/sedtmp.${pid} << EOF
         Run = "spinup"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
    # LUC
    cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #   do not calculate 13C because there is no 13C in the climate restart file
    #MCTEST
    # cable_user%YearEnd = 1889
    # cable_user%CASA_SPIN_ENDYEAR = 1869
    #MCTEST
    cat > ${tmp}/sedtmp.${pid} << EOF
        filename%restart_in            = ""
        cable_user%CLIMATE_fromZero    = .true.
        cable_user%YearStart           = 1860
        cable_user%YearEnd             = 1889
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .true.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .true.
        cable_user%CASA_SPIN_STARTYEAR = 1860
        cable_user%CASA_SPIN_ENDYEAR   = 1869
        cable_user%limit_labile        = .true.
        casafile%cnpipool              = ""
        cable_user%POP_fromZero        = .true.
        cable_user%POP_out             = "ini"
        cable_user%POP_restart_in      = ""
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "static"
        cable_user%c13o2               = .false.
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    imv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc
    cp ${mettype}_climate_rst.nc ${ClimateFile}
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
    cd ..
    cd ${pdir}
fi

# --------------------------------------------------------------------
# 2. First spinup phase from zero biomass
if [[ ${dofromzero} -eq 1 ]] ; then
    echo "2. First spinup from zero biomass"
    rid="zero_biomass"
    # Met forcing
    cat > ${tmp}/sedtmp.${pid} << EOF
         Run = "spinup"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml

    # LUC
    cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #MCTEST
    # cable_user%YearEnd = 1889
    # cable_user%CASA_SPIN_ENDYEAR = 1869
    # remove
    #   cable_user%CASA_OUT_FREQ
    #   output%averaging
    #MCTEST
    cat > ${tmp}/sedtmp.${pid} << EOF
        filename%restart_in               = ""
        cable_user%CLIMATE_fromZero       = .true.
        cable_user%YearStart              = 1860
        cable_user%YearEnd                = 1889
        icycle                            = 2
        spincasa                          = .false.
        cable_user%CASA_OUT_FREQ          = "monthly"
        cable_user%CASA_fromZero          = .true.
        cable_user%CASA_DUMP_READ         = .false.
        cable_user%CASA_DUMP_WRITE        = .true.
        cable_user%CASA_SPIN_STARTYEAR    = 1860
        cable_user%CASA_SPIN_ENDYEAR      = 1869
        cable_user%limit_labile           = .true.
        casafile%cnpipool                 = ""
        cable_user%POP_fromZero           = .true.
        cable_user%POP_out                = "ini"
        cable_user%POP_restart_in         = ""
        cable_user%POPLUC                 = .true.
        cable_user%POPLUC_RunType         = "static"
        cable_user%c13o2_restart_in_flux  = ""
        cable_user%c13o2_restart_in_pools = ""
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    imv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 3. Biomass into quasi-equilibrium with restricted N and P pools
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with restricted N and P pools"
    for ((iequi1=1; iequi1<=${nequi1}; iequi1++)) ; do
        # 3a. 30 year run starting from restart files
        echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
        #rid="spinup_limit_labile"
        rid="spinup_limit_labile${iequi1}"
        # Met forcing
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml

        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd = 1859
        # cable_user%CASA_SPIN_ENDYEAR = 1869
        #MCTEST
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1860
            cable_user%YearEnd             = 1889
            icycle                         = 2
            spincasa                       = .false.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .false.
            cable_user%CASA_DUMP_WRITE     = .true.
            cable_user%CASA_SPIN_STARTYEAR = 1860
            cable_user%CASA_SPIN_ENDYEAR   = 1869
            cable_user%limit_labile        = .true.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
        cd ..
        cd ${pdir}
        #
        # 3b. analytic quasi-equilibrium of biomass pools
        echo "   3b. Analytic solution of biomass pools"
        #rid="spinup_analytic_limit_labile"
        rid="spin_casa_limit_labile${iequi1}"
        # Met forcing
        if [[ "${mettype}" == "cru" ]] ; then
            cp ${rdir}/cru_${experiment}.nml ${rdir}/cru.nml
        elif [[ "${mettype}" == "plume" ]] ; then
            cp ${rdir}/plume_${experiment}.nml ${rdir}/plume.nml
        elif [[ "${mettype}" == "bios" ]] ; then
            cat > ${tmp}/sedtmp.${pid} << EOF
    	          Run = "spinup"
EOF
            applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        fi
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd = 1859
        # cable_user%CASA_SPIN_ENDYEAR = 1859
        #MCTEST
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1860
            cable_user%YearEnd             = 1889
            icycle                         = 12
            spincasa                       = .true.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .true.
            cable_user%CASA_DUMP_WRITE     = .false.
            cable_user%CASA_SPIN_STARTYEAR = 1860
            cable_user%CASA_SPIN_ENDYEAR   = 1889
            cable_user%limit_labile        = .true.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_casa_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
            cd ../outputs
            renameid ${rid} ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
            cd ..
        fi
        cd ${pdir}
    done
fi


# --------------------------------------------------------------------
# 4. Biomass into quasi-equilibrium without restricted N and P pools
if [[ ${doequi2} -eq 1 ]] ; then
    echo "4. Bring biomass into quasi-equilibrium without restricted N and P pools"
    for ((iequi2=1; iequi2<=${nequi2}; iequi2++)) ; do
        # 4a. 30 year run starting from restart files
        echo "   4a. 30 year spinup from accumulated biomass; iequi2=${iequi2}/${nequi2}"
        #rid="spinup"
        rid="spinup_nutrient_limited${iequi2}"
        # Met forcing
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd = 1859
        # cable_user%CASA_SPIN_ENDYEAR = 1869
        #MCTEST
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1860
            cable_user%YearEnd             = 1889
            icycle                         = 2
            spincasa                       = .false.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .false.
            cable_user%CASA_DUMP_WRITE     = .true.
            cable_user%CASA_SPIN_STARTYEAR = 1860
            cable_user%CASA_SPIN_ENDYEAR   = 1869
            cable_user%limit_labile        = .false.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
        cd ..
        cd ${pdir}
        #
        # 4b. analytic quasi-equilibrium of biomass pools
        echo "   4b. Analytic solution of biomass pools"
        #rid="spinup_analytic"
        rid="spin_casa_nutrient_limited${iequi2}"
        # Met forcing
        cat > ${tmp}/sedtmp.${pid} << EOF
            Run = "spinup"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml
        # LUC
        cp ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
        # Cable
        #MCTEST
        # cable_user%YearEnd = 1859
        # cable_user%CASA_SPIN_ENDYEAR = 1859
        #MCTEST
        cat > ${tmp}/sedtmp.${pid} << EOF
            cable_user%CLIMATE_fromZero    = .false.
            cable_user%YearStart           = 1860
            cable_user%YearEnd             = 1889
            icycle                         = 12
            spincasa                       = .true.
            cable_user%CASA_fromZero       = .false.
            cable_user%CASA_DUMP_READ      = .true.
            cable_user%CASA_DUMP_WRITE     = .false.
            cable_user%CASA_SPIN_STARTYEAR = 1860
            cable_user%CASA_SPIN_ENDYEAR   = 1889
            cable_user%limit_labile        = .false.
            cable_user%POP_fromZero        = .false.
            cable_user%POP_out             = "ini"
            cable_user%POPLUC              = .true.
            cable_user%POPLUC_RunType      = "static"
EOF
        applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
        # run model
        cd ${rdir}
        irm logs/log_cable.txt logs/log_out_cable.txt
        if [[ ${dompi} -eq 1 ]] ; then
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        else
            ./${iexe} > logs/log_out_cable.txt
        fi
        # save output
        renameid ${rid} ${mettype}.nml LUC.nml cable.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} ${mettype}_casa_rst.nc pop_${mettype}_ini.nc
        copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc
        if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
            cd ../outputs
            renameid ${rid} ${mettype}_out_casa.nc ${mettype}_out_casa_c13o2.nc
            cd ..
        fi
        cd ${pdir}
    done
fi

# --------------------------------------------------------------------
# 5a. First dynamic land use
if [[ ${doiniLUC} -eq 1 ]] ; then
    echo "5a. First dynamic land use"
    # Met forcing
    YearStart=1580
    YearEnd=1699
    cat > ${tmp}/sedtmp.${pid} << EOF
         Run = "spinup"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml

    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = ${YearStart}
         YearEnd   = ${YearEnd}
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml
    # Cable
    #MCTEST
    # cable_user%CASA_SPIN_ENDYEAR = 1859
    # cable_user%YearEnd             = $(( ${YearStart} + 1 ))
    #MCTEST
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero     = .false.
        cable_user%YearStart            = ${YearStart}
        cable_user%YearEnd              = ${YearEnd}
        icycle                          = 12
        spincasa                        = .false.
        cable_user%CASA_OUT_FREQ        = "annually"
        cable_user%CASA_fromZero        = .false.
        cable_user%CASA_DUMP_READ       = .true.
        cable_user%CASA_DUMP_WRITE      = .false.
        cable_user%CASA_SPIN_STARTYEAR  = 1860
        cable_user%CASA_SPIN_ENDYEAR    = 1889
        cable_user%limit_labile         = .false.
        cable_user%POP_fromZero         = .false.
        cable_user%POP_out              = "ini"
        cable_user%POPLUC               = .true.
        cable_user%POPLUC_RunType       = "init"
        cable_user%LUC_restart_in       = ""
        cable_user%c13o2_restart_in_LUC = ""
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_casa_rst.nc ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_LUC_rst.nc
    # cd ../outputs
    # renameid ${rid} ${mettype}_out_LUC.nc
    # cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 5b. Second full dynamic spinup
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "5b. Full dynamic spinup"
    # Met forcing
    YearStart=1700
    YearEnd=1899
    cat > ${tmp}/sedtmp.${pid} << EOF
         Run = "premet"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml

    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = ${YearStart}
         YearEnd   = ${YearEnd}
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml

    # [TODO] Need to change POPLUC and POPLUC_RunType to switch LUC on/off
    # Set POPLUC = .false. and POPLUC_RunType = 'static'  for no LUC
    # Set POPLUC = .true.  and POPLUC_RunType = 'restart' for LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1860
        cable_user%CASA_SPIN_ENDYEAR   = 1889
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "restart"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_LUC_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi


# --------------------------------------------------------------------
# 6. Final centennial run
if [[ ${dofinal} -eq 1 ]] ; then
    echo "6. Final centennial run"
    YearStart=1900
    YearEnd=2022
    cat > ${tmp}/sedtmp.${pid} << EOF
         Run = "standard"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/bios_${experiment}.nml ${rdir}/bios.nml

    rid=${YearStart}_${YearEnd}
    # LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
         YearStart = ${YearStart}
         YearEnd   = ${YearEnd}
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/LUC_${experiment}.nml ${rdir}/LUC.nml

    # [TODO] Need to change POPLUC and POPLUC_RunType to switch LUC on/off
    # Set POPLUC = .false. and POPLUC_RunType = 'static'  for no LUC
    # Set POPLUC = .true.  and POPLUC_RunType = 'restart' for LUC
    cat > ${tmp}/sedtmp.${pid} << EOF
        cable_user%CLIMATE_fromZero    = .false.
        cable_user%YearStart           = ${YearStart}
        cable_user%YearEnd             = ${YearEnd}
        icycle                         = 2
        spincasa                       = .false.
        cable_user%CASA_fromZero       = .false.
        cable_user%CASA_DUMP_READ      = .false.
        cable_user%CASA_DUMP_WRITE     = .false.
        cable_user%CASA_SPIN_STARTYEAR = 1860
        cable_user%CASA_SPIN_ENDYEAR   = 1889
        cable_user%limit_labile        = .false.
        cable_user%POP_fromZero        = .false.
        cable_user%POP_out             = "ini"
        cable_user%POPLUC              = .true.
        cable_user%POPLUC_RunType      = "restart"
EOF
    applysed ${tmp}/sedtmp.${pid} ${rdir}/cable_${experiment}.nml ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} ${mettype}.nml LUC.nml cable.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    copyid ${rid} ${mettype}_LUC_rst.nc pop_${mettype}_ini.nc
    copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_LUC_rst.nc
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    renameid ${rid} ${mettype}_out_casa_c13o2.nc
    cd ..
    cd ${pdir}
fi

# --------------------------------------------------------------------
# Finish
#
cd ${isdir}

t2=$(date +%s)
dt=$((t2-t1))
printf "\n"
if [[ ${dt} -lt 60 ]] ; then
    printf "Finished at %s   in %i seconds.\n" "$(date)" ${dt}
else
    dm=$(echo "(${t2}-${t1})/60." | bc -l)
    printf "Finished at %s   in %.2f minutes.\n" "$(date)" ${dm}
fi

exit
