INPUT=(
    '/g/data/rp23/experiments/2024-04-17_BIOS3-merge/sb8430/params_bios/def_soil_params.txt'
    '/g/data/rp23/experiments/2024-04-17_BIOS3-merge/sb8430/params_bios/def_veg_params.txt'
    '/g/data/rp23/experiments/2024-04-17_BIOS3-merge/sb8430/params_bios/pftlookup.csv'
)

# Use shell parameter expansion for runtime configurability - see
# https://www.gnu.org/software/bash/manual/html_node/Shell-Parameter-Expansion.html:
: "${EXE:=/home/189/sb8430/cable/offline/cable-mpi}"
: "${RESTART:=/g/data/rp23/experiments/2024-04-17_BIOS3-merge/sb8430/bios_act9test/restart_1700_1899}"
