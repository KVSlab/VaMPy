#!/bin/bash
CASES=(0003)
CONDITION=SR
for CASE in "${CASES[@]}"; do
  echo "Conevrting case ${CASE} , condiiton ${CONDITION}"
  pvpython combine_and_convert_energy_to_vtu.py --case $CASE --condition $CONDITION
  pvpython merge_energy_vtu.py --case $CASE --condition $CONDITION
  pvpython combine_and_convert_hemo_to_vtp.py --case $CASE --condition $CONDITION
  pvpython convert_brt_to_vtu.py --case $CASE --condition $CONDITION
done
