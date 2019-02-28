#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} object_creator EXE=${exe} L_SIZE=${l_size} W_SIZE=${w_size} H_SIZE=${h_size} DMETRIC=${dmetric} OBJECT_NAME=${object_name} OUT_DIR=${out_dir} DIM=${dim}

