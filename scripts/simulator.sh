#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} simulator EXE=${exe} SIMULATION=${simulation} FRAME=${frame} DT=${dt} LINE_SEARCH=${line_search} WEIGHT_LINE_SEARCH=${weight_line_search} GRAVITY=${gravity} DENSITY=${density} STIFFNESS=${stiffness} NEWTON_FASTMS=${newton_fastMS} OUT_DIR_SIMULATOR=${out_dir_simulator} OBJECT_NAME=${object_name} INPUT_OBJECT=${input_object} INPUT_CONSTRAINT=${input_constraint} FORCE_FUNCTION=${force_function} RADIUS=${radius} INTENSITY=${intensity}
