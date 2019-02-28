$(message "GET INTO makefile")
simulator:
	$(EXE) prog=simulator simulation=${SIMULATION} frame=${FRAME} dt=${DT} line_search=${LINE_SEARCH} weight_line_search=${WEIGHT_LINE_SEARCH} gravity=${GRAVITY} density=${DENSITY} stiffness=${STIFFNESS} newton_fastMS=${NEWTON_FASTMS} out_dir_simulator=${OUT_DIR_SIMULATOR} object_name=${OBJECT_NAME} input_object=${INPUT_OBJECT} input_constraint=${INPUT_CONSTRAINT} force_function=${FORCE_FUNCTION} radius=${RADIUS} intensity=${INTENSITY}
object_creator:
	$(EXE) prog=object_creator l_size=${L_SIZE} w_size=${W_SIZE} h_size=${H_SIZE} dmetric=${DMETRIC} object_name=${OBJECT_NAME} out_dir=${OUT_DIR} dim=${DIM}
