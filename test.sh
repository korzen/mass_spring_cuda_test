#! /bin/bash

<<EOF
#parameter for solvers

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic  
# if dynamic ,the num of frame
frame=500
# if dynamic ,the dt is needed
dt=0.01  # belong to object
#gravity
gravity=9.8
#density
density=10 # belong to object
#############common para################
# whether it uses line_search
line_search=1  #belong to object
# if using line_search, the weight
weight_line_search=1e-5  #belong to object
stiffness=600
newton_fastMS=newton
force_function=gravity

radius=0.11
intensity=100
##########################################################
project_dir=./
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/fast_mass_spring
base_output_dir=${project_dir}/example
##########################################################
object_name=cloth
out_dir_simulator=${base_output_dir}/${object_name}/${newton_fastMS}/
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk
input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv
mkdir -p ${out_dir_simulator}
source ${script_dir}/simulator.sh



#parameter for solvers

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic  
# if dynamic ,the num of frame
frame=500
# if dynamic ,the dt is needed
dt=0.01  # belong to object
#gravity
gravity=9.8
#density
density=10 # belong to object
#############common para################
# whether it uses line_search
line_search=0  #belong to object
# if using line_search, the weight
weight_line_search=1e-5  #belong to object
stiffness=600
newton_fastMS=fastMS_original
force_function=gravity

radius=0.11
intensity=100
##########################################################
project_dir=./
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/fast_mass_spring
base_output_dir=${project_dir}/example
##########################################################
object_name=cloth
out_dir_simulator=${base_output_dir}/${object_name}/${newton_fastMS}/
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk
input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv
mkdir -p ${out_dir_simulator}
source ${script_dir}/simulator.sh


#parameter for solvers

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic  
# if dynamic ,the num of frame
frame=500
# if dynamic ,the dt is needed
dt=0.01  # belong to object
#gravity
gravity=9.8
#density
density=10 # belong to object
#############common para################
# whether it uses line_search
line_search=0  #belong to object
# if using line_search, the weight
weight_line_search=1e-5  #belong to object
stiffness=600
newton_fastMS=fastMS_ChebyshevSIM
force_function=gravity

radius=0.11
intensity=100
##########################################################
project_dir=./
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/fast_mass_spring
base_output_dir=${project_dir}/example
##########################################################
object_name=cloth
out_dir_simulator=${base_output_dir}/${object_name}/${newton_fastMS}/
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk
input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv
mkdir -p ${out_dir_simulator}
source ${script_dir}/simulator.sh
EOF





<<EOF


EOF


#parameter for solvers

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic  
# if dynamic ,the num of frame
frame=10
# if dynamic ,the dt is needed
dt=0.01  # belong to object
#gravity
gravity=9.8
#density
density=10 # belong to object
#############common para################
# whether it uses line_search
line_search=1  #belong to object
# if using line_search, the weight
weight_line_search=1e-5  #belong to object
stiffness=600
newton_fastMS=newton
force_function=gravity

radius=0.11
intensity=0
##########################################################
project_dir=./
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/fast_mass_spring
base_output_dir=${project_dir}/example
##########################################################
object_name=cloth
out_dir_simulator=${base_output_dir}/${object_name}/${newton_fastMS}_1/
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk
input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv
mkdir -p ${out_dir_simulator}
source ${script_dir}/simulator.sh



#parameter for solvers

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic  
# if dynamic ,the num of frame
frame=1
# if dynamic ,the dt is needed
dt=0.01  # belong to object
#gravity
gravity=9.8
#density
density=10 # belong to object
#############common para################
# whether it uses line_search
line_search=0  #belong to object
# if using line_search, the weight
weight_line_search=1e-5  #belong to object
stiffness=600
newton_fastMS=fastMS_original
force_function=gravity

radius=0.11
intensity=0
##########################################################
project_dir=./
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/fast_mass_spring
base_output_dir=${project_dir}/example
##########################################################
object_name=cloth
out_dir_simulator=${base_output_dir}/${object_name}/${newton_fastMS}_1/
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk
input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv
mkdir -p ${out_dir_simulator}
source ${script_dir}/simulator.sh



#parameter for solvers

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic  
# if dynamic ,the num of frame
frame=1
# if dynamic ,the dt is needed
dt=0.01  # belong to object
#gravity
gravity=9.8
#density
density=10 # belong to object
#############common para################
# whether it uses line_search
line_search=0  #belong to object
# if using line_search, the weight
weight_line_search=1e-5  #belong to object
stiffness=600
newton_fastMS=fastMS_ChebyshevSIM
force_function=gravity

radius=0.11
intensity=0
##########################################################
project_dir=./
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/fast_mass_spring
base_output_dir=${project_dir}/example
##########################################################
object_name=cloth
out_dir_simulator=${base_output_dir}/${object_name}/${newton_fastMS}_1/
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk
input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv
mkdir -p ${out_dir_simulator}
source ${script_dir}/simulator.sh
