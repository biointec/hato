cmake_minimum_required(VERSION 3.0)

project(hato)

# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set the CPLEX root directory
set(CPLEX_ROOT_DIR  "../../cplex/CPLEX_Studio2211/"  CACHE PATH "Path to root Cplex installation")

# Set the CPLEX installation directory
set(CPLEX_DIR  "../../cplex/CPLEX_Studio2211/cplex"  CACHE PATH "Path to Cplex installation")


# Find CPLEX package configuration files
find_package(CPLEX REQUIRED)

# Include CPLEX headers
include_directories(${CPLEX_INCLUDE_DIRS})

# Link against CPLEX libraries
add_executable(hato hato.cpp generator.cpp cplexmodel.cpp )
target_link_libraries(hato ${CPLEX_LIBRARIES} dl)
install(TARGETS hato DESTINATION bin)
