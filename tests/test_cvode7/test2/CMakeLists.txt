#
# CMakeLists.txt
#
# If the SUNDIALS library is not installed in /usr or /usr/local,
# set CMAKE_INCLUDE_PATH and CMAKE_LIBRARY_PATH to the appropriate
# paths.
#

cmake_minimum_required(VERSION 3.22.1)
project(test_linearoscp_cv7 C)
add_executable(test_linearoscp_cv7 test_linearoscp_cv7.c linearoscp_cv7.c)

find_path(CVODE_INCLUDE_DIR cvode/cvode.h /usr/local/include /usr/include)
find_path(NVECTOR_SERIAL_INCLUDE_DIR nvector/nvector_serial.h /usr/local/include /usr/include)
find_path(SUNMATRIX_DENSE_INCLUDE_DIR sunmatrix/sunmatrix_dense.h /usr/local/include /usr/include)
find_path(SUNLINSOL_DENSE_INCLUDE_DIR sunlinsol/sunlinsol_dense.h /usr/local/include /usr/include)

find_library(SUNDIALS_CVODE_LIBRARY sundials_cvode /usr/local/lib /usr/lib)
find_library(SUNDIALS_NVECSERIAL_LIBRARY sundials_nvecserial /usr/local/lib /usr/lib)
find_library(SUNDIALS_CORE_LIBRARY sundials_core /usr/local/lib /usr/lib)
find_library(SUNDIALS_SUNMATRIXDENSE_LIBRARY sundials_sunmatrixdense /usr/local/lib /usr/lib)

include_directories(${CVODE_INCLUDE_DIR} ${NVECTOR_SERIAL_INCLUDE_DIR} ${SUNMATRIX_DENSE_INCLUDE_DIR} ${SUNLINSOL_DENSE_INCLUDE_DIR})

target_link_libraries(test_linearoscp_cv7 ${SUNDIALS_CVODE_LIBRARY} ${SUNDIALS_NVECSERIAL_LIBRARY} ${SUNDIALS_CORE_LIBRARY} ${SUNDIALS_SUNMATRIXDENSE_LIBRARY} -lm)
