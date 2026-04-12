cmake_minimum_required(VERSION 3.10)
project(nano_mezo_1 LANGUAGES CXX)
set(PROJECT_AUTHOR "Przemysław Ryś")

# --- Ustawienia kompilacji ---
set(CMAKE_BUILD_TYPE Debug) # Release
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

# --- DEBUG ---
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# --- RELEASE ---
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

# --- Katalogi ---
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(DATA_DIR ${PROJECT_SOURCE_DIR}/data)
set(PLOTS_DIR ${PROJECT_SOURCE_DIR}/plots)
set(PYTHON_DIR ${PROJECT_SOURCE_DIR}/python)

# --- Upewnij się, że istnieją ---
# file(MAKE_DIRECTORY ${DATA_DIR})
# file(MAKE_DIRECTORY ${PLOTS_DIR})
# file(MAKE_DIRECTORY ${PYTHON_DIR})

# --- Skrypty ---
set(PLOT_SCRIPT_CONTOUR ${PYTHON_DIR}/plot_contour.py)
set(PLOT_SCRIPT_ENERGIES ${PYTHON_DIR}/plot_energies.py)
set(PLOT_SCRIPT_MATRIX ${PYTHON_DIR}/plot_matrix.py)

# --- Pliki danych ---
set(DATA_BASIS_0  ${DATA_DIR}/basis_k0.dat)
set(DATA_BASIS_8  ${DATA_DIR}/basis_k8.dat)
set(DATA_BASIS_9  ${DATA_DIR}/basis_k9.dat)

set(DATA_MATRIX_H ${DATA_DIR}/matrix_H.dat)
set(DATA_MATRIX_S ${DATA_DIR}/matrix_S.dat)
# set(DATA_ENERGIES ${DATA_DIR}/energies.dat)
set(DATA_WF_0     ${DATA_DIR}/WF_state0.dat)
set(DATA_WF_1     ${DATA_DIR}/WF_state1.dat)
set(DATA_WF_2     ${DATA_DIR}/WF_state2.dat)
set(DATA_WF_3     ${DATA_DIR}/WF_state3.dat)
set(DATA_WF_4     ${DATA_DIR}/WF_state4.dat)
set(DATA_WF_5     ${DATA_DIR}/WF_state5.dat)
set(DATA_ENERGIES_OMEGA ${DATA_DIR}/energies_omega.dat)

set(DATA_WF_0_zad5     ${DATA_DIR}/WF_state0_wy.dat)
set(DATA_WF_1_zad5     ${DATA_DIR}/WF_state1_wy.dat)
set(DATA_WF_2_zad5     ${DATA_DIR}/WF_state2_wy.dat)
set(DATA_WF_3_zad5     ${DATA_DIR}/WF_state3_wy.dat)
set(DATA_WF_4_zad5     ${DATA_DIR}/WF_state4_wy.dat)
set(DATA_WF_5_zad5     ${DATA_DIR}/WF_state5_wy.dat)

# --- Pliki wykresów ---
set(PLOT_BASIS_0  ${PLOTS_DIR}/zad_1/basis_k0.png)
set(PLOT_BASIS_8  ${PLOTS_DIR}/zad_1/basis_k8.png)
set(PLOT_BASIS_9  ${PLOTS_DIR}/zad_1/basis_k9.png)

set(PLOT_MATRIX_H  ${PLOTS_DIR}/zad_2/matrix_H.png)
set(PLOT_MATRIX_S  ${PLOTS_DIR}/zad_2/matrix_S.png)
# set(PLOT_ENERGIES ${PLOTS_DIR}/zad_2/energies.png)
set(PLOT_WF_0     ${PLOTS_DIR}/zad_3/WF_state0.png)
set(PLOT_WF_1     ${PLOTS_DIR}/zad_3/WF_state1.png)
set(PLOT_WF_2     ${PLOTS_DIR}/zad_3/WF_state2.png)
set(PLOT_WF_3     ${PLOTS_DIR}/zad_3/WF_state3.png)
set(PLOT_WF_4     ${PLOTS_DIR}/zad_3/WF_state4.png)
set(PLOT_WF_5     ${PLOTS_DIR}/zad_3/WF_state5.png)

set(PLOT_ENERGIES_OMEGA ${PLOTS_DIR}/zad_4/energies_omega.png)
set(PLOT_WF_0_zad5     ${PLOTS_DIR}/zad_5/WF_state0_wy.png)
set(PLOT_WF_1_zad5     ${PLOTS_DIR}/zad_5/WF_state1_wy.png)
set(PLOT_WF_2_zad5     ${PLOTS_DIR}/zad_5/WF_state2_wy.png)
set(PLOT_WF_3_zad5     ${PLOTS_DIR}/zad_5/WF_state3_wy.png)
set(PLOT_WF_4_zad5     ${PLOTS_DIR}/zad_5/WF_state4_wy.png)
set(PLOT_WF_5_zad5     ${PLOTS_DIR}/zad_5/WF_state5_wy.png)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

## --- Eigen (header-only) ---
find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})

# --- GSL/LAPACK/BLAAS niepotrzebne (zakomentowane) ---
# find_package(GSL REQUIRED)
# target_link_libraries(${PROJECT_NAME} PRIVATE gsl gslcblas m)
# target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
# find_package(LAPACK REQUIRED)

# [CXX] std::string path = std::string(PROJECT_ROOT_DIR) + "/data/file.dat";
target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# --- Po kompilacji: uruchomienie programu i generowanie wykresu ---
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}

    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_BASIS_0}"  "${PLOT_BASIS_0}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_BASIS_8}"  "${PLOT_BASIS_8}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_BASIS_9}"  "${PLOT_BASIS_9}"

    COMMAND python3 "${PLOT_SCRIPT_MATRIX}" "${DATA_MATRIX_H}"  "${PLOT_MATRIX_H}"
    COMMAND python3 "${PLOT_SCRIPT_MATRIX}" "${DATA_MATRIX_S}"  "${PLOT_MATRIX_S}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_0}"     "${PLOT_WF_0}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_1}"     "${PLOT_WF_1}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_2}"     "${PLOT_WF_2}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_3}"     "${PLOT_WF_3}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_4}"     "${PLOT_WF_4}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_5}"     "${PLOT_WF_5}"
    COMMAND python3 "${PLOT_SCRIPT_ENERGIES}" "${DATA_ENERGIES_OMEGA}" "${PLOT_ENERGIES_OMEGA}"

    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_0_zad5}" "${PLOT_WF_0_zad5}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_1_zad5}" "${PLOT_WF_1_zad5}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_2_zad5}" "${PLOT_WF_2_zad5}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_3_zad5}" "${PLOT_WF_3_zad5}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_4_zad5}" "${PLOT_WF_4_zad5}"
    COMMAND python3 "${PLOT_SCRIPT_CONTOUR}" "${DATA_WF_5_zad5}" "${PLOT_WF_5_zad5}"
)