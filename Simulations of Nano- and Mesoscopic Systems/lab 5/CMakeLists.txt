cmake_minimum_required(VERSION 3.10)
project(nano_mezo_5 LANGUAGES CXX)
set(PROJECT_AUTHOR "Przemysław Ryś")

# --- Ustawienia kompilacji ---
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

# --- OpenMP ---
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# --- Katalogi ---
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(DATA_DIR ${PROJECT_SOURCE_DIR}/data)
set(PLOTS_DIR ${PROJECT_SOURCE_DIR}/plots)
set(PYTHON_DIR ${PROJECT_SOURCE_DIR}/python)

file(MAKE_DIRECTORY ${DATA_DIR})
file(MAKE_DIRECTORY ${PLOTS_DIR})

# --- Skrypty Pythona ---
set(PLOT_SCRIPT_ZAD1 ${PYTHON_DIR}/plot_zad1.py)
set(PLOT_SCRIPT_ZAD2 ${PYTHON_DIR}/plot_zad2.py)
set(PLOT_SCRIPT_ZAD3 ${PYTHON_DIR}/plot_zad3.py)
set(PLOT_SCRIPT_ZAD4 ${PYTHON_DIR}/plot_zad4.py)
set(PLOT_SCRIPT_ZAD5 ${PYTHON_DIR}/plot_zad5.py)

# --- Pliki danych ---
set(DATA_ZAD1_DISP ${DATA_DIR}/zad1_dispersion.dat)
set(DATA_ZAD2_U_E1 ${DATA_DIR}/zad2_u_E1_)
set(DATA_ZAD2_U_E2 ${DATA_DIR}/zad2_u_E2_)
set(DATA_ZAD3_WF1 ${DATA_DIR}/zad3_wf_200.dat)
set(DATA_ZAD3_WF2 ${DATA_DIR}/zad3_wf_400.dat)
set(DATA_ZAD4_COND ${DATA_DIR}/zad4_conductance.dat)
set(DATA_ZAD4_POT ${DATA_DIR}/zad4_potential.dat)
set(DATA_ZAD5_WF1 ${DATA_DIR}/zad5_wf_-130.dat)
set(DATA_ZAD5_WF2 ${DATA_DIR}/zad5_wf_-110.dat)
set(DATA_ZAD5_WF3 ${DATA_DIR}/zad5_wf_-85.dat)

# --- Pliki wykresów ---
set(PLOT_ZAD1_DISP ${PLOTS_DIR}/zad_1/dispersion.png)
set(PLOT_ZAD1_DISP_E1 ${PLOTS_DIR}/zad_1/dispersion_E1.png)
set(PLOT_ZAD1_DISP_E2 ${PLOTS_DIR}/zad_1/dispersion_E2.png)
set(PLOT_ZAD2_E1_MODES ${PLOTS_DIR}/zad_2/modes_E1.png)
set(PLOT_ZAD2_E2_MODES ${PLOTS_DIR}/zad_2/modes_E2.png)
set(PLOT_ZAD3_WF1 ${PLOTS_DIR}/zad_3/wf_E02.png)
set(PLOT_ZAD3_WF2 ${PLOTS_DIR}/zad_3/wf_E04.png)
set(PLOT_ZAD4_T ${PLOTS_DIR}/zad_4/conductance_T.png)
set(PLOT_ZAD4_R ${PLOTS_DIR}/zad_4/conductance_R.png)
set(PLOT_ZAD4_TR ${PLOTS_DIR}/zad_4/conductance_TR.png)
set(PLOT_ZAD4_POT ${PLOTS_DIR}/zad_4/potential.png)
set(PLOT_ZAD5_WF1 ${PLOTS_DIR}/zad_5/wf_G0.png)
set(PLOT_ZAD5_WF2 ${PLOTS_DIR}/zad_5/wf_G1.png)
set(PLOT_ZAD5_WF3 ${PLOTS_DIR}/zad_5/wf_G2.png)

# --- Źródła ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Eigen (header-only) ---
find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen m)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})

target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# --- Po kompilacji: uruchomienie programu i generowanie wykresów ---
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}
    COMMAND python3 "${PLOT_SCRIPT_ZAD1}" "${DATA_ZAD1_DISP}" "${PLOT_ZAD1_DISP}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD2}" "${DATA_ZAD2_U_E1}" "${DATA_DIR}/zad2_modes_E1.dat" "${DATA_ZAD1_DISP}" "0.2" "${PLOTS_DIR}/zad_2/E1"
    COMMAND python3 "${PLOT_SCRIPT_ZAD2}" "${DATA_ZAD2_U_E2}" "${DATA_DIR}/zad2_modes_E2.dat" "${DATA_ZAD1_DISP}" "0.4" "${PLOTS_DIR}/zad_2/E2"
    COMMAND python3 "${PLOT_SCRIPT_ZAD3}" "${DATA_ZAD3_WF1}" "0.2" "${PLOT_ZAD3_WF1}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD3}" "${DATA_ZAD3_WF2}" "0.4" "${PLOT_ZAD3_WF2}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD4}" "${DATA_ZAD4_COND}" "${DATA_ZAD4_POT}" "${PLOTS_DIR}/zad_4/conductance" "${PLOT_ZAD4_POT}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD5}" "${DATA_ZAD5_WF1}" "-1.3" "0" "${PLOT_ZAD5_WF1}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD5}" "${DATA_ZAD5_WF2}" "-1.1" "G_0" "${PLOT_ZAD5_WF2}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD5}" "${DATA_ZAD5_WF3}" "-0.85" "2G_0" "${PLOT_ZAD5_WF3}"
    COMMENT "Uruchamianie symulacji i generowanie wykresów"
)
