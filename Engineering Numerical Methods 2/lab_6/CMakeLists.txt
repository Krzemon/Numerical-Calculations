cmake_minimum_required(VERSION 3.10)
project(lab_06)
set(PROJECT_AUTHOR "Przemysław Ryś")

# --- Ustawienia kompilacji ---
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")

# --- Katalogi projektu ---
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(DATA_DIR ${PROJECT_SOURCE_DIR}/data)
set(PLOTS_DIR ${PROJECT_SOURCE_DIR}/plots)
set(PYTHON_DIR ${PROJECT_SOURCE_DIR}/python)

# --- Pliki wynikowe ---
set(DATA_FILE_ZAD1 ${DATA_DIR}/zad_1.dat)
set(DATA_FILE_ZAD2_E ${DATA_DIR}/zad_2_global_E.dat)
set(DATA_FILE_ZAD2_O ${DATA_DIR}/zad_2_global_O.dat)

set(DATA_FILE_ZAD2_MODE_0 ${DATA_DIR}/mode_0.dat)
set(DATA_FILE_ZAD2_MODE_1 ${DATA_DIR}/mode_1.dat)
set(DATA_FILE_ZAD2_MODE_2 ${DATA_DIR}/mode_2.dat)
set(DATA_FILE_ZAD2_MODE_3 ${DATA_DIR}/mode_3.dat)
set(DATA_FILE_ZAD2_MODE_4 ${DATA_DIR}/mode_4.dat)
set(DATA_FILE_ZAD2_MODE_5 ${DATA_DIR}/mode_5.dat)
set(DATA_FILE_ZAD2_MODE_6 ${DATA_DIR}/mode_6.dat)
set(DATA_FILE_ZAD2_MODE_7 ${DATA_DIR}/mode_7.dat)
set(DATA_FILE_ZAD2_MODE_8 ${DATA_DIR}/mode_8.dat)
set(DATA_FILE_ZAD2_MODE_9 ${DATA_DIR}/mode_9.dat)

# set(DATA_FILE_ZAD2_MODES ${DATA_DIR}/zad_2_modes.dat)
# set(DATA_FILE_ZAD3 ${DATA_DIR}/zad_3.dat)
# set(DATA_FILE_ZAD6 ${DATA_DIR}/zad_6.dat)

set(PLOT_FILE_ZAD1 ${PLOTS_DIR}/zad_1_grid.png)
set(PLOT_FILE_ZAD2_E ${PLOTS_DIR}/zad_2_E_heatmap.png)
set(PLOT_FILE_ZAD2_O ${PLOTS_DIR}/zad_2_O_heatmap.png)
set(PLOT_FILE_ZAD2_MODE_0 ${PLOTS_DIR}/zad_2_mode_0.png)
set(PLOT_FILE_ZAD2_MODE_1 ${PLOTS_DIR}/zad_2_mode_1.png)
set(PLOT_FILE_ZAD2_MODE_2 ${PLOTS_DIR}/zad_2_mode_2.png)
set(PLOT_FILE_ZAD2_MODE_3 ${PLOTS_DIR}/zad_2_mode_3.png)
set(PLOT_FILE_ZAD2_MODE_4 ${PLOTS_DIR}/zad_2_mode_4.png)
set(PLOT_FILE_ZAD2_MODE_5 ${PLOTS_DIR}/zad_2_mode_5.png)
set(PLOT_FILE_ZAD2_MODE_6 ${PLOTS_DIR}/zad_2_mode_6.png)
set(PLOT_FILE_ZAD2_MODE_7 ${PLOTS_DIR}/zad_2_mode_7.png)
set(PLOT_FILE_ZAD2_MODE_8 ${PLOTS_DIR}/zad_2_mode_8.png)
set(PLOT_FILE_ZAD2_MODE_9 ${PLOTS_DIR}/zad_2_mode_9.png)

# set(PLOT_FILE_ZAD2_MODES ${PLOTS_DIR}/zad_2_modes.png)
# set(PLOT_FILE_ZAD3 ${PLOTS_DIR}/zad_3_grid.png)
# set(PLOT_FILE_ZAD6 ${PLOTS_DIR}/zad_6_grid.png)

# --- Skrypty wykresow ---
set(PLOT_SCRIPT_GRID ${PYTHON_DIR}/plot_grid.py)
set(PLOT_SCRIPT_HEATMAP ${PYTHON_DIR}/plot_heatmap.py)
set(PLOT_SCRIPT_MODES ${PYTHON_DIR}/plot_modes.py)

# set(PLOT_SCRIPT_ZAD3 ${PYTHON_DIR}/plot_zad3.py)
# set(PLOT_SCRIPT_ZAD6 ${PYTHON_DIR}/plot_zad6.py)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

find_package(LAPACK REQUIRED)
# --- Linkowanie ---
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE gsl gslcblas m ${LAPACK_LIBRARIES})
target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# --- Po kompilacji: uruchomienie programu i generowanie wykresu ---
add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}
    COMMAND python3 ${PLOT_SCRIPT_GRID} ${DATA_FILE_ZAD1} ${PLOT_FILE_ZAD1} "Siatka węzłów"
    COMMAND python3 ${PLOT_SCRIPT_HEATMAP} ${DATA_FILE_ZAD2_E} ${PLOT_FILE_ZAD2_E} "Mapa ciepła macierzy sztywności E"
    COMMAND python3 ${PLOT_SCRIPT_HEATMAP} ${DATA_FILE_ZAD2_O} ${PLOT_FILE_ZAD2_O} "Mapa ciepła macierzy całek przekrywania O"

    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_0} ${PLOT_FILE_ZAD2_MODE_0} "mod 0"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_1} ${PLOT_FILE_ZAD2_MODE_1} "mod 1"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_2} ${PLOT_FILE_ZAD2_MODE_2} "mod 2"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_3} ${PLOT_FILE_ZAD2_MODE_3} "mod 3"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_4} ${PLOT_FILE_ZAD2_MODE_4} "mod 4"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_5} ${PLOT_FILE_ZAD2_MODE_5} "mod 5"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_6} ${PLOT_FILE_ZAD2_MODE_6} "mod 6"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_7} ${PLOT_FILE_ZAD2_MODE_7} "mod 7"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_8} ${PLOT_FILE_ZAD2_MODE_8} "mod 8"
    COMMAND python3 ${PLOT_SCRIPT_MODES} ${DATA_FILE_ZAD2_MODE_9} ${PLOT_FILE_ZAD2_MODE_9} "mod 9"


    # COMMAND python3 ${PLOT_SCRIPT_ZAD3} ${DATA_FILE_ZAD3} ${PLOT_FILE_ZAD3} "Drugi mod własny"
    # COMMAND python3 ${PLOT_SCRIPT_ZAD6} ${DATA_FILE_ZAD6} ${PLOT_FILE_ZAD6} "Mapa u x,y"
)
