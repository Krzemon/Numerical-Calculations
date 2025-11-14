cmake_minimum_required(VERSION 3.10)
project(lab_03)
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
set(DATA_FILE_ZAD1 ${DATA_DIR}/zad_1_matrix.dat)
set(DATA_FILE_ZAD2E_M5 ${DATA_DIR}/zad_2_E_M5.dat)
set(DATA_FILE_ZAD2U_M5 ${DATA_DIR}/zad_2_psi_M5_alpha_1_4.dat)
set(DATA_FILE_ZAD3E_M10 ${DATA_DIR}/zad_3_E_M10.dat)
set(DATA_FILE_ZAD3E_M30 ${DATA_DIR}/zad_3_E_M30.dat)
set(DATA_FILE_ZAD3U_M10 ${DATA_DIR}/zad_3_psi_M10_alpha_1_4.dat)
set(DATA_FILE_ZAD3U_M30 ${DATA_DIR}/zad_3_psi_M30_alpha_1_4.dat)

set(PLOT_FILE_ZAD1 ${PLOTS_DIR}/zad_1_heatmap.png)
set(PLOT_FILE_ZAD2E_M5 ${PLOTS_DIR}/zad_2_E_M5.png)
set(PLOT_FILE_ZAD2U_M5 ${PLOTS_DIR}/zad_2_psi_M5.png)
set(PLOT_FILE_ZAD3E_M10 ${PLOTS_DIR}/zad_3_E_M10.png)
set(PLOT_FILE_ZAD3E_M30 ${PLOTS_DIR}/zad_3_E_M30.png)
set(PLOT_FILE_ZAD3U_M10 ${PLOTS_DIR}/zad_3_psi_M10.png)
set(PLOT_FILE_ZAD3U_M30 ${PLOTS_DIR}/zad_3_psi_M30.png)

# --- Skrypty wykresow ---
set(PLOT_SCRIPT_ZAD1 ${PYTHON_DIR}/plot_zad1.py)
set(PLOT_SCRIPT_ZAD2E_M5 ${PYTHON_DIR}/plot_zad2_E_M5.py)
set(PLOT_SCRIPT_ZAD2U_M5 ${PYTHON_DIR}/plot_zad2_psi_M5.py)
set(PLOT_SCRIPT_ZAD3E_M10 ${PYTHON_DIR}/plot_zad3_E_M10.py)
set(PLOT_SCRIPT_ZAD3E_M30 ${PYTHON_DIR}/plot_zad3_E_M30.py)
set(PLOT_SCRIPT_ZAD3U_M10 ${PYTHON_DIR}/plot_zad3_psi_M10.py)
set(PLOT_SCRIPT_ZAD3U_M30 ${PYTHON_DIR}/plot_zad3_psi_M30.py)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Linkowanie ---
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE gsl gslcblas m)
target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# --- Po kompilacji: uruchomienie programu i generowanie wykresu ---
add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}
    COMMAND python3 ${PLOT_SCRIPT_ZAD1} ${DATA_FILE_ZAD1} ${PLOT_FILE_ZAD1} "Mapa ciepła macierzy S"
    COMMAND python3 ${PLOT_SCRIPT_ZAD2E_M5} ${DATA_FILE_ZAD2E_M5} ${PLOT_FILE_ZAD2E_M5} "E_mu_alpha, M=5"
    COMMAND python3 ${PLOT_SCRIPT_ZAD2U_M5} ${DATA_FILE_ZAD2U_M5} ${PLOT_FILE_ZAD2U_M5} "Funkcje własne u_mu_x, M=5"
    COMMAND python3 ${PLOT_SCRIPT_ZAD3E_M10} ${DATA_FILE_ZAD3E_M10} ${PLOT_FILE_ZAD3E_M10} "E_mu_alpha, M=10"
    COMMAND python3 ${PLOT_SCRIPT_ZAD3E_M30} ${DATA_FILE_ZAD3E_M30} ${PLOT_FILE_ZAD3E_M30} "E_mu_alpha, M=30"
    COMMAND python3 ${PLOT_SCRIPT_ZAD3U_M10} ${DATA_FILE_ZAD3U_M10} ${PLOT_FILE_ZAD3U_M10} "Funkcje własne u_mu_x, M=10"
    COMMAND python3 ${PLOT_SCRIPT_ZAD3U_M30} ${DATA_FILE_ZAD3U_M30} ${PLOT_FILE_ZAD3U_M30} "Funkcje własne u_mu_x, M=30"
)