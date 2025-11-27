cmake_minimum_required(VERSION 3.10)
project(lab_04)
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
set(DATA_FILE_ZAD2 ${DATA_DIR}/zad_2.dat)
set(DATA_FILE_ZAD3 ${DATA_DIR}/zad_3.dat)
set(DATA_FILE_ZAD4 ${DATA_DIR}/zad_4.dat)
set(DATA_FILE_ZAD8_3 ${DATA_DIR}/zad_8_3.dat)
set(DATA_FILE_ZAD8_10 ${DATA_DIR}/zad_8_10.dat)

set(PLOT_FILE_ZAD2 ${PLOTS_DIR}/zad_2_nodes.png)
set(PLOT_FILE_ZAD3 ${PLOTS_DIR}/zad_3_S_matrix.png)
set(PLOT_FILE_ZAD4 ${PLOTS_DIR}/zad_4_S_matrix.png)
set(PLOT_FILE_ZAD8_3 ${PLOTS_DIR}/zad_8_n3_solution.png)
set(PLOT_FILE_ZAD8_10 ${PLOTS_DIR}/zad_8_n10_solution.png)

# --- Skrypty wykresow ---
set(PLOT_SCRIPT_ZAD2 ${PYTHON_DIR}/plot_zad2.py)
set(PLOT_SCRIPT_HEATMAP ${PYTHON_DIR}/plot_heatmap.py)
set(PLOT_SCRIPT_ZAD8 ${PYTHON_DIR}/plot_zad8.py)

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
    COMMAND python3 ${PLOT_SCRIPT_ZAD2} ${DATA_FILE_ZAD2} ${PLOT_FILE_ZAD2} "Siatka elementów z węzłami"
    COMMAND python3 ${PLOT_SCRIPT_HEATMAP} ${DATA_FILE_ZAD3} ${PLOT_FILE_ZAD3} "Mapa ciepła macierzy S przed modyfikacją"
    COMMAND python3 ${PLOT_SCRIPT_HEATMAP} ${DATA_FILE_ZAD4} ${PLOT_FILE_ZAD4} "Mapa ciepła macierzy S po modyfikacji"
    COMMAND python3 ${PLOT_SCRIPT_ZAD8} ${DATA_FILE_ZAD8_3} ${PLOT_FILE_ZAD8_3}
    COMMAND python3 ${PLOT_SCRIPT_ZAD8} ${DATA_FILE_ZAD8_10} ${PLOT_FILE_ZAD8_10}
)