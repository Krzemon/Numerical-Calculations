cmake_minimum_required(VERSION 3.10)
project(lab_05)
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

set(PLOT_FILE_ZAD1 ${PLOTS_DIR}/zad_1_grid.png)

# --- Skrypty wykresow ---
set(PLOT_SCRIPT_ZAD1 ${PYTHON_DIR}/plot_zad1_grid.py)

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
    COMMAND python3 ${PLOT_SCRIPT_ZAD1} ${DATA_FILE_ZAD1} ${PLOT_FILE_ZAD1} "Siatka elementów z węzłami"
)