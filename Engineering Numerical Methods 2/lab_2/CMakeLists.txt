cmake_minimum_required(VERSION 3.10)
project(lab_02)
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

# --- Pliki wynikowe ---
set(DATA_FILE_2 ${DATA_DIR}/zad_2.dat)
set(DATA_FILE_3 ${DATA_DIR}/zad_3.dat)
set(DATA_FILE_4 ${DATA_DIR}/zad_4.dat)
set(DATA_FILE_5_1 ${DATA_DIR}/zad_5_1.dat)
set(DATA_FILE_5_2 ${DATA_DIR}/zad_5_2.dat)
set(DATA_FILE_5_3 ${DATA_DIR}/zad_5_3.dat)

set(PLOT_FILE_U ${PLOTS_DIR}/rozwiazanie_u.png)
set(PLOT_FILE_2 ${PLOTS_DIR}/zad_2.png)
set(PLOT_FILE_3 ${PLOTS_DIR}/zad_3.png)
set(PLOT_FILE_4 ${PLOTS_DIR}/zad_4.png)
set(PLOT_FILE_5_1 ${PLOTS_DIR}/zad_5_1.png)
set(PLOT_FILE_5_2 ${PLOTS_DIR}/zad_5_2.png)
set(PLOT_FILE_5_3 ${PLOTS_DIR}/zad_5_3.png)

# --- Skrypty wykresow ---
set(PLOT_SCRIPT_U ${PROJECT_SOURCE_DIR}/python/plot_u.py)
set(PLOT_SCRIPT_DELTA_U ${PROJECT_SOURCE_DIR}/python/plot_delta_u.py)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Linkowanie ---
target_link_libraries(${PROJECT_NAME} PRIVATE gsl gslcblas m)
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# --- Po kompilacji: uruchomienie programu i generowanie wykresu ---
add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}"
    COMMAND python3 "${PLOT_SCRIPT_U}" "${DATA_FILE_2}" "${PLOT_FILE_U}"
    COMMAND python3 "${PLOT_SCRIPT_DELTA_U}" "${DATA_FILE_2}" "${PLOT_FILE_2}" "Metoda kolokacji - u_d - u_n dla bazy 1" 
    COMMAND python3 "${PLOT_SCRIPT_DELTA_U}" "${DATA_FILE_3}" "${PLOT_FILE_3}" "Metoda najmniejszych kwadratów - u_d - u_n dla bazy 1" 
    COMMAND python3 "${PLOT_SCRIPT_DELTA_U}" "${DATA_FILE_4}" "${PLOT_FILE_4}" "Metoda Galerkina - u_d - u_n dla bazy 1" 
    COMMAND python3 "${PLOT_SCRIPT_DELTA_U}" "${DATA_FILE_5_1}" "${PLOT_FILE_5_1}" "Metoda kolokacji - u_d - u_n dla bazy 2" 
    COMMAND python3 "${PLOT_SCRIPT_DELTA_U}" "${DATA_FILE_5_2}" "${PLOT_FILE_5_2}" "Metoda najmniejszych kwadratów - u_d - u_n dla bazy 2" 
    COMMAND python3 "${PLOT_SCRIPT_DELTA_U}" "${DATA_FILE_5_3}" "${PLOT_FILE_5_3}" "Metoda Galerkina - u_d - u_n dla bazy 2" 
)