cmake_minimum_required(VERSION 3.10)
project(lab_01)
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
set(DATA_FILE1 ${DATA_DIR}/zad_1.dat)
set(DATA_FILE2 ${DATA_DIR}/zad_2.dat)
set(DATA_FILE30 ${DATA_DIR}/zad_3_0.dat)
set(DATA_FILE31 ${DATA_DIR}/zad_3_1.dat)
set(DATA_FILE40 ${DATA_DIR}/zad_4_0.dat)
set(DATA_FILE41 ${DATA_DIR}/zad_4_1.dat)
set(DATA_FILE42 ${DATA_DIR}/zad_4_2.dat)

set(PLOT_FILE1 ${PLOTS_DIR}/zad_1.png)
set(PLOT_FILE2 ${PLOTS_DIR}/zad_2.png)
set(PLOT_FILE30 ${PLOTS_DIR}/zad_3_0.png)
set(PLOT_FILE31 ${PLOTS_DIR}/zad_3_1.png)
set(PLOT_FILE40 ${PLOTS_DIR}/zad_4_0.png)
set(PLOT_FILE41 ${PLOTS_DIR}/zad_4_1.png)
set(PLOT_FILE42 ${PLOTS_DIR}/zad_4_2.png)

# --- Skrypty wykresów ---
set(PLOT_SCRIPT1 ${PROJECT_SOURCE_DIR}/python/plot_zad1.py)
set(PLOT_SCRIPT2 ${PROJECT_SOURCE_DIR}/python/plot_zad2.py)
set(PLOT_SCRIPT3 ${PROJECT_SOURCE_DIR}/python/plot_zad3.py)
set(PLOT_SCRIPT40 ${PROJECT_SOURCE_DIR}/python/plot_zad4_0.py)
set(PLOT_SCRIPT41 ${PROJECT_SOURCE_DIR}/python/plot_zad4_1.py)
set(PLOT_SCRIPT42 ${PROJECT_SOURCE_DIR}/python/plot_zad4_2.py)

# --- Źródła ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Linkowanie: GSL + CBLAS + libm ---
target_link_libraries(${PROJECT_NAME} PRIVATE gsl gslcblas m)
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})
target_compile_definitions(${PROJECT_NAME} PRIVATE PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# --- Tworzenie katalogow ---
file(MAKE_DIRECTORY ${DATA_DIR})
file(MAKE_DIRECTORY ${PLOTS_DIR})

# --- Po kompilacji: uruchomienie programu i generowanie wykresu ---
add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}"
    COMMAND python3 "${PLOT_SCRIPT1}" "${DATA_FILE1}" "${PLOT_FILE1}"
    COMMAND python3 "${PLOT_SCRIPT2}" "${DATA_FILE2}" "${PLOT_FILE2}"
    COMMAND python3 "${PLOT_SCRIPT3}" "${DATA_FILE30}" "${PLOT_FILE30}"
    COMMAND python3 "${PLOT_SCRIPT3}" "${DATA_FILE31}" "${PLOT_FILE31}"
    COMMAND python3 "${PLOT_SCRIPT40}" "${DATA_FILE40}" "${PLOT_FILE40}"
    COMMAND python3 "${PLOT_SCRIPT42}" "${DATA_FILE42}" "${PLOT_FILE42}"
)
