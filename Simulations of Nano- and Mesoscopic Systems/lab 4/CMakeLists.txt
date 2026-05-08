cmake_minimum_required(VERSION 3.10)
project(nano_mezo_4 LANGUAGES CXX)
set(PROJECT_AUTHOR "Przemysław Ryś")

# --- Ustawienia kompilacji ---
set(CMAKE_BUILD_TYPE Release) # Debug
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

# --- DEBUG ---
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# --- RELEASE ---
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

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
set(PLOT_SCRIPT_ZAD1_2  ${PYTHON_DIR}/plot_zad1_2.py)
set(PLOT_SCRIPT_ZAD3    ${PYTHON_DIR}/plot_zad3.py)
set(PLOT_SCRIPT_ZAD4    ${PYTHON_DIR}/plot_zad4.py)

# --- Pliki danych ---
set(DATA_ZAD1   ${DATA_DIR}/zad1_single_barrier.dat)
set(DATA_ZAD2   ${DATA_DIR}/zad2_single_barrier_var_mass.dat)
set(DATA_ZAD3_TE ${DATA_DIR}/zad3_rtd_TE.dat)
set(DATA_ZAD3_IV ${DATA_DIR}/zad3_rtd_IV.dat)
set(DATA_ZAD4_EN  ${DATA_DIR}/zad4_En.dat)
set(DATA_ZAD4_GE  ${DATA_DIR}/zad4_GE.dat)
set(DATA_ZAD4_GVG ${DATA_DIR}/zad4_GVg.dat)

# --- Pliki wykresów ---
set(PLOT_ZAD1_DIR   ${PLOTS_DIR}/zad_1)
set(PLOT_ZAD2_DIR   ${PLOTS_DIR}/zad_2)
set(PLOT_ZAD3_TE ${PLOTS_DIR}/zad_3/TE.png)
set(PLOT_ZAD3_IV ${PLOTS_DIR}/zad_3/IV.png)
set(PLOT_ZAD4_EN  ${PLOTS_DIR}/zad_4/En.png)
set(PLOT_ZAD4_GE  ${PLOTS_DIR}/zad_4/GE.png)
set(PLOT_ZAD4_GVG ${PLOTS_DIR}/zad_4/GVg.png)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

## --- Eigen (header-only) ---
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen m OpenMP::OpenMP_CXX)
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

    COMMAND python3 "${PLOT_SCRIPT_ZAD1_2}" "${DATA_ZAD1}" "${DATA_ZAD2}" "${PLOT_ZAD1_DIR}" "${PLOT_ZAD2_DIR}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD3}"   "${DATA_ZAD3_TE}" "${DATA_ZAD3_IV}" "${PLOT_ZAD3_TE}" "${PLOT_ZAD3_IV}"
    COMMAND python3 "${PLOT_SCRIPT_ZAD4}"   "${DATA_ZAD4_EN}" "${DATA_ZAD4_GE}" "${DATA_ZAD4_GVG}" "${PLOT_ZAD4_EN}" "${PLOT_ZAD4_GE}" "${PLOT_ZAD4_GVG}"
    
)
