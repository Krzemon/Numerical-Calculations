cmake_minimum_required(VERSION 3.10)
project(nano_mezo_3 LANGUAGES CXX)
set(PROJECT_AUTHOR "Przemysław Ryś")

# --- Ustawienia kompilacji ---
# -----------------------------
set(CMAKE_BUILD_TYPE Release) # Debug
# -----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

# --- DEBUG ---
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# --- RELEASE ---
# -----------------------------
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
# -----------------------------

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
set(PLOT_SCRIPT_ENERGIES   ${PYTHON_DIR}/plot_energies.py)
set(PLOT_SCRIPT_WAVEFUNCTIONS ${PYTHON_DIR}/plot_wavefunctions.py)
set(PLOT_SCRIPT_EVOLUTION     ${PYTHON_DIR}/plot_evolution.py)
set(PLOT_SCRIPT_EVOLUTION_MAP ${PYTHON_DIR}/plot_evolution_map.py)

# --- Pliki danych ---
set(DATA_ENERGIES_4STATES ${DATA_DIR}/energies_4states.dat)
set(DATA_WAVEFUNCTIONS ${DATA_DIR}/wavefunctions.dat)
set(DATA_EVOLUTION_RESONANCE ${DATA_DIR}/evolution_resonance.dat)
set(DATA_EVOLUTION_95PERCENT ${DATA_DIR}/evolution_95percent.dat)
set(DATA_EVOLUTION_MAP ${DATA_DIR}/map_p1_t_F.dat)

# --- Pliki wykresów ---
set(PLOT_ENERGIES_4STATES ${PLOTS_DIR}/zad_1/energies_4states.png)
set(PLOT_ENERGIES_2STATES ${PLOTS_DIR}/zad_1/energies_2states.png)
set(PLOT_WAVEFUNCTIONS ${PLOTS_DIR}/zad_2/wavefunctions.png)
set(PLOT_EVOLUTION_RESONANCE ${PLOTS_DIR}/zad_3/evolution_resonance.png)
set(PLOT_EVOLUTION_95PERCENT ${PLOTS_DIR}/zad_3/evolution_95percent.png)
set(PLOT_EVOLUTION_MAP ${PLOTS_DIR}/zad_4/map_p1_t_F.png)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

## --- Eigen (header-only) ---
find_package(Eigen3 REQUIRED)
# -----------------------------
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen m OpenMP::OpenMP_CXX)
# -----------------------------
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

    COMMAND python3 "${PLOT_SCRIPT_ENERGIES}" "${DATA_ENERGIES_4STATES}" "${PLOT_ENERGIES_4STATES}" "${PLOT_ENERGIES_2STATES}"
    COMMAND python3 "${PLOT_SCRIPT_WAVEFUNCTIONS}" "${DATA_WAVEFUNCTIONS}" "${PLOT_WAVEFUNCTIONS}"
    COMMAND python3 "${PLOT_SCRIPT_EVOLUTION}" "${DATA_EVOLUTION_RESONANCE}" "${DATA_EVOLUTION_95PERCENT}" "${PLOT_EVOLUTION_RESONANCE}" "${PLOT_EVOLUTION_95PERCENT}"
    COMMAND python3 "${PLOT_SCRIPT_EVOLUTION_MAP}" "${DATA_EVOLUTION_MAP}" "${PLOT_EVOLUTION_MAP}"
    
)
