cmake_minimum_required(VERSION 3.10)
project(nano_mezo_2 LANGUAGES CXX)
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
set(PLOT_SCRIPT_CONVERGENCE ${PYTHON_DIR}/plot_convergence.py)
set(PLOT_SCRIPT_E_VS_A ${PYTHON_DIR}/plot_E_vs_a.py)
set(PLOT_SCRIPT_DENSITY ${PYTHON_DIR}/plot_density.py)

# --- Pliki danych ---
set(DATA_CONVERGENCE ${DATA_DIR}/energies_k.dat)
set(DATA_E_VS_A_IT ${DATA_DIR}/E_vs_a_IT.dat)
set(DATA_E_VS_A_HF ${DATA_DIR}/E_vs_a_HF.dat)
set(DATA_DENSITY_A30 ${DATA_DIR}/density_a30.dat)
set(DATA_DENSITY_A60 ${DATA_DIR}/density_a60.dat)

# --- Pliki wykresów ---
set(PLOT_CONVERGENCE ${PLOTS_DIR}/zad_1/energies_k.png)
set(PLOT_E_VS_A ${PLOTS_DIR}/zad_2_4/E_vs_a.png)
set(PLOT_DENSITY_A30 ${PLOTS_DIR}/zad_3/density_a30.png)
set(PLOT_DENSITY_A60 ${PLOTS_DIR}/zad_3/density_a60.png)

# --- Zrodla ---
file(GLOB SOURCES "${SRC_DIR}/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})

## --- Eigen (header-only) ---
find_package(Eigen3 REQUIRED)
# target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen m)
target_link_libraries(${PROJECT_NAME} PRIVATE m)
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

    COMMAND python3 "${PLOT_SCRIPT_CONVERGENCE}" "${DATA_CONVERGENCE}" "${PLOT_CONVERGENCE}"
    COMMAND python3 "${PLOT_SCRIPT_E_VS_A}" "${DATA_E_VS_A_IT}" "${DATA_E_VS_A_HF}" "${PLOT_E_VS_A}"
    COMMAND python3 "${PLOT_SCRIPT_DENSITY}" "${DATA_DENSITY_A30}" "${PLOT_DENSITY_A30}" --mode raw # --mode raw
    COMMAND python3 "${PLOT_SCRIPT_DENSITY}" "${DATA_DENSITY_A60}" "${PLOT_DENSITY_A60}" --mode raw # --mode interp --n 200
)