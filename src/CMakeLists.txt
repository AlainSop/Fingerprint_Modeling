# Create variable containing all the CXX sources
set(SOURCES
        ${PROJECT_SOURCE_DIR}/src/Image.cpp ${PROJECT_SOURCE_DIR}/src/Restoration.cpp ${PROJECT_SOURCE_DIR}/src/GeometricalWarps.cpp
        ${PROJECT_SOURCE_DIR}/src/LinearFiltering.cpp ${PROJECT_SOURCE_DIR}/src/MorphologicalFiltering.cpp ${PROJECT_SOURCE_DIR}/src/RegistrationOptimization.cpp
        ${PROJECT_SOURCE_DIR}/src/Simulation.cpp
        )

# Add all the sources to the project
add_library(${PROJECT_NAME} ${SOURCES})

