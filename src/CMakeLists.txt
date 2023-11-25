CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

# PROJECT(MkBEM)
# SET(proj_name "MkBEM")
# SET(CMAKE_CXX_STANDARD 20)

SET(STATIC_LIB ${proj_name}_static)
SET(DYNAMIC_LIB ${proj_name})

SET( EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )
SET( LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build )
# SET( LIBRARY_OUTPUT_PATH /usr/local/lib )
SET( BINARY ${CMAKE_PROJECT_NAME} )

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.h *.cpp)
file(GLOB GLLIB "/opt/homebrew/lib/libGL.dylib")
file(GLOB GLUTLIB "/opt/homebrew/lib/libglut.dylib")
file(GLOB MkLIB "/usr/local/lib/libMkLIB.dylib")

# set sources
SET (INC_FILES_FOR_LIB
        bemdata.h
        inputunit.h
        bemunit.h
)

SET (SRC_FILES_FOR_LIB
        bemdata.cpp
        inputunit.cpp
        bemunit.cpp
)

# SET(SRC_FILES ${SCOURCES}) 

ADD_DEFINITIONS(-DMKDEBUG)

#build
# ADD_EXECUTABLE(${BINARY}_RUN ${SRC_FILES})
ADD_LIBRARY (${DYNAMIC_LIB} SHARED ${SRC_FILES_FOR_LIB})
ADD_LIBRARY (${STATIC_LIB} ${SRC_FILES_FOR_LIB})

## Windows
# TARGET_LINK_LIBRARIES(${BINARY} glu32 freeglut)

## macOS
TARGET_LINK_LIBRARIES(${DYNAMIC_LIB} ${GLLIB} ${GLUTLIB} )
TARGET_LINK_LIBRARIES(${STATIC_LIB} ${GLLIB} ${GLUTLIB})

## Linux
# TARGET_LINK_LIBRARIES(${BINARY} GL glut)

