# CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

# PROJECT(MkLib Tester)
# SET(proj_name "CMake MkLib Tester")
# SET(CMAKE_CXX_STANDARD 20)

# ADD_SUBDIRECTORY(src)
# ADD_SUBDIRECTORY(tests)
# # ADD_SUBDIRECTORY(extern/googletest)

# SET (srcs test.cpp )

# ADD_EXECUTABLE(test ${srcs})
# TARGET_LINK_LIBRARIES(test ../build/MkLib glu32 freeglut)
SET(STATIC_LIB ${proj_name}_static)
SET(DYNAMIC_LIB ${proj_name})

SET(TESTER ${CMAKE_PROJECT_NAME}_test)
# SET(TESTER test)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR} )

# file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)
# SET(SOURCES ${TEST_SOURCES})
SET(TEST_SOURCES test.cpp)

add_executable(${TESTER} ${TEST_SOURCES})
add_test(NAME ${TESTER} COMMAND ${TESTER})
target_link_libraries(${TESTER} PUBLIC ${STATIC_LIB} )


