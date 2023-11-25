CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

SET(proj_name "MkBEM")
PROJECT(${proj_name})
SET( STATIC_LIB ${CMAKE_PROJECT_NAME}_static )
SET( DYNAMIC_LIB ${CMAKE_PROJECT_NAME} )

SET(CMAKE_CXX_STANDARD 20)

# -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg
SET(CMAKE_CXX_FLAGS -pg)
SET(CMAKE_EXE_LINKER_FLAGS  -pg)

# SET (INC_FILES_FOR_EXE
#         src/mainunit.h
#         src/bemdata.h
#         src/bemunit.h
#         src/inputunit.h
# )

# SET (SRC_FILES_FOR_EXE
#         src/mainunit.cpp
#         src/bemdata.cpp
#         src/bemunit.cpp
#         src/inputunit.cpp
# )

SET (INC_FILES_FOR_LIB
        src/bemdata.h
        src/inputunit.h
        src/bemunit.h
)

SET (SRC_FILES_FOR_LIB
        src/bemdata.cpp
        src/inputunit.cpp
        src/bemunit.cpp
)
# for mac
INCLUDE_DIRECTORIES(src) 
# INCLUDE_DIRECTORIES(/Users/myungkyusong/Documents/Develop/01_THMBEMSPH/01_MkLIB/src)
# INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../MkInclude) # obsolete
INCLUDE_DIRECTORIES(/opt/homebrew/include)
INCLUDE_DIRECTORIES(/usr/local/include)
# for windows
# INCLUDE_DIRECTORIES(src)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

LINK_DIRECTORIES(/opt/homebrew/lib)
LINK_DIRECTORIES(/usr/local/lib)
# LINK_DIRECTORIES(/Users/myungkyusong/Documents/Develop/01_THMBEMSPH/01_MkLIB)
# LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/../MkInclude) # obsolete

set_target_properties(${DYNAMIC_LIB} PROPERTIES PUBLIC_HEADER "${INC_FILES}")
# set_target_properties(${STATIC_LIB} PROPERTIES PUBLIC_HEADER "${INC_FILES}")

INSTALL(TARGETS ${DYNAMIC_LIB}
        # LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/../MkLIB
        # PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/../MkInclude
        LIBRARY DESTINATION /usr/local/lib
        PUBLIC_HEADER DESTINATION /usr/local/include
)

INSTALL(TARGETS ${STATIC_LIB}
        # LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/../MkLIB
        # PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/../MkInclude
        LIBRARY DESTINATION /usr/local/lib
        PUBLIC_HEADER DESTINATION /usr/local/include

)

