CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

SET(proj_name "MkBEM")
PROJECT(${proj_name})

SET(CMAKE_CXX_STANDARD 20)

# -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg
SET(CMAKE_CXX_FLAGS -pg)
SET(CMAKE_EXE_LINKER_FLAGS -pg)

SET (INC_FILES
        src/MkAbstract.hpp
        src/MkArray.hpp 
        src/MkColor.hpp
        src/MkCircle.hpp
        src/MkCylinder.hpp
        src/MkCube.hpp
        src/MkLine.hpp 
        src/MkMatrix.hpp
        src/MkMisc.hpp 
        src/MkPlane.hpp
        src/MkPoint.hpp 
        src/MkPointData.hpp 
        src/MkPolygon.hpp
        src/MkRect.hpp
        src/MkShape.hpp
        src/MkSphere.hpp
        src/MkTriangle.hpp    

)

# for mac
INCLUDE_DIRECTORIES(src) 
# INCLUDE_DIRECTORIES(/Users/myungkyusong/Documents/Develop/01_THMBEMSPH/01_MkLIB/src)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../MkInclude)
INCLUDE_DIRECTORIES(/opt/homebrew/include)

# for windows
# INCLUDE_DIRECTORIES(src)

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

LINK_DIRECTORIES(/opt/homebrew/lib)
# LINK_DIRECTORIES(/Users/myungkyusong/Documents/Develop/01_THMBEMSPH/01_MkLIB)
LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/../MkInclude)

set_target_properties(MkLIB PROPERTIES PUBLIC_HEADER "${INC_FILES}")
INSTALL(TARGETS MkLIB
        LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/../MkLib
        PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/../MkInclude
)

