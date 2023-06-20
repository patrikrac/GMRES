cmake_minimum_required( VERSION 3.10 )

set (CMAKE_CXX_STANDARD 17)

project(proj LANGUAGES CXX)

SET( APP_EXE proj)

SET( TEST_EXE test)

ADD_EXECUTABLE ( ${APP_EXE} 
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    ${CMAKE_SOURCE_DIR}/include/vector_op.hpp
    ${CMAKE_SOURCE_DIR}/include/vector_space.hpp
    ${CMAKE_SOURCE_DIR}/include/dense_matrix.hpp
    ${CMAKE_SOURCE_DIR}/include/map_matrix.hpp
    ${CMAKE_SOURCE_DIR}/include/lu.hpp
    ${CMAKE_SOURCE_DIR}/src/solvers.cpp
    )


ADD_EXECUTABLE ( ${TEST_EXE} 
    ${CMAKE_SOURCE_DIR}/src/test.cpp
    ${CMAKE_SOURCE_DIR}/include/vector_op.hpp
    ${CMAKE_SOURCE_DIR}/include/vector_space.hpp
    ${CMAKE_SOURCE_DIR}/include/dense_matrix.hpp
    ${CMAKE_SOURCE_DIR}/include/map_matrix.hpp
    ${CMAKE_SOURCE_DIR}/include/lu.hpp
    ${CMAKE_SOURCE_DIR}/src/solvers.cpp
    )

target_include_directories(${APP_EXE} PRIVATE  
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ) 

target_include_directories(${TEST_EXE} PRIVATE  
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ) 

target_compile_options(${APP_EXE}
    PRIVATE -O3
)

target_compile_options(${TEST_EXE}
    PRIVATE -O3
)
