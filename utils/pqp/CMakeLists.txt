cmake_minimum_required(VERSION 3.20)

# Distance_OBB.cpp used Eigen
find_package(Eigen3 CONFIG REQUIRED)

add_library(pqp STATIC
    BV.cpp
    Build.cpp
    PQP.cpp
    TriDist.cpp
    Distance_OBB.cpp
)

target_include_directories(pqp PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

# Link Eigen3（header-only）
target_link_libraries(pqp PUBLIC Eigen3::Eigen)