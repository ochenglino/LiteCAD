# Prohibit CGAL from using GMP/MPFR (this project only uses the floating-point kernel, which is incompatible with MSVC)
set(CGAL_DISABLE_GMP TRUE)

# Find dependencies installed via vcpkg
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

# ANN precompiled library (dynamic link, ANN.dll required at runtime)
set(ANN_LIB "${CMAKE_SOURCE_DIR}/utils/ann/lib/ANN.lib")
set(ANN_DLL "${CMAKE_SOURCE_DIR}/utils/ann/lib/ANN.dll")

if(EXISTS "${ANN_LIB}")
    add_library(ann SHARED IMPORTED)
    set_target_properties(ann PROPERTIES
        IMPORTED_IMPLIB "${ANN_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/utils/ann/include"
    )
    if(EXISTS "${ANN_DLL}")
        set_target_properties(ann PROPERTIES
            IMPORTED_LOCATION "${ANN_DLL}"
        )
    endif()
    set(ANN_TARGET ann)
else()
    set(ANN_TARGET "")
    message(WARNING "ANN.lib not found at utils/ann/lib/ANN.lib. ANN linking may fail.")
endif()

add_executable(LiteCAD
    main.cpp
    BaseModel.cpp
    Point3D.cpp
    Inputed_Triangle_Mesh_Reconstruction.cpp
    Reconstruction_Basedon_Paricles_LBFGS.cpp
)

target_include_directories(LiteCAD PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/utils/pqp
    ${CMAKE_SOURCE_DIR}/utils/ann/include
    ${CMAKE_SOURCE_DIR}/utils/lbfgs
)

target_link_libraries(LiteCAD PRIVATE
    pqp
    CGAL::CGAL
    Eigen3::Eigen
    ${ANN_TARGET}
    OpenMP::OpenMP_CXX 
)

#compile lbfgs sources
target_sources(LiteCAD PRIVATE
    ${CMAKE_SOURCE_DIR}/utils/lbfgs/lbfgs.cpp
)

# MSVC set
target_compile_definitions(LiteCAD PRIVATE
    CGAL_DO_NOT_USE_MPZF
    NOMINMAX
    _HAS_STD_BYTE=0
)

# copy ANN.dll to executable file directory
if(EXISTS "${ANN_DLL}")
    add_custom_command(TARGET LiteCAD POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${ANN_DLL}"
            "$<TARGET_FILE_DIR:LiteCAD>/ANN.dll"
    )
endif()