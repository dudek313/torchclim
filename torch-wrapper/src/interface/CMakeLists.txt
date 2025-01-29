
enable_language(C)

#needed for libtorch ans stc
enable_language(CXX)

enable_language(Fortran)


# verify the compatibility of the C/Fortran and C++/Fortran compilers
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

include_directories(${ONNXRUNTIME_ROOTDIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_library(torch-plugin SHARED torch-plugin.f90)
#add_library(torch-plugin STATIC torch-plugin.f90)
set_property(TARGET torch-plugin PROPERTY CXX_STANDARD 17)

target_include_directories(torch-plugin PUBLIC
    "${ONNXRUNTIME_ROOTDIR}/include"                           # Pre-built package
    "${ONNXRUNTIME_ROOTDIR}/include/onnxruntime"               # Linux local install to /usr/local
    "${ONNXRUNTIME_ROOTDIR}/include/onnxruntime/core/session"  # Windows local install
)

target_link_directories(torch-plugin PUBLIC "${ONNXRUNTIME_ROOTDIR}/lib")

set(ONNXRUNTIME_LIBRARIES "libonnxruntime.so")

target_link_libraries(torch-plugin PUBLIC
    ${ONNXRUNTIME_LIBRARIES}
)

target_sources(torch-plugin
    PUBLIC
      torch-plugin.f90
    PRIVATE
    #torch-wrap-cdef.f90
      torch-wrap.cpp

  )
