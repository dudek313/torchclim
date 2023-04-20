
enable_language(C)

#needed for libtorch ans stc
enable_language(CXX)

enable_language(Fortran)


# verify the compatibility of the C/Fortran and C++/Fortran compilers
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

include_directories(
	    "${TORCH_INCLUDE_DIRS}" # libtorch headers
	    )



add_library(torch-plugin SHARED torch-plugin.f90)
#add_library(torch-plugin STATIC torch-plugin.f90)
set_property(TARGET torch-plugin PROPERTY CXX_STANDARD 14)

target_link_libraries(torch-plugin
    PRIVATE
    "${TORCH_LIBRARIES}"
  )

target_sources(torch-plugin
    PUBLIC
      torch-plugin.f90
    PRIVATE
    #torch-wrap-cdef.f90
      torch-wrap.cpp

  )
