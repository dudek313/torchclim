# needed to find the backtrace library
#enable_language(C)
# needed for the random number generator
#enable_language(CXX)
#enable_language(Fortran)

# verify the compatibility of the C/Fortran and C++/Fortran compilers
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


add_subdirectory(interface)
#add_subdirectory(utils)

link_directories(interface)


add_executable(plugin-test plugin-test.f90)
set_property(TARGET plugin-test PROPERTY CXX_STANDARD 14)
#include_directories(interface)

target_include_directories(plugin-test
	PUBLIC
	  interface) 

#add_library(torch-plugin STATIC SHARED torch-plugin.f90)



#find_package(Torch REQUIRED)


#target_link_libraries(torch-plugin
#  PRIVATE
#  "${TORCH_LIBRARIES}"
#  )

target_sources(plugin-test
    PRIVATE
      cam4_profile.f90
  )


target_link_libraries(plugin-test
  PRIVATE
    torch-plugin
  )


