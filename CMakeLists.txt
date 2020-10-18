cmake_minimum_required(VERSION 3.0)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    message("ccache used.")
endif() 

project(palabos)

enable_language(CXX)

# To link palabos in directories other than this
cmake_policy(SET CMP0079 NEW)

add_subdirectory(palaboslib)
add_subdirectory(cylinder2d)