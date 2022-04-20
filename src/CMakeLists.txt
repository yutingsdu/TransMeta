cmake_minimum_required(VERSION 2.6)
project (gingko_graph)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#Warning: The followings need your settings...(we have shown some examples below them.)
#set(BOOST_LIB_DIR <Path_To_BOOST_LIB_DIR>)
set(BOOST_LIB_DIR /rhome/jliu/bigdata/local/boost/lib)
#set(BOOST_INCLUDE_DIR <Path_To_BOOST_INCLUDE_DIR>)
set(BOOST_INCLUDE_DIR /rhome/jliu/bigdata/local/boost/include)
#set(BAMTOOLS_LIB_DIR <Path_To_BAMTOOLS_LIB_DIR>)
set(BAMTOOLS_LIB_DIR /rhome/jliu/bigdata/yuting/Software/bamtools/lib64)
#set(BAMTOOLS_INCLUDE_DIR <Path_To_BAMTOOLS_INCLUDE_DIR>)
set(BAMTOOLS_INCLUDE_DIR /rhome/jliu/bigdata/yuting/Software/bamtools/include/bamtools/)

include_directories(. ${BOOST_INCLUDE_DIR} ${BAMTOOLS_INCLUDE_DIR})
link_directories(${BOOST_LIB_DIR} ${BAMTOOLS_LIB_DIR})

add_executable(gingko_path_search gingko_path_search.cpp ./unitepathsearch/junction_graph.cc
		./unitepathsearch/QuadProg++.cc  ./unitepathsearch/simplify_graph.cc  )
add_executable(${PROJECT_NAME} gingko_graph.cc assembler.cpp simplify_graph.cc junction_graph.cc QuadProg++.cc)
add_executable(gingko_individual final_individual.cc)
add_executable(gingko_merge final_merge.cc)
add_executable(gingko_cover get_covered_reference.cc)

target_link_libraries(${PROJECT_NAME} boost_math_c99 bamtools boost_program_options boost_regex boost_graph boost_thread boost_system boost_filesystem z)
