cmake_minimum_required(VERSION 3.6.0)
project(sxwnl-cpp)
set(CMAKE_VERBOSE_MAKEFILE on)
SET(CMAKE_BUILD_TYPE "Release")
# 编译参数
set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -std=c++1z -fexceptions")
set(CMAKE_CXX_FLAGS $ENV{CXXFLAGS})
message("CMAKE_CXX_FLAGS:" ${CMAKE_CXX_FLAGS})

include_directories(
	${PROJECT_SOURCE_DIR}/eph
	${PROJECT_SOURCE_DIR}/lunar
	${PROJECT_SOURCE_DIR}/mylib
)
file(GLOB SRCS
	"${PROJECT_SOURCE_DIR}/eph/*"
	"${PROJECT_SOURCE_DIR}/lunar/*"
	"${PROJECT_SOURCE_DIR}/mylib/*"
)
# 测试
# add_executable(test0 test/test.cpp ${SRCS})
add_executable(test1 test/test1.cpp ${SRCS})