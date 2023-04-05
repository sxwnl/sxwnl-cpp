cmake_minimum_required(VERSION 3.6.0)
project(sxwnl-cpp)
set(CMAKE_VERBOSE_MAKEFILE on)
set(PRIVATE_CXX_FLAGS
	"-std=c++1z"
	"-fexceptions"
	"-g"
	"-Os"
	)

# set compile flags
add_compile_options(${PRIVATE_CXX_FLAGS})

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