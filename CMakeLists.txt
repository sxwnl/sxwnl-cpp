cmake_minimum_required(VERSION 3.6.0)
project(MultipleFiles)
include_directories(
	${PROJECT_SOURCE_DIR}/eph
	${PROJECT_SOURCE_DIR}/lunar
	${PROJECT_SOURCE_DIR}/
)
file(GLOB SRCS
	"${PROJECT_SOURCE_DIR}/eph/*"
	"${PROJECT_SOURCE_DIR}/lunar/*"
	"${PROJECT_SOURCE_DIR}/*"
)
#add_executable(MultipleFiles ${SRCS})
add_executable(test0 test/test.cpp ${SRCS})
add_executable(test1 test/test1.cpp ${SRCS})

