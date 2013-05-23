CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(colhck-tests)

INCLUDE_DIRECTORIES(
  ${kazmath_SOURCE_DIR}/src
  ${glhck_SOURCE_DIR}/include
  ${GLFW_SOURCE_DIR}/include
)

ADD_EXECUTABLE(coltest coltest.c)
TARGET_LINK_LIBRARIES(coltest glhck glfw ${GLFW_LIBRARIES})
