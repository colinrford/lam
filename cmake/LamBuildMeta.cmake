# LamBuildMeta.cmake
# Shared build metadata: git hash, debug flag, and boolean mapping helpers.

if(_LAM_BUILD_META_INCLUDED)
  return()
endif()
set(_LAM_BUILD_META_INCLUDED TRUE)

# Git hash
find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
else()
  set(GIT_HASH "unknown")
endif()

# Debug flag
if(CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
  set(IS_DEBUG "true")
else()
  set(IS_DEBUG "false")
endif()

# Boolean-to-string helper for config file generation
macro(map_bool_to_string VAR RESULT)
  if(${VAR})
    set(${RESULT} "true")
  else()
    set(${RESULT} "false")
  endif()
endmacro()
