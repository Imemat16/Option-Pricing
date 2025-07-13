# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Debug")
  file(REMOVE_RECURSE
  "BinaryOptionsPricing_autogen"
  "CMakeFiles\\BinaryOptionsPricing_autogen.dir\\AutogenUsed.txt"
  "CMakeFiles\\BinaryOptionsPricing_autogen.dir\\ParseCache.txt"
  )
endif()
