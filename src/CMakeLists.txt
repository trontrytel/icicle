## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
#  @section DESCRIPTION
#    a sanity check to ensure all header files are self-contained
#    and may be safely included multiple times

execute_process(
  COMMAND "ls" "-1" "src"
  COMMAND "fgrep" ".hpp" 
  COMMAND "awk" "{printf(\"%s;\",substr($1,1,index($1,\".hpp\")-1))}" 
  OUTPUT_VARIABLE HPPS
)
include(CheckCXXSourceCompiles)
set(TMP CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -I${CMAKE_CURRENT_SOURCE_DIR}")
foreach(HPP ${HPPS})
  check_cxx_source_compiles("#include <${HPP}.hpp>\n#include <${HPP}.hpp>\nint main() {}" ${HPP}_hpp)
  if (NOT ${HPP}_hpp)
    message(WARNING 
      "header file not self-contained or not protected against multiple-inclusion: 
      src/${HPP}.hpp (c.f. CMakeFiles/CMakeError.log)"
    )
  endif()
endforeach()
set(CMAKE_REQUIRED_FLAGS ${TMP})
