# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /gpfs1/home/gzpan2/cmake-3.26.4-linux-x86_64/bin/cmake

# The command to remove a file.
RM = /gpfs1/home/gzpan2/cmake-3.26.4-linux-x86_64/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gzpan2/scratch/virus/wxd/app/seqGraph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2

# Include any dependencies generated for this target.
include CMakeFiles/eResult.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/eResult.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/eResult.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/eResult.dir/flags.make

CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o: CMakeFiles/eResult.dir/flags.make
CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o: /home/gzpan2/scratch/virus/wxd/app/seqGraph/utils/e_result_reads.cpp
CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o: CMakeFiles/eResult.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gzpan2/scratch/virus/wxd/app/seqGraph/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o -MF CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o.d -o CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o -c /home/gzpan2/scratch/virus/wxd/app/seqGraph/utils/e_result_reads.cpp

CMakeFiles/eResult.dir/utils/e_result_reads.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/eResult.dir/utils/e_result_reads.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gzpan2/scratch/virus/wxd/app/seqGraph/utils/e_result_reads.cpp > CMakeFiles/eResult.dir/utils/e_result_reads.cpp.i

CMakeFiles/eResult.dir/utils/e_result_reads.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/eResult.dir/utils/e_result_reads.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gzpan2/scratch/virus/wxd/app/seqGraph/utils/e_result_reads.cpp -o CMakeFiles/eResult.dir/utils/e_result_reads.cpp.s

# Object files for target eResult
eResult_OBJECTS = \
"CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o"

# External object files for target eResult
eResult_EXTERNAL_OBJECTS =

eResult: CMakeFiles/eResult.dir/utils/e_result_reads.cpp.o
eResult: CMakeFiles/eResult.dir/build.make
eResult: CMakeFiles/eResult.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gzpan2/scratch/virus/wxd/app/seqGraph/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable eResult"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/eResult.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/eResult.dir/build: eResult
.PHONY : CMakeFiles/eResult.dir/build

CMakeFiles/eResult.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/eResult.dir/cmake_clean.cmake
.PHONY : CMakeFiles/eResult.dir/clean

CMakeFiles/eResult.dir/depend:
	cd /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gzpan2/scratch/virus/wxd/app/seqGraph /home/gzpan2/scratch/virus/wxd/app/seqGraph /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2 /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2 /home/gzpan2/scratch/virus/wxd/app/seqGraph/build2/CMakeFiles/eResult.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/eResult.dir/depend

