# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /gpfs1/home/xuedowang2/app/cmake/bin/cmake

# The command to remove a file.
RM = /gpfs1/home/xuedowang2/app/cmake/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3

# Include any dependencies generated for this target.
include CMakeFiles/matching.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/matching.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/matching.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/matching.dir/flags.make

CMakeFiles/matching.dir/main.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/main.cpp.o: ../main.cpp
CMakeFiles/matching.dir/main.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/matching.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/main.cpp.o -MF CMakeFiles/matching.dir/main.cpp.o.d -o CMakeFiles/matching.dir/main.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/main.cpp

CMakeFiles/matching.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/main.cpp > CMakeFiles/matching.dir/main.cpp.i

CMakeFiles/matching.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/main.cpp -o CMakeFiles/matching.dir/main.cpp.s

CMakeFiles/matching.dir/src/Edge.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/Edge.cpp.o: ../src/Edge.cpp
CMakeFiles/matching.dir/src/Edge.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/matching.dir/src/Edge.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/Edge.cpp.o -MF CMakeFiles/matching.dir/src/Edge.cpp.o.d -o CMakeFiles/matching.dir/src/Edge.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Edge.cpp

CMakeFiles/matching.dir/src/Edge.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/Edge.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Edge.cpp > CMakeFiles/matching.dir/src/Edge.cpp.i

CMakeFiles/matching.dir/src/Edge.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/Edge.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Edge.cpp -o CMakeFiles/matching.dir/src/Edge.cpp.s

CMakeFiles/matching.dir/src/EndPoint.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/EndPoint.cpp.o: ../src/EndPoint.cpp
CMakeFiles/matching.dir/src/EndPoint.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/matching.dir/src/EndPoint.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/EndPoint.cpp.o -MF CMakeFiles/matching.dir/src/EndPoint.cpp.o.d -o CMakeFiles/matching.dir/src/EndPoint.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/EndPoint.cpp

CMakeFiles/matching.dir/src/EndPoint.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/EndPoint.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/EndPoint.cpp > CMakeFiles/matching.dir/src/EndPoint.cpp.i

CMakeFiles/matching.dir/src/EndPoint.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/EndPoint.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/EndPoint.cpp -o CMakeFiles/matching.dir/src/EndPoint.cpp.s

CMakeFiles/matching.dir/src/Exception.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/Exception.cpp.o: ../src/Exception.cpp
CMakeFiles/matching.dir/src/Exception.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/matching.dir/src/Exception.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/Exception.cpp.o -MF CMakeFiles/matching.dir/src/Exception.cpp.o.d -o CMakeFiles/matching.dir/src/Exception.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Exception.cpp

CMakeFiles/matching.dir/src/Exception.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/Exception.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Exception.cpp > CMakeFiles/matching.dir/src/Exception.cpp.i

CMakeFiles/matching.dir/src/Exception.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/Exception.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Exception.cpp -o CMakeFiles/matching.dir/src/Exception.cpp.s

CMakeFiles/matching.dir/src/Graph.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/Graph.cpp.o: ../src/Graph.cpp
CMakeFiles/matching.dir/src/Graph.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/matching.dir/src/Graph.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/Graph.cpp.o -MF CMakeFiles/matching.dir/src/Graph.cpp.o.d -o CMakeFiles/matching.dir/src/Graph.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Graph.cpp

CMakeFiles/matching.dir/src/Graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/Graph.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Graph.cpp > CMakeFiles/matching.dir/src/Graph.cpp.i

CMakeFiles/matching.dir/src/Graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/Graph.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Graph.cpp -o CMakeFiles/matching.dir/src/Graph.cpp.s

CMakeFiles/matching.dir/src/Junction.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/Junction.cpp.o: ../src/Junction.cpp
CMakeFiles/matching.dir/src/Junction.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/matching.dir/src/Junction.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/Junction.cpp.o -MF CMakeFiles/matching.dir/src/Junction.cpp.o.d -o CMakeFiles/matching.dir/src/Junction.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Junction.cpp

CMakeFiles/matching.dir/src/Junction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/Junction.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Junction.cpp > CMakeFiles/matching.dir/src/Junction.cpp.i

CMakeFiles/matching.dir/src/Junction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/Junction.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Junction.cpp -o CMakeFiles/matching.dir/src/Junction.cpp.s

CMakeFiles/matching.dir/src/Vertex.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/Vertex.cpp.o: ../src/Vertex.cpp
CMakeFiles/matching.dir/src/Vertex.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/matching.dir/src/Vertex.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/Vertex.cpp.o -MF CMakeFiles/matching.dir/src/Vertex.cpp.o.d -o CMakeFiles/matching.dir/src/Vertex.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Vertex.cpp

CMakeFiles/matching.dir/src/Vertex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/Vertex.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Vertex.cpp > CMakeFiles/matching.dir/src/Vertex.cpp.i

CMakeFiles/matching.dir/src/Vertex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/Vertex.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Vertex.cpp -o CMakeFiles/matching.dir/src/Vertex.cpp.s

CMakeFiles/matching.dir/src/Weight.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/Weight.cpp.o: ../src/Weight.cpp
CMakeFiles/matching.dir/src/Weight.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/matching.dir/src/Weight.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/Weight.cpp.o -MF CMakeFiles/matching.dir/src/Weight.cpp.o.d -o CMakeFiles/matching.dir/src/Weight.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Weight.cpp

CMakeFiles/matching.dir/src/Weight.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/Weight.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Weight.cpp > CMakeFiles/matching.dir/src/Weight.cpp.i

CMakeFiles/matching.dir/src/Weight.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/Weight.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/Weight.cpp -o CMakeFiles/matching.dir/src/Weight.cpp.s

CMakeFiles/matching.dir/src/matching.cpp.o: CMakeFiles/matching.dir/flags.make
CMakeFiles/matching.dir/src/matching.cpp.o: ../src/matching.cpp
CMakeFiles/matching.dir/src/matching.cpp.o: CMakeFiles/matching.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/matching.dir/src/matching.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/matching.dir/src/matching.cpp.o -MF CMakeFiles/matching.dir/src/matching.cpp.o.d -o CMakeFiles/matching.dir/src/matching.cpp.o -c /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/matching.cpp

CMakeFiles/matching.dir/src/matching.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matching.dir/src/matching.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/matching.cpp > CMakeFiles/matching.dir/src/matching.cpp.i

CMakeFiles/matching.dir/src/matching.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matching.dir/src/matching.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/src/matching.cpp -o CMakeFiles/matching.dir/src/matching.cpp.s

# Object files for target matching
matching_OBJECTS = \
"CMakeFiles/matching.dir/main.cpp.o" \
"CMakeFiles/matching.dir/src/Edge.cpp.o" \
"CMakeFiles/matching.dir/src/EndPoint.cpp.o" \
"CMakeFiles/matching.dir/src/Exception.cpp.o" \
"CMakeFiles/matching.dir/src/Graph.cpp.o" \
"CMakeFiles/matching.dir/src/Junction.cpp.o" \
"CMakeFiles/matching.dir/src/Vertex.cpp.o" \
"CMakeFiles/matching.dir/src/Weight.cpp.o" \
"CMakeFiles/matching.dir/src/matching.cpp.o"

# External object files for target matching
matching_EXTERNAL_OBJECTS =

matching: CMakeFiles/matching.dir/main.cpp.o
matching: CMakeFiles/matching.dir/src/Edge.cpp.o
matching: CMakeFiles/matching.dir/src/EndPoint.cpp.o
matching: CMakeFiles/matching.dir/src/Exception.cpp.o
matching: CMakeFiles/matching.dir/src/Graph.cpp.o
matching: CMakeFiles/matching.dir/src/Junction.cpp.o
matching: CMakeFiles/matching.dir/src/Vertex.cpp.o
matching: CMakeFiles/matching.dir/src/Weight.cpp.o
matching: CMakeFiles/matching.dir/src/matching.cpp.o
matching: CMakeFiles/matching.dir/build.make
matching: CMakeFiles/matching.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX executable matching"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matching.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/matching.dir/build: matching
.PHONY : CMakeFiles/matching.dir/build

CMakeFiles/matching.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/matching.dir/cmake_clean.cmake
.PHONY : CMakeFiles/matching.dir/clean

CMakeFiles/matching.dir/depend:
	cd /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3 /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3 /gpfs1/scratch/ResearchGroups/cs_shuaicli/wxd/app/seqGraph/build3/CMakeFiles/matching.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/matching.dir/depend

