# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.26

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2023.2\bin\cmake\win\x64\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2023.2\bin\cmake\win\x64\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\halin\Documents\GitHub\multigrid

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/multigrid.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/multigrid.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/multigrid.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/multigrid.dir/flags.make

CMakeFiles/multigrid.dir/main.cpp.obj: CMakeFiles/multigrid.dir/flags.make
CMakeFiles/multigrid.dir/main.cpp.obj: C:/Users/halin/Documents/GitHub/multigrid/main.cpp
CMakeFiles/multigrid.dir/main.cpp.obj: CMakeFiles/multigrid.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/multigrid.dir/main.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/multigrid.dir/main.cpp.obj -MF CMakeFiles\multigrid.dir\main.cpp.obj.d -o CMakeFiles\multigrid.dir\main.cpp.obj -c C:\Users\halin\Documents\GitHub\multigrid\main.cpp

CMakeFiles/multigrid.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multigrid.dir/main.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\halin\Documents\GitHub\multigrid\main.cpp > CMakeFiles\multigrid.dir\main.cpp.i

CMakeFiles/multigrid.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multigrid.dir/main.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\halin\Documents\GitHub\multigrid\main.cpp -o CMakeFiles\multigrid.dir\main.cpp.s

CMakeFiles/multigrid.dir/multigrid3d.cpp.obj: CMakeFiles/multigrid.dir/flags.make
CMakeFiles/multigrid.dir/multigrid3d.cpp.obj: C:/Users/halin/Documents/GitHub/multigrid/multigrid3d.cpp
CMakeFiles/multigrid.dir/multigrid3d.cpp.obj: CMakeFiles/multigrid.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/multigrid.dir/multigrid3d.cpp.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/multigrid.dir/multigrid3d.cpp.obj -MF CMakeFiles\multigrid.dir\multigrid3d.cpp.obj.d -o CMakeFiles\multigrid.dir\multigrid3d.cpp.obj -c C:\Users\halin\Documents\GitHub\multigrid\multigrid3d.cpp

CMakeFiles/multigrid.dir/multigrid3d.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/multigrid.dir/multigrid3d.cpp.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\halin\Documents\GitHub\multigrid\multigrid3d.cpp > CMakeFiles\multigrid.dir\multigrid3d.cpp.i

CMakeFiles/multigrid.dir/multigrid3d.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/multigrid.dir/multigrid3d.cpp.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.2\bin\mingw\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\halin\Documents\GitHub\multigrid\multigrid3d.cpp -o CMakeFiles\multigrid.dir\multigrid3d.cpp.s

# Object files for target multigrid
multigrid_OBJECTS = \
"CMakeFiles/multigrid.dir/main.cpp.obj" \
"CMakeFiles/multigrid.dir/multigrid3d.cpp.obj"

# External object files for target multigrid
multigrid_EXTERNAL_OBJECTS =

multigrid.exe: CMakeFiles/multigrid.dir/main.cpp.obj
multigrid.exe: CMakeFiles/multigrid.dir/multigrid3d.cpp.obj
multigrid.exe: CMakeFiles/multigrid.dir/build.make
multigrid.exe: CMakeFiles/multigrid.dir/linkLibs.rsp
multigrid.exe: CMakeFiles/multigrid.dir/objects1.rsp
multigrid.exe: CMakeFiles/multigrid.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable multigrid.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\multigrid.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/multigrid.dir/build: multigrid.exe
.PHONY : CMakeFiles/multigrid.dir/build

CMakeFiles/multigrid.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\multigrid.dir\cmake_clean.cmake
.PHONY : CMakeFiles/multigrid.dir/clean

CMakeFiles/multigrid.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\halin\Documents\GitHub\multigrid C:\Users\halin\Documents\GitHub\multigrid C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug C:\Users\halin\Documents\GitHub\multigrid\cmake-build-debug\CMakeFiles\multigrid.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/multigrid.dir/depend

