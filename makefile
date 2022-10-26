##########################################
#           Editable options             #
##########################################

# Compiler options
CC=g++
CFLAGS= -c -Wall -g
LDFLAGS=-g
EXECUTABLE_NAME=bb

# Folders
SRC=./Blackbox_files/src
BIN=./9_dfo
# OBJ=./Blackbox_files/unix_intermediate/obj
OBJ=./9_dfo
INC=./Blackbox_files/include

# Files
SOURCE_FILES=\
	program_options.cpp\
	blackbox_model.cpp

#########################################################
#     Fixed statements do not edit unless necessary    	#
#########################################################
EXECUTABLE_FILES = $(EXECUTABLE_NAME:%=$(BIN)/%)
OBJECT_FILES     = $(SOURCE_FILES:%.cpp=$(OBJ)/%.o)
# ^^^ A more succinct expression for $(OBJECT_FILES), using
#     http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

build: $(EXECUTABLE_FILES)

clean:
	rm -r -f $(BIN)/$(EXECUTABLE_NAME)
	@# ^^^ I don't recommend suppressing the echoing of the command using @

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean

$(EXECUTABLE_FILES): $(OBJECT_FILES)
	@$(CC) $(LDFLAGS) -o $@ $^
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
	@echo "Build successful!"

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJECT_FILES): $(OBJ)/%.o: $(SRC)/%.cpp
	@echo Compiling $<
	@# ^^^ "compile a .cpp file" to create a .o file.
	@mkdir -p $(@D)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029
	@$(CC) $(CFLAGS) -I$(INC) -o $@ $<
	@# ^^^ Use $(CFLAGS), not $(LDFLAGS), when compiling.