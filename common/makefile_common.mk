#COMPILER MODE C++17
CXX=g++ -std=c++20


#create folders
dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

#COMPILER & LINKER FLAGS
CXXFLAG=-O3 -mavx2 -mfma
LDFLAG=-O3

#COMMIT TRACING
COMMIT_VERS=$(shell git rev-parse --short HEAD)
COMMIT_DATE=$(shell git log -1 --format=%cd --date=short)
CXXFLAG+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\"
CXXFLAG+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"

# DYNAMIC LIBRARIES # Standard libraries are still dynamic in static exe
DYN_LIBS_FOR_STATIC=-lz -lpthread -lbz2 -llzma -lcurl -lcrypto -ldeflate
# Non static exe links with all libraries
DYN_LIBS=$(DYN_LIBS_FOR_STATIC) -lboost_iostreams -lboost_program_options -lboost_serialization -lhts

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

NAME=$(shell basename $(CURDIR))
BFILE=bin/$(NAME)
EXEFILE=bin/$(NAME)_static

# Only search for libraries if goals != clean
ifeq (,$(filter clean,$(MAKECMDGOALS)))

#################################
# HTSLIB for static compilation #
#################################
# These are the default paths when installing htslib from source
HTSSRC=/usr/local
HTSLIB_INC=$(HTSSRC)/include/htslib
HTSLIB_LIB=$(HTSSRC)/lib/libhts.a

##########################################
# Boost libraries for static compilation #
##########################################
BOOST_INC=/usr/include

# If not set by user command, search for it
BOOST_LIB_IO?=$(shell whereis libboost_iostreams | grep -o '\S*\.a\b')
ifneq ($(suffix $(BOOST_LIB_IO)),.a)
    # If not found check default path
    ifeq ($(wildcard /usr/local/lib/libboost_iostreams.a),)
        # File does not exist
        $(warning libboost_iostreams.a not found, you can specify it with "make BOOST_LIB_IO=/path/to/lib...")
    else
        # File exists, set the variable
        BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
    endif
endif

# If not set by user command, search for it
BOOST_LIB_PO?=$(shell whereis libboost_program_options | grep -o '\S*\.a\b')
ifneq ($(suffix $(BOOST_LIB_PO)),.a)
    # If not found check default path
    ifeq ($(wildcard /usr/local/lib/libboost_program_options.a),)
        # File does not exist
        $(warning libboost_program_options.a not found, you can specify it with "make BOOST_LIB_PO=/path/to/lib...")
    else
        # File exists, set the variable
        BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
    endif
endif

# If not set by user command, search for it
BOOST_LIB_SE?=$(shell whereis libboost_serialization | grep -o '\S*\.a\b')
ifneq ($(suffix $(BOOST_LIB_SE)),.a)
    # If not found check default path
    ifeq ($(wildcard /usr/local/lib/libboost_serialization.a),)
        # File does not exist
        $(warning libboost_serialization.a not found, you can specify it with "make BOOST_LIB_SE=/path/to/lib...")
    else
        # File exists, set the variable
        BOOST_LIB_SE=/usr/local/lib/libboost_serialization.a
    endif
endif

# Endif makefile goals != clean
endif

#CONDITIONAL PATH DEFINITON
desktop: $(BFILE)

olivier: HTSSRC=$(HOME)/Tools
olivier: HTSLIB_INC=$(HTSSRC)/htslib-1.15
olivier: HTSLIB_LIB=$(HTSSRC)/htslib-1.15/libhts.a
olivier: BOOST_INC=/usr/include
olivier: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
olivier: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
olivier: $(BFILE)

rgc: HTSSRC=/mnt/efs_v2/agds_methods/users/olivier.delaneau/LIBS
rgc: HTSLIB_INC=$(HTSSRC)/htslib-1.18
rgc: HTSLIB_LIB=$(HTSSRC)/htslib-1.18/libhts.a
rgc: BOOST_INC=/usr/include
rgc: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
rgc: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
rgc: $(BFILE)

laptop: HTSSRC=$(HOME)/Tools
laptop: HTSLIB_INC=$(HTSSRC)/htslib-1.10
laptop: HTSLIB_LIB=$(HTSSRC)/htslib-1.10/libhts.a
laptop: BOOST_INC=/usr/include
laptop: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
laptop: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
laptop: $(BFILE)

debug: CXXFLAG=-g -mavx2 -mfma
debug: LDFLAG=-g
debug: CXXFLAG+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\"
debug: CXXFLAG+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
debug: HTSSRC=$(HOME)/Tools
debug: HTSLIB_INC=$(HTSSRC)/htslib-1.15
debug: HTSLIB_LIB=$(HTSSRC)/htslib-1.15/libhts.a
debug: BOOST_INC=/usr/include
debug: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
debug: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
debug: $(BFILE)

static_exe: CXXFLAG=-O2 -mavx2 -mfma -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
static_exe: HTSLIB_INC=$(HTSSRC)/htslib-1.18
static_exe: HTSLIB_LIB=$(HTSSRC)/htslib-1.18/libhts.a
static_exe: BOOST_INC=/usr/include
static_exe: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
static_exe: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
static_exe: LDFLAG=-O2
static_exe: $(EXEFILE)

# static desktop Robin
static_exe_robin_desktop: CXXFLAG=-O2 -mavx2 -mfma -D__COMMIT_ID__=\"$(COMMIT_VERS)\" -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
static_exe_robin_desktop: LDFLAG=-O2
static_exe_robin_desktop: HTSSRC=/home/robin/Dropbox/LIB
static_exe_robin_desktop: HTSLIB_INC=$(HTSSRC)/htslib_minimal
static_exe_robin_desktop: HTSLIB_LIB=$(HTSSRC)/htslib_minimal/libhts.a
static_exe_robin_desktop: BOOST_INC=/usr/include
static_exe_robin_desktop: BOOST_LIB_IO=$(HTSSRC)/boost/lib/libboost_iostreams.a
static_exe_robin_desktop: BOOST_LIB_PO=$(HTSSRC)/boost/lib/libboost_program_options.a
static_exe_robin_desktop: $(EXEFILE)

dnanexus: BOOST_INC=/usr/include
dnanexus: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
dnanexus: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
dnanexus: HTSLIB_INC=/usr/local/include/
dnanexus: HTSLIB_LIB=/usr/local/lib/libhts.a
dnanexus: $(BFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(EXEFILE): $(OFILE)
	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS_FOR_STATIC)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean:
	rm -f obj/*.o $(BFILE) $(EXEFILE)
