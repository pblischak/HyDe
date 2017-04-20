EXE = hyde_cpp
OBJ = MbRandom.o HyDe.o main.o

CXX = g++
CXXFLAGS = -Wall -O3 -g -std=c++11

# Set this to `yes` if you have OpenMP.
OPENMP ?= no

ifeq ($(strip $(OPENMP)), yes)
	CXXFLAGS += -fopenmp
endif

.PHONY : clean install uninstall

$(EXE) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)
	@printf "\nTo install hyde system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\n Removing object (*.o) files and executable...\n\n"
	@rm $(OBJ) $(EXE)

install :
	@printf "\n Copying executable to /usr/local/bin...\n\n"
	@cp $(EXE) /usr/local/bin
	@printf " To uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\n Removing executable from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/$(EXE)
