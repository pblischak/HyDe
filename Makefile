all : src/hyde_cpp

.PHONY : install clean

src/hyde_cpp : src/main.cpp src/HyDe.cpp src/MbRandom.cpp src/Makefile
	@cd src; make

clean :
	@cd src; make clean

install :
	@cp ./src/hyde
