all : src/hyde hyde

.PHONY : install clean

src/hyde : src/main.cpp src/HyDe.cpp src/MbRandom.cpp src/Makefile
	@cd src; make

hyde : hyde/hyde.py hyde/__init__.py setup.py
	pip install .

clean :
	@cd src; make clean

install :
	@cp ./src/hyde
