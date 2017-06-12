all : src/hyde_cpp

.PHONY : install clean clean_module test

src/hyde_cpp : src/main.cpp src/HyDe.cpp src/MbRandom.cpp src/Makefile
	@cd src; make

clean :
	@cd src; make clean

clean_module :
	rm -rf build dist HyDe.egg-info

install :
	@cd src; make install

test :
	hyde_cpp -h
	run_hyde.py -i examples/snake-data.txt -m examples/snake-map.txt -n 52 -t 7 -s 8466 -o out
	cd test; python test.py
