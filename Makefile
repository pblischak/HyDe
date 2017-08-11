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
	@printf "**** Testing run_hyde.py (full analysis). ****\n"
	run_hyde.py -i examples/snake-data.txt -m examples/snake-map.txt -n 52 -t 7 -s 8466 -o out
	@printf "\n**** Testing run_hyde.py (using triples). ****\n"
	run_hyde.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing bootstrap_hyde.py. ****\n"
	bootstrap_hyde.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing individual_hyde.py. ****\n"
	individual_hyde.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing HyDe Python API (phyde module). ****\n\n"
	cd test; python test.py
