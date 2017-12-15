all : test test_threads

.PHONY: test test_threads

test :
	@printf "**** Testing run_hyde.py (full analysis). ****\n"
	run_hyde.py -i examples/snake-data.txt -m examples/snake-map.txt -n 52 -t 7 -s 8466 -o out
	@printf "\n**** Testing run_hyde.py (using triples). ****\n"
	run_hyde.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing run_hyde.py (phylip format). ****\n"
	run_hyde.py -i test/data-phylip.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing bootstrap_hyde.py. ****\n"
	bootstrap_hyde.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing individual_hyde.py. ****\n"
	individual_hyde.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing HyDe Python API (phyde module). ****\n\n"
	cd test; python test.py

test_threads :
	@printf "**** Testing run_hyde_mp.py (multithreaded; full analysis). ****\n"
	run_hyde_mp.py -i examples/snake-data.txt -m examples/snake-map.txt -n 52 -t 7 -s 8466 -o out
	@printf "\n**** Testing run_hyde_mp.py (multithreaded; using triples). ****\n"
	run_hyde_mp.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing run_hyde_mp.py (multithreaded; phylip format). ****\n"
	run_hyde_mp.py -i test/data-phylip.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing bootstrap_hyde_mp.py (multithreaded). ****\n"
	bootstrap_hyde_mp.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
	@printf "\n**** Testing individual_hyde_mp.py (multithreaded). ****\n"
	individual_hyde_mp.py -i test/data.txt -m test/map.txt -o out -tr test/triples.txt -n 16 -t 4 -s 50000
