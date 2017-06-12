from __future__ import print_function
import hyde as hd
import sys

print("\n**** Test 1: Data import. ****")
data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
print("**** Good. ****")

print("\n**** Test 2: Running test_triple(). ****")
data.test_triple('sp1', 'sp2', 'sp3')
print("**** Good. ****")

print("\n**** Test 3: Running hyde analysis with bootstrapping. ****")
hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootReps=20)
print("**** Good. ****")

print("\n**** Test 4: Reading in bootstrap reps. ****")
boots = hd.Bootstrap("hyde-boot.txt")
print("**** Good. ****")

print("\n**** Test 5: Getting gamma values from each boot rep. ****")
g = boots.gammas()
print(g)
print("**** Good. ****")
