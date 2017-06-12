from __future__ import print_function
import hyde as hd
import sys

print("**** Test 1: Data import...", end='')
data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
print("Good. ****")

print("**** Test 2: Running test_triple()...", end='')
data.test_triple('sp1', 'sp2', 'sp3')
print("Good. ****")

print("**** Test 3: Running hyde analysis with bootstrapping...", end='')
hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootRep=20)
print("Good. ****")

print("**** Test 4: Reading in bootstrap reps...", end='')
boots = hd.Bootstrap("hyde-boot.txt")
print("Good. ****")

print("**** Test 5: Getting gamma values from each boot rep...", end='')
g = boots.gammas()
print(g)
print("Good. ****")
