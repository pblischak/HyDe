from __future__ import print_function
import hyde as hd
import sys

print("**** Test 1: Data import...", end='')
data = hd.HydeData("./examples/snake-data.txt", "./examples/snake-map.txt", "out", 52, 7, 8466)
print("Good. ****")

print("**** Test 2: Running test_triple()...", end='')
data.test_triple('sms', 'smc', 'smi')
print("Good. ****")

print("**** Test 3: Reading in bootstrap reps...", end='')
boots = hd.Bootstrap("hyde-boot.txt")
print("Good. ****")

print("**** Test 4: Getting gamma values from each boot rep...", end='')
g = boots.gammas()
print(g)
print("Good. ****")
