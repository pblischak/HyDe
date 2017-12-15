#!/usr/bin/env python

from __future__ import print_function
import phyde as hd
print("HyDe version: ", hd.__version__, sep='')

print("\n**** Test 1: Data import. ****")
data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
print("**** Good. ****")

print("\n**** Test 2: Running test_triple(). ****")
print(data.test_triple('sp1', 'sp2', 'sp3'))
print("**** Good. ****")

print("\n**** Test 3: Running test_individuals(). ****")
print(data.test_individuals('sp1', 'sp2', 'sp3'))
print("**** Good. ****")

print("\n**** Test 4: Running bootstrap_triple(). ****")
print(data.bootstrap_triple('sp1', 'sp2', 'sp3', 20))
print("**** Good. ****")

print("\n**** Test 5: Reading in HyDe results. ****")
res2 = hd.HydeResult("../hyde-out.txt")
print("**** Good. ****")

print("\n**** Test 6: Testing ABBA-BABA on HydeResult object. ****")
print(res2.abba_baba("sp1", "sp2", "sp3"))
print("**** Good. ****")

print("\n**** Test 7: Reading in Bootstrap results. ****")
boot = hd.Bootstrap("../hyde-boot.txt")
print(boot("Gamma", "sp1", "sp2", "sp3"))
print("**** Good. ****")

print("\n**** Test 8: Testing ABBA-BABA on Bootstrap object. ****")
print(boot.abba_baba("sp1", "sp2", "sp3"))
print("**** Good. ****")
