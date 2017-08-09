from __future__ import print_function
import phyde as hd

print("\n**** Test 1: Data import. ****")
data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
print("**** Good. ****")

print("\n**** Test 2: Running test_triple(). ****")
data.test_triple('sp1', 'sp2', 'sp3')
print("**** Good. ****")

print("\n**** Test 3: Running test_individuals(). ****")
data.test_individuals('sp1', 'sp2', 'sp3')
print("**** Good. ****")

print("\n**** Test 4: Running hyde without bootstrapping. ****")
res = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000)
print(res.res[res.triples[0]])
print("**** Good. ****")

print("\n**** Test 5: Reading in HyDe results. ****")
res2 = hd.HydeResult("hyde-out.txt")
print("**** Good. ****")

print("\n**** Test 6: Testing ABBA-BABA on HydeResult object. ****")
print(res.abba_baba("sp1", "sp2", "sp3"))
print("**** Good. ****")

#print("\n**** Test 10: Testing visualization. ****")
#p = hd.viz.density(boots, "Gamma", "sp1", "sp2", "sp3")
#print("**** Good. ****")
