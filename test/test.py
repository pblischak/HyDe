from __future__ import print_function
import phyde as hd

print("\n**** Test 1: Data import. ****")
data = hd.HydeData("data.txt", "map.txt", "out", 16, 4, 50000)
print("**** Good. ****")

print("\n**** Test 2: Running test_triple(). ****")
data.test_triple('sp1', 'sp2', 'sp3')
print("**** Good. ****")

print("\n**** Test 3: Running hyde without bootstrapping. ****")
res = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000)
print(res.res[res.triples[0]])
print("**** Good. ****")

print("\n**** Test 4: Running hyde analysis with bootstrapping. ****")
res, boot = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootReps=20)
print(res.res[res.triples[0]])
print(boot.breps[boot.triples[0]])
print("**** Good. ****")

print("\n**** Test 6: Reading in HyDe results. ****")
boots = hd.HydeResult("hyde-out.txt")
print("**** Good. ****")

print("\n**** Test 7: Reading in bootstrap reps. ****")
boots = hd.Bootstrap("hyde-boot.txt")
print("**** Good. ****")

print("\n**** Test 8: Getting gamma values from each boot rep. ****")
g = boots.gamma('sp1', 'sp2', 'sp3')
print(g)
print("**** Good. ****")

print("\n**** Test 9: Testing callable Bootstrap method. ****")
aaaa = boots('AAAA', 'sp1', 'sp2', 'sp3')
print(aaaa)
print("**** Good. ****")

print("\n**** Test 10: Testing visualization. ****")
hd.viz.hist(boots, 'Gamma')
print("**** Good. ****")
