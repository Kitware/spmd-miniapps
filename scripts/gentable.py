#! /usr/bin/python

import sys

nthreads = []
avgtime = []

times = []

testlog = open(sys.argv[1], "rb")

for line in testlog:
  if line.startswith("Running with "):
    words = line.split()
    nthreads.append(int(words[2]))
  elif line.startswith("done in "):
    words = line.split()
    times.append(float(words[2]))
  elif line.startswith("--------") and times:
    #avgtime.append(float(sum(times))/float(len(times)))
    avgtime.append(min(times))
    times = []

heading = "| nthreads | actual time (sec) | ideal time (sec) | speedup |"
print heading
print "-"*len(heading)
row_format = "|{0:>10}|{1:>19.4f}|{2:>18.4f}|{3:>9.4f}|"
for i in range(0, len(nthreads)):
  row = [nthreads[i], avgtime[i], avgtime[0]/nthreads[i],
         avgtime[0]/avgtime[i]]
  print row_format.format(*row)

