import sys
import os
import re

with open(sys.argv[1]) as f:
	for l in f:
		l = l.rstrip("\n")
		ll = l.split("\t")
		chrx = ll[0]
		t1 = int(ll[1])
		t2 = int(ll[2])
		mid = (t1+t2)/2
		sys.stdout.write("%s\t%d\t%d\n" % (chrx, mid, mid+1))
