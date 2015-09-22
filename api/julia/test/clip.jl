#!/usr/bin/env julia

using m8r

m8r.init()
inp = m8r.input("in")

n1 = m8r.histint(inp,"n1")
n2 = m8r.leftsize(inp,1)

clip = m8r.get("clip")

print(clip)