from __future__ import print_function

# fix for python2
try: input = raw_input
except NameError: pass

string = input('Input a string: ')

for char in string:
    print(char)

