import string, sys

def old_version():
    version = string.replace(sys.version,"+","")
    version = string.split(string.split(version)[0], ".")
    old = map(lambda x: int(x[0]), version) < [2, 2, 0]
    return old
