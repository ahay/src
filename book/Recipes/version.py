import string, sys

def old_version():
    version = sys.version.replace("+","")
    version = version.split()[0].split(".")
    old = [int(x[0]) for x in version] < [2, 2, 0]
    return old
