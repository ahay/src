'''
Created on Feb 11, 2015

@author: mlai
'''
import os

def doDebug():
    import pydevd;pydevd.settrace()

def appendToFile(file_name, string_to_append):
    touch(file_name)
    with open(file_name,'a') as file_object:
        file_object.write(string_to_append)
    
def touch(file_name):
    with open(file_name, 'a'):
            #None argument sets the access and modified times to the current time.
        os.utime(file_name, None)
        
def printFileToStdOut(file_name):
    with open(file_name,'r') as file_object:
        print(file_object.read())
          
if __name__ == "__main__":
    my_file='/tmp/test123'
    appendToFile(my_file,'a string to append\n') 
    printFileToStdOut(my_file)       