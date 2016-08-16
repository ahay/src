'''
Created on Mar 26, 2015

@author: mlai
'''
import os
import numpy as np
import collections
import cPickle as pickle
# try except if you have rsf
import rsf.api 
import m8r
######
#rsfapi = None
#from utility import _rsf

def makeDirectoryIfNotExist(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def readInRsfFileAsNumpyArrayAndHeader(rsf_file_name):
    output = collections.namedtuple('rsf_file',
                                          ['data',
                                           'header'])
    output.header = convertRsfHeaderIntoDictionary(rsf_file_name)
    array_dimensions = getArrayDimensionsFromRsfHeader(output.header)
    output.data = np.squeeze(np.zeros(array_dimensions,'f'))
    rsf_input = rsfapi.Input(rsf_file_name)
    rsf_input.read(output.data)
    return output

def getArrayDimensionsFromRsfHeader(rsf_header_dictionary):
    output = []
    current_dimension = 1
    current_dimension_size = rsf_header_dictionary.get("n" + str(current_dimension))
    while current_dimension_size:
        output.append(current_dimension_size)
        current_dimension = current_dimension + 1
        current_dimension_size = rsf_header_dictionary.get("n" + str(current_dimension))
    return output


def save2DRsfFileaAsNumpyArray(rsf_file):
    rsf_input = rsf.api.Input(rsf_file)
    n1 = rsf_input.int("n1")
    n2 = rsf_input.int("n2")   
    data_by_trace_time = np.zeros((n1,n2),'f')
    rsf_input.read(data_by_trace_time)
    np.save(rsf_file,data_by_trace_time)


def save2DRsfFileaAsTextFile(rsf_file):
    rsf_input = rsf.api.Input(rsf_file)
    n1 = rsf_input.int("n1")
    n2 = rsf_input.int("n2")
    data_by_trace_time = np.zeros((n1,n2),'f')
    rsf_input.read(data_by_trace_time)
    np.savetxt(rsf_file, data_by_trace_time, delimiter=',')


def convertRsfHeaderIntoDictionary(rsf_file_name):
    rsf_input = rsf.api.Input(rsf_file_name)
    output = dict()

            # get dimensional info from header
    current_dimension = 1
    header_for_given_dimension = readRsfHeaderForGivenDimension(rsf_input, current_dimension)
    while header_for_given_dimension.get("n" + str(current_dimension)):
        output.update(header_for_given_dimension)
        current_dimension += 1
        header_for_given_dimension = readRsfHeaderForGivenDimension(rsf_input, current_dimension)

           # get other header info
    output.update(readRsfHeaderForDataTypeInformation(rsf_input))

    return output

def readRsfHeaderForDataTypeInformation(rsf_input):
    output = dict()
    output['in'] = rsf_input.string('in')
    output['esize'] = rsf_input.int('esize')
    output['type'] = rsf_input.type
    output['form'] = rsf_input.form    
    return output    
    
def readRsfHeaderForGivenDimension(rsf_input,dimension_int):    
    output = dict()

                #create keys
    dimension_int_str = str(dimension_int)   
    num_sample_key = "n" + dimension_int_str
    delta_sample_key = "d" + dimension_int_str
    origin_sample_key = "o" + dimension_int_str
    label_key = "label" + dimension_int_str
    unit_key = "unit" + dimension_int_str

                #get and store key value pairs
    output[num_sample_key] = rsf_input.int(num_sample_key)
    output[delta_sample_key] = rsf_input.float(delta_sample_key)
    output[origin_sample_key] = rsf_input.float(origin_sample_key)
    output[label_key] = rsf_input.string(label_key)
    output[unit_key] = rsf_input.string(unit_key)
    
    return output

def save_numpy_array_as_rsf_file(numpy_array, target_file_name):
    out = m8r.File(numpy_array)
    out.sfin()



def saveObject(obj, filename):
    output = open(filename, 'wb')
    pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)    
    output.close()
    
if __name__ == '__main__':
    
    pass
