from __future__ import division  #force floating point division
import numpy as np
import math
import collections
import logging
import sys
import os
import time
import math_operation as math_op
import miscellaneous_helper as misc

                                     
def recoverDomainWeightedGainViaMedianBalancing(value_by_signal_domain, 
                                                list_domain,
                                                initial_domain_exponent_power,
                                                max_iteration,
                                                delta_exponent_power_tolerance,
                                                output_directory, 
                                                print_to_stdout = 0,
                                                logging_level = logging.DEBUG):
    """Recovers a domain weighted gain to correct signal via median balancing
       
        This function finds a scalar that is the exponent power for a domain 
       exponentiated weight.  This weight is used for gain correction on a signal 
       that decays exponentially over time e.g. a seismic trace.  
       
       If the given function is f(x), then the gain corrected function is 
         g(x) = x^(power) * f(x)
       
       This functions finds the power value in the above equation.
         
       The algorithm is described by Jon  Clarebout and Zhiming Li in their 
       Stanford Exploration Project Report 42 article 'Definition of Time
       Gain Power'
         http://sepwww.stanford.edu/theses/sep42/
       
       Keyword arguments:  
       value_by_signal_domain  -- values of signal corresponding  to list_domain
       list_domain -- assumed to be sorted in ascending order
       initial_domain_exponent_power -- initial guess at what the exponent power 
         should be in the domain exponent weighting   
       max_iteration -- the maximum number of iterations allowed 
       delta_exponent_power_tolerance -- algorithm will exit if the change in exponent 
         values is smaller than this tolerance   
         
       Output:
         NamedTuple "DomainExponentPowerFromMedianMatching"
            exit_code -- 0 = tolerance reached 
                         1 = max iterations reached
                         2 = median of partition 2 is zero
            domain_exponent_power
            iteration_count
            initial_domain_exponent
            delta_exponent_power_tolerance   
    """
    assert list_domain.ndim == 1
    assert value_by_signal_domain.ndim == 2
    
                #setup logging
    unique_to_function_call_logger = misc.createUniqueToFunctionCallLogger()         
                                        
                #-->log to file
    current_function_name = misc.getCallingFunctionName()                        
    log_file = output_directory + \
               os.path.sep + \
               current_function_name + \
               unique_to_function_call_logger.name + \
               ".log"            
    file_console_handler = logging.FileHandler(log_file)
    misc.setupHandlerAndAddToLogger(file_console_handler, 
                                    unique_to_function_call_logger,
                                    logging_level)      
 
                #-->log to stdout
    if (print_to_stdout):            
        stdout_console_handler = logging.StreamHandler(sys.stdout)
        misc.setupHandlerAndAddToLogger(stdout_console_handler, 
                                        unique_to_function_call_logger,
                                        logging_level)           

    
                #set up return value
    output = collections.namedtuple('DomainExponentPowerFromMedianMatching', 
                                          ['exit_code', 
                                           'domain_exponent_power',
                                           'iteration_count',
                                           'initial_domain_exponent'
                                           'delta_exponent_power_tolerance',
                                           'iteration_information']) 
    output.exit_code = 'uninitialized_exit_code'
    output.domain_exponent_power = initial_domain_exponent_power;
    output.iteration_count = 0
    output.initial_domain_exponent_power = initial_domain_exponent_power   
    output.delta_exponent_tolerance = delta_exponent_power_tolerance
    
                #set up values for algorithm
                #--> set up two partitions              
    num_sample = np.size(list_domain)
    half_num_sample = math.floor(num_sample/2)
    partition_one_domain = list_domain[0:half_num_sample] 
    partition_two_domain = list_domain[half_num_sample:num_sample]          
                #--> set up median ratio scaling factor using partion's endpoints
    t_a = partition_one_domain[0]
    t_b = partition_one_domain[-1]
    t_c = partition_two_domain[0]
    t_d = partition_two_domain[-1]
    step_size_scaling = np.log( (t_c/t_b) * (t_d/t_a) ) 
    abs_value_by_signal_domain = np.abs(value_by_signal_domain)

                #setup iteration information capture
    output.iteration_information = collections.namedtuple('IterationInformation',
                                                   ['iteration_count',
                                                    'delta_domain_exponent_power',
                                                    'domain_exponent_power',
                                                    'median_1_x',
                                                    'median_1_y',
                                                    'median_1_location_percentile',
                                                    'median_2_x',
                                                    'median_2_y',
                                                    'median_2_location_percentile',
                                                    'rate',
                                                    'semblance'])
    for name in output.iteration_information._fields:
        setattr(output.iteration_information,name,-777*np.ones(max_iteration))            
    
                #log inputs to the algorithm
    unique_to_function_call_logger.info("------------------------------------------------------------------")
    unique_to_function_call_logger.info("initial_power=%g, max_iteration=%g, tolerance=%g, partition_1=[%g,%g], partition_2=[%g,%g]  " % (initial_domain_exponent_power,
                                                                max_iteration,
                                                                delta_exponent_power_tolerance,
                                                                t_a,
                                                                t_b,
                                                                t_c,
                                                                t_d))
    unique_to_function_call_logger.info("------------------------------------------------------------------")
    unique_to_function_call_logger.info("iteration | delta | exponent_power |           median_1(x,y)    [Pctl]       |           median_2(x,y)    [Pctl]        |    rate    |  semblance   ")

    delta_domain_exponent_power = np.inf 
    for iteration_count in range(max_iteration):    
        if  np.abs(delta_domain_exponent_power) < delta_exponent_power_tolerance:
            output.exit_code = 0
            truncateUninitializedIterationInformation(output.iteration_information, output.iteration_count)
            return output
        
                        #add comment on what you are doing?
        weighted_abs_value_by_signal_domain = math_op.weightSignalByDomainExponentiated(
                                                list_domain,
                                                abs_value_by_signal_domain,
                                                output.domain_exponent_power)
        partition_1_median_point = math_op.findLowerMedianTraceIndexTimeIndexValueOfFamilyOfTrace(
                                             weighted_abs_value_by_signal_domain[:,0:half_num_sample]) 
        partition_2_median_point = math_op.findLowerMedianTraceIndexTimeIndexValueOfFamilyOfTrace(
                                             weighted_abs_value_by_signal_domain[:,half_num_sample:num_sample])
        partition_1_median = partition_1_median_point.value
        partition_2_median = partition_2_median_point.value
        if partition_2_median == 0:
            output.exit_code = 2
            output.domain_exponent_power = 0
            truncateUninitializedIterationInformation(output.iteration_information, output.iteration_count)
            return output 
        log_partition_median_ratio = np.log(partition_1_median/partition_2_median)
        delta_domain_exponent_power =  log_partition_median_ratio/step_size_scaling
                    # break into helper function, throw divide by zero exception
        
        domain_exponent_power = output.domain_exponent_power + delta_domain_exponent_power
        new_old_domain_exponent_power_ratio = domain_exponent_power / output.domain_exponent_power  
        output.domain_exponent_power = domain_exponent_power
          
                            #calculate iteration info
        median_1_location = partition_one_domain[partition_1_median_point.time_index] 
        median_2_location = partition_two_domain[partition_2_median_point.time_index]  
        median_1_location_percentile =  math_op.calculateValuePercentile(t_a,
                                                                           t_b,
                                                                           median_1_location)
        median_2_location_percentile =  math_op.calculateValuePercentile(t_c,
                                                                           t_d,
                                                                           median_2_location)
        semblance = math_op.calculateSemblance(weighted_abs_value_by_signal_domain)
 
                            #store and log iteration info
        output.iteration_information.iteration_count[iteration_count] = iteration_count
        output.iteration_information.delta_domain_exponent_power[iteration_count] = delta_domain_exponent_power
        output.iteration_information.domain_exponent_power[iteration_count] = output.domain_exponent_power
        output.iteration_information.median_1_x[iteration_count] =    median_1_location
        output.iteration_information.median_1_y[iteration_count] = partition_1_median_point.value
        output.iteration_information.median_1_location_percentile[iteration_count] = median_1_location_percentile
        output.iteration_information.median_2_x[iteration_count] =    median_2_location
        output.iteration_information.median_2_y[iteration_count] = partition_2_median_point.value
        output.iteration_information.median_2_location_percentile[iteration_count] = median_2_location_percentile        
        output.iteration_information.rate[iteration_count] = new_old_domain_exponent_power_ratio
        output.iteration_information.semblance[iteration_count] = semblance
        
        unique_to_function_call_logger.info("%3d | %15g | %10g | (%10g,%15g)[%10.6g]| (%10g,%15g)[%10.6g] | %10.6g |%10.6g" % (iteration_count, 
                                                         delta_domain_exponent_power, 
                                                         output.domain_exponent_power,
                                                         median_1_location,
                                                         partition_1_median_point.value,
                                                         np.round(median_1_location_percentile,3),
                                                         median_2_location,
                                                         partition_2_median_point.value,
                                                         np.round(median_2_location_percentile,3),
                                                         np.round(new_old_domain_exponent_power_ratio,6),
                                                         np.round(semblance,6)))    
        output.iteration_count = output.iteration_count + 1

    output.exit_code = 1
    truncateUninitializedIterationInformation(output.iteration_information, output.iteration_count)
    return output                                                                      
                
        
                    
def truncateUninitializedIterationInformation(iteration_information, iteration_count):    
    for name in iteration_information._fields:         
        value_array = getattr(iteration_information, name) 
        setattr(iteration_information,name,value_array[0:iteration_count])                                                                    
               
if __name__ == '__main__':
    pass



               
         
 
 
