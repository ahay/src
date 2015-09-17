#! /usr/bin/env python
'Do median balancing.'

# Copyright (C) 2007-2010 Ioan Vlad
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os, sys
from rsf.prog import RSFROOT
import os, sys
from ConfigParser import SafeConfigParser
import logging
import numpy as np
import m8r
import rsf.api
import rsf.user.median_balance as mb
import plot_helper
import math_operation as math_op
import file_operation

try: # Give precedence to local version
    import ivlad, m8rex, sf
except: # Use distributed version
    import rsf.user.ivlad as ivlad
    import rsf.user.m8rex as m8rex
    import rsf.user.sf as sf




def main(par):
    output_rsf_file_name = 'test'
    treat_traces_individually = 0
    inp = par.string('inp') # input file
    out = par.string('out') # output file


# def main(par):
#
#
#     inp = par.string('inp') # input file
#     out = par.string('out') # output file
#
#     verb = par.bool('verb', False) # if y, print system commands, outputs
#     pclip = par.float('pclip',99)  # percentile clip
#
#     if pclip < 0 or pclip > 100:
#         raise m8rex.ParamOutOfRange('pclip',0,100)
#
#     prog_nm_root = os.path.join(RSFROOT, 'bin', 'sf')
#     sfquantile = prog_nm_root + 'quantile'
#     sfclip     = prog_nm_root + 'clip'
#
#     clip = ivlad.getout('sfquantile', 'pclip='+str(pclip), inp, verb, True)
#     sf.clip(inp, out, clip, verb)
#
#     return ivlad.unix_success

    ###################################################################
    #                USER CONFIG: choose run parameters
    ###################################################################
    # -->location of config file
    config_file = "user_settings.cfg"
    # -->change these values for code speedup
    logging_level = logging.DEBUG
    print_to_stdout = 1
    ###################################################################

    # read in config file
    parser = SafeConfigParser()
    list_successfully_read_file = parser.read(config_file)
    assert list_successfully_read_file, "parser failed to read file = %s" % config_file
    figure_output_settings = plot_helper.getFigureOutputSettings(parser)
    max_iteration = parser.getint('median_balance_algo', 'max_iteration')
    delta_exponent_power_tolerance = parser.getfloat('median_balance_algo', 'delta_exponent_power_tolerance')
    initial_exponent_power = parser.getfloat('median_balance_algo', 'initial_exponent_power')

    # setup output directory
    output_directory = os.path.abspath(os.path.dirname(sys.argv[0]))

    # read in rsf data file as numpy array
    #data_file = "data.rsf"
    #rsf_input = rsf.api.Input(data_file)
    rsf_input = inp
    num_time_sample = rsf_input.int("n1")
    origin_time_sample = rsf_input.float("o1")
    delta_time_sample = rsf_input.float("d1")
    num_trace = rsf_input.int("n2")
    origin_trace = rsf_input.float("o2")
    delta_trace = rsf_input.float("d2")

    domain_time_sample = math_op.calculateLinearlySpacedValues(
        origin_time_sample,
        delta_time_sample,
        num_time_sample)

    data_by_trace_time = np.zeros((num_trace, num_time_sample), 'f')
    rsf_input.read(data_by_trace_time)

    processed_data_by_trace_time = np.zeros((num_trace, num_time_sample), 'f')
    if treat_traces_individually:
        list_power_by_traceindex = np.zeros(num_trace)
        list_iterationcount_by_traceindex = np.zeros(num_trace)
        list_traceindex = range(num_trace)
        for trace_index in list_traceindex:
            trace = data_by_trace_time[[trace_index], :]

            # normalize values
            normalized_trace = math_op.divideEachRowByItsMaxAbsValue(trace)

            # setup output directory
            local_output_directory = output_directory + os.sep + "individual" + str(trace_index)
            file_operation.makeDirectoryIfNotExist(local_output_directory)

            median_balance_output = \
                mb.recoverDomainWeightedGainViaMedianBalancing(normalized_trace,
                                                               domain_time_sample,
                                                               initial_exponent_power,
                                                               max_iteration,
                                                               delta_exponent_power_tolerance,
                                                               local_output_directory,
                                                               print_to_stdout=1,
                                                               logging_level=logging.DEBUG)
            weighted_trace = math_op.weightSignalByDomainExponentiated(
                domain_time_sample,
                np.squeeze(trace),
                median_balance_output.domain_exponent_power)
            processed_data_by_trace_time[trace_index, :] = weighted_trace
            list_power_by_traceindex[trace_index] = median_balance_output.domain_exponent_power
            list_iterationcount_by_traceindex[trace_index] = median_balance_output.iteration_count

            # END FOR

        plot_helper.plotPowerByTraceIndex(list_traceindex,
                                          list_power_by_traceindex,
                                          output_directory,
                                          'power_by_trace_index',
                                          figure_output_settings)

        plot_helper.plotIterationCountByTraceIndex(list_traceindex,
                                                   list_iterationcount_by_traceindex,
                                                   output_directory,
                                                   'iterationcount_by_trace_index',
                                                   figure_output_settings)
    else:
        # treat traces as familiy of traces
        # --> setup output directory
        local_output_directory = output_directory + os.sep + "family"
        file_operation.makeDirectoryIfNotExist(local_output_directory)

        median_balance_output = \
            mb.recoverDomainWeightedGainViaMedianBalancing(data_by_trace_time,
                                                           domain_time_sample,
                                                           initial_exponent_power,
                                                           max_iteration,
                                                           delta_exponent_power_tolerance,
                                                           local_output_directory,
                                                           print_to_stdout=print_to_stdout,
                                                           logging_level=logging_level)
        weighted_trace = math_op.weightSignalByDomainExponentiated(
            domain_time_sample,
            data_by_trace_time,
            median_balance_output.domain_exponent_power)
        processed_data_by_trace_time = weighted_trace

            # output data as rsf
    #output = rsf.api.Output(output_rsf_file_name)
    output = out
    output.put("n1", num_time_sample)
    output.put("o1", origin_time_sample)
    output.put("d1", delta_time_sample)
    output.put("n2", num_trace)
    output.put("o2", origin_trace)
    output.put("d2", delta_trace)
    output.put("label1", rsf_input.string("label1"))
    output.put("label2", rsf_input.string("label2"))
    output.put("unit1", rsf_input.string("unit1"))
    output.put("unit2", rsf_input.string("unit2"))

    output.write(processed_data_by_trace_time)
    return ivlad.unix_success







if __name__ == '__main__':
    ivlad.run(balance_amplitude_via_median_balancing, cpar=['inp','out'])
    # balance_amplitude_via_median_balancing(output_rsf_file_name='family_gained.rsf',
    #                                        treat_traces_individually=0)



