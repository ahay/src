'''
Created on Mar 5, 2015

@author: mlai
'''
from __future__ import division  #force floating point division
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import shutil
import matplotlib.cm as cm
import matplotlib
import json

def plotGroupsOfColumnDataInSubplots(data_by_row_column, 
                                   num_trace, 
                                   num_subplot_col):
                #split data into groups of columns for visualization
    list_data_by_groupofcol = np.array_split(data_by_row_column,
                                             num_trace,
                                             1)
    assert num_trace == len(list_data_by_groupofcol)
    
    num_subplot_row = math.ceil(num_trace/num_subplot_col)
    fig, axes_array = plt.subplots(int(num_subplot_row), 
                                   int(num_subplot_col), 
                                   sharex=True, 
                                   sharey=True)
    for index_group_of_col_data in xrange(num_trace):
        groupofcol = list_data_by_groupofcol[index_group_of_col_data]
        subplot_row_index = math.floor(index_group_of_col_data/num_subplot_col)
        subplot_col_index =  index_group_of_col_data % num_subplot_col
        axes = axes_array[subplot_row_index,subplot_col_index ]
        axes.plot(groupofcol)
        axes.set_title('num col = %d, max = %g abs-max = %g' % (groupofcol.shape[1], 
                                                np.max(groupofcol),
                                                np.max(np.abs(groupofcol))))
    plt.show() 
     
def plotSeismicShotRecord(list_time,
                   list_offset,
                   data_by_time_offset,
                   output_directory,
                   output_file_base_name,
                   figure_output_settings):
    figure = plt.figure()
    plt.pcolormesh(list_offset,
                   list_time,
                   np.transpose(data_by_time_offset),
                   cmap = cm.gray)
    
    
                #make figure look like typical seismic trace plo
    figure.axes[0].xaxis.set_label_position('top') 
    figure.axes[0].xaxis.tick_top()
    plt.gca().invert_yaxis()
    figure.axes[0].axis('tight')
    
    plt.xlabel('Offset (km)')
    plt.ylabel('Time (s)')
    
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file,**figure_output_settings)
    
        
def plotSignalMagnitudePhase(signal_support,
                             value_by_signal_range,
                             list_signal_description,
                             output_directory,
                             output_file_base_name,
                             figure_output_settings):
    figure, axes_array = plt.subplots(3,1)
    num_signal = value_by_signal_range.shape[0]
    for idx_signal in range(num_signal):
        value_by_range = value_by_signal_range[idx_signal,:]
        axes = axes_array[0]
        line, = axes.plot(signal_support,
                          value_by_range, 
                          label=list_signal_description[idx_signal])      
        axes.set_title("value_by_signal_range")
        legend = axes.legend(loc='upper right', shadow=True)

        fft_signal = np.fft.fftshift(np.fft.fft(value_by_range))
        axes = axes_array[1]    
        axes.plot(np.abs(fft_signal))
        axes.set_title("fft magnitude")
                
        phase = np.unwrap(np.angle(fft_signal))
        axes = axes_array[2]    
        axes.plot(phase)
        axes.set_title("fft phase")
    
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file, **figure_output_settings)
    
def saveFigureToFile(figure,
                     output_file,
                     figure_width_in = 6,
                     figure_height_in = 4,
                     figure_font_size = 6,
                     list_figure_image_type_extension = ['jpg'],
                     figure_dpi = 500):
    
    figure.set_size_inches(figure_width_in, figure_height_in)
    plt.rc('font', size=figure_font_size)
    plt.tight_layout()
    for ext in list_figure_image_type_extension:    
        plt.savefig(output_file + "." + ext,
               dpi=figure_dpi)
    plt.close(figure)  

def plotSemblanceVsIterationCount(list_iteration_count,
                  list_semblance,
                  output_directory,
                  output_file_base_name,
                  figure_output_settings):          

    figure = plt.figure()
    plt.plot(list_iteration_count,
              list_semblance)  
    plt.xlabel('iteration_count')
    plt.ylabel('semblance')
    
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file,**figure_output_settings)

def plotSemblanceVsDomainPower(list_domain_power,
                              list_semblance,
                              output_directory,
                              output_file_base_name,
                              figure_output_settings):          

    figure = plt.figure()
    plt.plot(list_domain_power,
              list_semblance)  
    plt.xlabel('domain_power')
    plt.ylabel('semblance')
    max_location_index = np.argmax(list_semblance)
    max_semblance_value = list_semblance[max_location_index]
    max_semblance_location = list_domain_power[max_location_index]
    min_location_index = np.argmin(list_semblance)
    min_semblance_value = list_semblance[min_location_index]
    min_semblance_location = list_domain_power[min_location_index]
    
    plt.title("max=(%g,%g),min=(%g,%g)" % (max_semblance_location,
                                           max_semblance_value,
                                           min_semblance_location,
                                           min_semblance_value))
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file,**figure_output_settings)
        
def plotDifferenceBetweenTwoSignal(signal_support,
                             value_by_signal_range,
                             list_signal_description,
                             output_directory,
                             output_file_base_name,
                             figure_output_settings):          
    assert value_by_signal_range.ndim == 2
    assert value_by_signal_range.shape[0] == 2
    figure, axes_array = plt.subplots(3,1)

    axes = axes_array[0]
    axes.plot(signal_support,
              value_by_signal_range[0,:])  
    axes.set_title(list_signal_description[0])

    axes = axes_array[1]    
    axes.plot(signal_support,
              value_by_signal_range[1,:])  
    axes.set_title(list_signal_description[1])
            
    axes = axes_array[2]
    abs_difference = np.abs(value_by_signal_range[0,:] - value_by_signal_range[1,:])
    axes.plot(signal_support,
              abs_difference)      
    axes.set_title("absolute difference")
    
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file,**figure_output_settings)
    
def plotPowerByTraceIndex(list_traceindex,
                          list_power_by_traceindex, 
                          output_directory,
                          output_file_base_name,
                          figure_output_settings):
    figure = plt.figure()
    plt.plot(list_traceindex,
              list_power_by_traceindex)  
    plt.xlabel('Trace Index')
    plt.ylabel('Power Gain')
    
    plt.title('avg=%g min=%g max=%g var=%g' %(np.mean(list_power_by_traceindex),
                                              np.min(list_power_by_traceindex),
                                              np.max(list_power_by_traceindex),
                                              np.var(list_power_by_traceindex)))
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file,**figure_output_settings)

def plotIterationCountByTraceIndex(list_traceindex,
                          list_iterationcount_by_traceindex, 
                          output_directory,
                          output_file_base_name,
                          figure_output_settings):
    figure = plt.figure()
    plt.plot(list_traceindex,
              list_iterationcount_by_traceindex)  
    plt.xlabel('Trace Index')
    plt.ylabel('Iteration Count')
    
    plt.title('avg=%g min=%g max=%g var=%g' %(np.mean(list_iterationcount_by_traceindex),
                                              np.min(list_iterationcount_by_traceindex),
                                              np.max(list_iterationcount_by_traceindex),
                                              np.var(list_iterationcount_by_traceindex)))
    output_file = output_directory + os.path.sep + output_file_base_name 
    saveFigureToFile(figure,output_file,**figure_output_settings)
     
        
def getFigureOutputSettings(parser):
    section = 'figure_output'
    figure_image_type_extension_comma_sep = parser.get(section,'list_figure_image_type_extension')
    list_figure_image_type_extension_readin = figure_image_type_extension_comma_sep.split(',')
    return dict(figure_width_in = parser.getfloat(section, 'figure_width_in'),
         figure_height_in = parser.getfloat(section, 'figure_height_in'),
         figure_font_size = parser.getfloat(section, 'figure_font_size'),
         figure_dpi = parser.getfloat(section, 'figure_dpi'),
         list_figure_image_type_extension = list_figure_image_type_extension_readin)
    
def copyDisplayAllDirAndJpgPhpAsIndexPhp(output_directory, display_all_dir_and_jpg_php):
    output_file = output_directory + os.path.sep + "index.php"
    shutil.copy2(display_all_dir_and_jpg_php, output_file)            