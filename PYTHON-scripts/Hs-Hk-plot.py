"""
The program has a dual scope: 

1. Drawing and saving 4 plots in pdf format:  

    (a) Resolution & Relevance plot: same colors are indicative of Hs-Hk points come out from same number of retained sites and different mappings.
        Moreover a zoom of the region of interest (where the slope is close to -1) is also shown on the same plot;      
    (b) Zoom of Relevance & Resolution curve (Hk-Hs) in the windows where the slope is -1; 
    (c) Slope against an increasing index (1 to N points); 
    (d) Histogram of Frequencies, that is the number of sites with more occourrences. 


2. Calculating the optimal number of sites of a biomolcule starting from an atomistic trajectory, 
   such that the information is the highest possible, after decimating atoms. It corresponds to the number of sites 
   with more occourrences obtained in the (d) plot. 
"""

                   #####################################################################################
#####################  1. Importing libraries, parsing arguments, and checking if errors are present  ####################################
                   #####################################################################################


##### In remoto per fare in modo di salvare il file PDF serve questo:

import matplotlib
matplotlib.use('Agg')

#####

# 1.1 Importing main libraries 
#import matplotlib
import math 
import os 
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

from operator import itemgetter
from collections import Counter

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


# 1.2 Finding the path of the  main folder (usually "PrOpRe") after searching for 'PYTHON-scripts' folder. 
#     Then, add /lib in order to find our libraries. 

desired_folder_name = "PYTHON-scripts"
current_directory = os.getcwd()
desired_path = None

while True:
    if desired_folder_name in os.listdir(current_directory):
        desired_path = current_directory
        break
    elif current_directory == os.path.dirname(current_directory):
        print("ERROR. 'PYTHON-script' folder has not been found. Please, check it out...\n")
        quit()
    else:
        current_directory = os.path.dirname(current_directory)

python_modules_path = desired_path + "/lib"
sys.path.append(python_modules_path)


# 1.3 Importing user-libraries 
from inp_out import * 
from general import *
from check_errors import * 
from plot_func import * 


# 1.4 Input Arguments -------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False) 

group_in=parser.add_argument_group("Required Arguments") 

group_in.add_argument('option',        help = argparse.SUPPRESS)                                                        # Mandatory			                                             
group_in.add_argument('-f', '--file',  dest='DataFile',  action='store', metavar = 'FILE', help = argparse.SUPPRESS)    # Mandatory
group_in.add_argument('-h', '--help',          action='help',        help = argparse.SUPPRESS)                                  # Optional
group_in.add_argument('-d', '--DensityPoints', dest='DensityPoints', metavar = 'INT',  help = argparse.SUPPRESS)                # Optional 
group_in.add_argument('-s', '--SlopeRange',    dest='SlopeRange',    metavar = 'STR',  help = argparse.SUPPRESS)                # Optional 
group_in.add_argument('-w', '--NumberWindows', dest='NumberWindows', metavar = 'INT',  help = argparse.SUPPRESS)                # Optional 
# --------------------------------------------------------------------------------------------------------------------------------------

    

if __name__ == '__main__':
    
    # 1.5 Printing on terminal that this code is running 
    print("\n####################################################\n")
    print("'{}' running...".format(os.path.basename(sys.argv[0]))) 
    print("---------------------------------\n")
    
    
    # 1.6 Printing help message if the script does not have any arguments  
    if len(sys.argv)==1:
        print_usage_HsHkplot()  
        quit()
        
    
    # 1.7 Printing help message if the script does not present valid arguments: in particular, the first argument must be 'density' or 'bin' according the user choice. 
    checking_accepted_tasks_HsHkplot()
    check_argv_errors_HsHkplot() 
    
    
    # 1.8 Printing error and help message if code presents not allowed arguments
    checking_valid_arguments_HsHkplot(parser)


    # 1.9 Parsing arguments
    args           = parser.parse_args()
    
    option         = args.option            # Mandatory
    DataFile       = args.DataFile          # Mandatory 
    DensityPoints  = args.DensityPoints     # Optional ('density' option only)
    SlopeRange     = args.SlopeRange        # Optional
    NumberWindows  = args.NumberWindows     # Optional ('bin' option only)
    
    
    # 1.10 Checking if mandatory files are present 
    mandatory_files_present_HsHkplot(DataFile)


    # 1.11 Printing, on terminal which option has been set between 'density' and 'bin'
    print("\n● '{}' option set. 10% completed.\n".format(sys.argv[1]))
    
    
    
                       ##################################################################
    #####################  2. Reading 'Hs-Hk-Nsites.txt' File and checking for errors  ####################################
                       ##################################################################
    
    """
    Reading 'Hs-Hk-Nsites.txt' file in order to split the Resolution (Hs), Relevance (Hk) and the number of sites (N)
    in the 1st, 2nd and 3rd row, respectively. Morover, eight error-checks are performed in order to ensure that the file 
    contains the right elements. 
    """
    
    
    # 2.1 Opening 'DataFile' 
    f_data = open(DataFile)
    
    
    # 2.2 Checking that 'DataFile' is actually found
    checking_file_found(DataFile)
    
    
    # 2.3 Checking that 'DataFile' is not empty
    check_empty_file(DataFile)
    
    
    # 2.4 Checking that 'DataFile' contains exactly 3 rows
    checking_three_rows(f_data, DataFile)
    
    
    # 2.5 Checking that each row is not empty
    checking_empty_rows_HsHkplot(f_data, DataFile)
    
    
    # 2.6 Checking that the 1st row (Resolution) contains float/int numbers, and returns the Resolution (Hs)
    Hs = checking_1st_row_int_float(f_data)
    
    
    # 2.7 Checking that the 2nd row (Relevance) contains float/int numbers, and return the Relevance (Hk)
    Hk = checking_2nd_row_int_float(f_data)
    
    
    # 2.8 Checking that the 3rd row (Nsites) contains only int numbers, and returns Nsites (N)
    N  = checking_3rd_row_int(f_data)
    
    
    # 2.9 Checking that each row contains the same number of elements 
    checking_number_elements(Hs, Hk, N)   
    
      
    # 2.10 If everything went fine, printing on terminal that the file was correctly read, without errors
    print("\n● '-f/--file {}' set. File containing the information about Hs, Hk and Nsites correctly read... 20% completed.\n".format(os.path.basename(DataFile))) 
    
    
    
    
                       ################################################################################################
    #####################  3. Checking for errors for optional arguments (DensityPoints, NumberWindows, SlopeRange)  ####################################
                       ################################################################################################
    
    """
    Checking for errors for optional arguments, that are 'DensityPoints', 'NumberWindows', and 'SlopeRange'. 
    Each function has been analysed in details in the user library 'lib/check_errors.py'
    In particular: 
       a) 'DensityPoints' flag must be an integer number higher than 0 and it is accepted only for 'density' option (details in 3.1);
       b) 'NumberWindows' flag must be an integer number higher than 0 and it is accepted only for 'bin' option (details in 3.2); 
       c) 'SlopeRange' flag can be defined as percentage or the 'closest' string is also accepted (details in 3.3).                  
    """
    
    
    # 3.1 Checking if the optional argument "DensityPoints" is set. An error occurs if a float, string or negative (incl. zero) is inserted. Only an integer is accepted. 
    #     Moreover, this argument is accepted only if 'density' option is set. 
    DensityPoints = checking_errors_DensityPoints_opt_arg(DensityPoints)  
    
    
    # 3.2 Checking if the optional argument "NumberWindows" is set. An error occurs if a float, string or negative (incl. zero) is inserted. Only an integer is accepted. 
    #     Moreover, this argument is accepted only if 'bin' option is set. 
    NumberWindows = checking_errors_NumberWindows_opt_arg(NumberWindows) 
    
    
    # 3.3 Checking if the optional argument "SlopeRange" is set. Two ways of defining this argument are feasible:
    #         A) If 'closest' string is set [-s/--SlopeRange closest], then the closest value of slope to -1 is token.
    #         B) If a percentage is defined, then the rightmost value of slope in the range (-1 ± <INT>%) is token.   
    #     This optional argument can be used for both 'bin' and 'option' arguments.   
    SlopeRange = checking_errors_SlopeRange_opt_arg(SlopeRange)     
     
    
    
    
                       ########################################################
    #####################  4. Printing summary of arguments and option used  ####################################
                       ########################################################
    """
    A summary of arguments and options is printed on screen in order to check that everything is going fine. 
    Moreover, it is useful for the user as double check that the arguments and options used in the code are exactly what requested.                    
    """
    
    # 4.1 Printing a summary for the option and the arguments used in this code 
    print_summary_HsHk_plot(DataFile, DensityPoints, NumberWindows, SlopeRange)
 


   
                      #################################################################
    #####################  5. Computing the average values of Resolution & Relevance  ####################################
                      ################################################################# 
    
    
    if (option.strip() == "density"): 
        
        """
        IF 'DENSITY' OPTION IS SET (based on same density points)
        -----------------------------------------------------------   
                 
                       
        In order to compute the average values of Resolution (Hs_avg) and Relevance (Hk_avg), let us take a number 'X' of intervals
        such that, each one contains the same number of points (50 is the default value). If using the default values of "ResRel.py" program 
        the total number of points is about 10'000 (check such value calculating len(Hs)), therefore the number of Windows is about 200.
    
        For instance: tot_points (len(Hs))=10350.  Windows=50. => N_intervals = len(Hs)/Windows = 10350/50 = 207. 
        'N_intervals' is the number of intervals having the same density of points (50 as default). 
                       
        In each one of these intervals it is possible possible to compute:
                       
            A) Hs_avg: arithmetic average of the Resolution points (Hs) that fall in such interval.  
            B) Hk_avg: arithmetic average of the Relevance points (Hk) that fall in such interval.                                                       
        """
    
   
        
        # 5.1 Defining the total number of points (tot_points) a storing the values of Hs, Hk, Nsites in 'Hs_Hk_N' list of list.  
        #     Hs_Hk_N = [[Hs1, Hk1, N1], [Hs2, Hk2, N2], [...]]
        
        tot_points = len(Hs)                
        Hs_Hk_N    = [[Hs[i], Hk[i], N[i]] for i in range(tot_points)]
            
      
        # 5.2 Sorting 1st column, that is the Resolution (Hs), in order to have ascensing values of Relevance.  
        Hs_Hk_N = sorted(Hs_Hk_N, key=itemgetter(0)) 	                                            

    
    
        # 5.3 Defining the number of intervals having the same density of points (N_intervals). 
        #     The latter is defined as the number down to the nearest integer bacause having float interval is meaningless.   
        N_intervals = math.floor(tot_points/DensityPoints)    
    
    
        # 5.4 Computing 'Hs_avg' and 'Hk_avg' after calculating 'Hs_interval' and 'Hk_interval' that are the values of Hs and Hk in every interval.
        #     'np.mean' allows to calculate the average value of an np.array. 
        Hs_avg = []
        Hk_avg = []
    
        for i in range(N_intervals):         
            Hs_interval = [x[0] for x in Hs_Hk_N[i*DensityPoints : (i*DensityPoints+DensityPoints)]]
            Hs_interval = np.array(Hs_interval)
            Hs_interval = np.mean(Hs_interval)
            Hs_avg.append(Hs_interval)
        
        
            Hk_interval = [x[1] for x in Hs_Hk_N[i*DensityPoints : (i*DensityPoints+DensityPoints)]]
            Hk_interval = np.array(Hk_interval)
            Hk_interval = np.mean(Hk_interval)
            Hk_avg.append(Hk_interval)  
    
        # 5.5. Inseriting the value "0" as first element of the np.array (position "0") for both 'Hs_avg' and 'Hk_avg'
        Hs_avg.insert(0,0) 
        Hk_avg.insert(0,0)   
    
  
    ############################################################################################################################
    
    if (option.strip() == "bin"):
       
        """
        IF 'BIN' OPTION IS SET (based on same lenght bin)
        --------------------------------------------------
                       
        The calculation of the average values of Resolution (Hs_avg) and Relevance (Hk_avg) can be also done taking in account 
        a bin with identical lenght. In practice, it is possible to divide the Resolution (Hs) in a number 'NumberWindows' of windows
        such that each one has identical lenght (default value is 50 intervals). Do not care about the number of points that fall in that interval: 
        indeed, the density of points is different.
                       
        Taking in account that the Resolution goes between 0 and 1, then the 'bin' is equals to 1/NumberWindows intervals. 
        If using the default value, that is 50 intervals, the bin is 1/50 = 0.02. In practice, the first windows goes between 0 and 0.02, 
        the second windows goes between 0.02 and 0.04 and so on, until 0.98-1.00. 
                       
        In each of such intervals it is possible to compute:
                       
            A) Hs_avg: central value of the interval. Example: interval: 0.70-0.72 => Hs_avg = 0.71 
            B) Hk_avg: arithmetic average of the Relevance points (Hk) corresponding at the Resolution (Hs) points in each interval.   
                       
        Note that there could be some intervals, especially for low values of Resolution, that do not contain points. 
        In such cases, 'Hs_avg' is the central value of the interval, while 'Hk_avg' is 'NaN' (not a number). 
        Therefore, we define the function 'nan_helper(y)' whose scope is to fill 'nan' with a linear interpolation of points. 
        Example:  y=[1, nan, nan, 2, 3] => y=[1, 1.33, 1.66, 2, 3]                                                    
        """
                       
         
                      
        # 5.a Defining the total number of points (tot_points) a storing the values of Hs, Hk, Nsites in 'Hs_Hk_N' list of list.  
        #     Hs_Hk_N = [[Hs1, Hk1, N1], [Hs2, Hk2, N2], [...]]
    
        tot_points = len(Hs)
        Hs_Hk_N    = [[Hs[i], Hk[i], N[i]] for i in range(tot_points)]
        
        
        # 5.b Defining the lenght of each interval (bin), namely the ratio between the Resolution lenght and the number of windows (number_windows)
        bin = 1/NumberWindows 
    
    
        # 5.c Computing 'Hs_avg' and 'Hk_avg' after calculating 'Hs_interval' and 'Hk_interval' that are the values of Hs and Hk in each interval.  
        #     The former is the central value of the interval; the latter is the arithmetic average of the Relevance points (Hk) 
        #     corresponding at the Resolution (Hs) points in each interval. 
    
        Hs_avg = []
        Hk_avg = []
    
        for i in range(NumberWindows):
            Hs_avg.append(bin*(i+(1/2)))

            Hk_interval = [x[1] for x in Hs_Hk_N if i*bin <= x[0] < i*bin + bin]  
            Hk_interval = np.array(Hk_interval)
            Hk_interval = np.mean(Hk_interval)      
            Hk_avg.append(Hk_interval)     
        
   
        # 5.d As there could be some intervals, especially for low values of Resolution, that do not contain points, then 'Hk_avg' is 'NaN' (not a number).
        #     Therefore, we define the function 'nan_helper(y)' whose scope is to fill 'nan' with a linear interpolation of points. 
        #     Example:  y=[1, nan, nan, 2, 3] => y=[1, 1.33, 1.66, 2, 3]   
   
   
        ## 5.d.1 Transforming 'Hk_avg' list in np.array 
        Hk_avg = np.array(Hk_avg) 
   
        ## 5.d.2 Defining nan_helper(y) for substituting 'nan' values with the linear interpolation of the array values. Ex. y=[1,nan,nan,2,3] => y=[1, 1.33, 1.66, 2, 3] 
        def nan_helper(y): 
            return np.isnan(y), lambda z: z.nonzero()[0]
       
        ## 5.d.3 Substituting 'nan' values with the linear interpolation of the array values.  
        nans, x      = nan_helper(Hk_avg)
        Hk_avg[nans] = np.interp(x(nans), x(~nans), Hk_avg[~nans])
    
        ## 5.d.4 Re-transforming 'Hk_avg' np.array in list again, for plotting afterwards 'Hs-Hk' data 
        Hk_avg = list(Hk_avg)
    
    
    
                       #############################################################
    #####################  6. Computing the derivative of Resolution & Relevance  ####################################
                       #############################################################
      
                       
    if (option.strip() == "density"): 
        
        """
        IF 'DENSITY' OPTION IS SET (based on same density points)
        -----------------------------------------------------------
  
        After computing the average values for Resolution and Relevance (Hs_avg and Hk_avg), it is possible to compute the slope of the line 
        that passes for each couple of two consecutive points calculating Delta(Y)/Delta(X). According with the value of 'N_intervals'
        the number of Hs_avg-Hk_avg points is the same, therefore it is possible to calculate a number 'N_intervals-1' of slopes.
                       
        For instance, if N_intervals assumes a value of 201, there will be 201 Hs_avg-Hk_avg points, and thus 200 values of slopes 
        each couple of consecutive points. 
                       
        slope = u = DeltaY / DeltaX =  
                  = (Hk_avg[i]-Hk_avg[i-1])/(Hs_avg[i]-Hs_avg[i-1]) 
                       
        Please, note that two particular cases are possible: indeed, especially for high values of resolution, close to 1, there is, often, a superposition of points:
                       
            a) Delta(Y) = 0 && Delta(X) not 0  ==>  It means that two consecutive values of Relevance are identical. In this case, slope = 0, that is not correct.
                                                    In this case, thus, we substitute slope=0 with slope = NaN. 
                       
            b) Delta(Y) = 0 = Delta(X)         ==>  It means that two consecutive values for both Relevance and Resolution are identical. In this case, 
                                                    the slope does not assume an acceptable value (0/0). In this case slope = NaN.     
                       
        Therefore, we use the function 'nan_helper(y)' whose scope is to fill 'nan' with a linear interpolation of points. 
        Example:  y=[1, nan, nan, 2, 3] => y=[1, 1.33, 1.66, 2, 3]
                       
        This choice makes completely sense since these particular cases happens for Resolution values very close to 1. Indeed, we are interested
        in slopes close to -1, that occur for values of resolutions around 0.70. We use 'nan_helper(y)' function, just for trace the slope plot 
        as smooth as possible.                        
        """

    
        # 6.1 Computing the number of Hs-Hk-avg points. Such number is also equal to 'N_intervals+1' (as the element 0.0 has been added as 1st element of array in 2.5)
        number_points = len(Hs_avg) 
   
    
        # 6.2 Computing the slope for each couple of consecutive points, and storing the values in 'slope' list 
        slope = [] 
        for i in range(1, number_points):
            slope.append((Hk_avg[i]-Hk_avg[i-1])/(Hs_avg[i]-Hs_avg[i-1]))
    
    
        # 6.3 Substituting the value of 0.0 with NaN, since if slope = 0 means that two consecutive values of Relevance (Hk_avg) are identical [Delta(Hk_avg)==0]
        #     Slope = 0 is not correct beacuse is due to overlapped Relevance values.
        for i in range(len(slope)):
            if(slope[i]==0.0):
                slope[i] = np.nan
    
    
        # 6.4 Transforming the 'slope' list in np.array. 
        slope = np.array(slope) 
    
    
        # 6.5 Defining nan_helper(y) for substituting 'nan' values with the linear interpolation of the array values. Ex. y=[1,nan,nan,2,3] => y=[1, 1.33, 1.66, 2, 3] 
        def nan_helper(y):
            return np.isnan(y), lambda z: z.nonzero()[0]

    
        # 6.6 Substituting 'nan' values with the linear interpolation of the array values.  
        nans, x= nan_helper(slope)
        slope[nans]= np.interp(x(nans), x(~nans), slope[~nans])
    
    
        # 6.7 Re-transforming 'slope' np.array in list again, for plotting afterwards 'slope-index' data 
        slope = list(slope)
    
    
        # 6.8 Keeping only the slope values higher than -4. Indeed, under such value the slope plot is getting too much noisy, and becomes meaningless (we are interested in slopes around -1)  
        slope2 = []
        for i in slope:             
            if(i>-4.0):
                slope2.append(i)
            else:
                break 
        
        slope = slope2
            

    ############################################################################################################################
    
    if (option.strip() == "bin"):
       
        """
        IF 'BIN' OPTION IS SET (based on same lenght bin)
        --------------------------------------------------

        After computing the average values for Resolution and relevance (Hs_avg and Hk_avg) it is possible to compute the slope of the line 
        that passes for each couple of two consecutive points calculating Delta(Y)/Delta(X). According with the value of 'number_windows'
        the number of Hs_avg-Hk_avg points is the same, therefore it is possible to calculate a number 'number_windows' of slopes.
                       
        For instance, if 'number_windows' assumes a value of 200, there will be 200 Hs_avg-Hk_avg points, and thus 200 values of slopes 
        each couple of consecutive points. 
                       
        slope = u = DeltaY / DeltaX =  
                  = (Hk_avg[i]-Hk_avg[i-1])/(Hs_avg[i]-Hs_avg[i-1]) 
                       
        NOTE: In this case, at difference with the calculation of the slope based on same density of points, Hs_avg and Hk_avg assume always different values;
              indeed, a linear interpolation of values has been performed if NaN have been found. Therefore, the ratio between Delta(Hk_avg) and Delta(Hs_avg)
              cannot assume form of indetermination. 
        """
        number_points = len(Hs_avg)

        # 6.a Computing the slope for each couple of consecutive points, and storing the values in 'slope' list 
        slope = [] 
        for i in range(1, number_points):
            slope.append((Hk_avg[i]-Hk_avg[i-1])/(Hs_avg[i]-Hs_avg[i-1]))
  
    
    print("\n● The calculation of derivative of Resolution & Relevance correctly done... 50% completed")
    
                           ###########################################################
    #########################  7. Finding best interval where slope is close to -1  #######################################
                           ###########################################################
    
    """
    In order to find the best interval where the slope is close to -1, as first it is necessary to give a definition of such interval.
                           
    Two ways are possible: 
    
        a) in 'slope' list, taking the closest value of slope to -1. 
        b) in 'slope' list, defining a range close to -1 (default: from -1.10 to -0.90), and taking the rightmost value in such range (namely with higher Resolution)  
                           
    Both methods are equally correct and the user can choose between them. However, by default, we propose the option (b) and a range [-1.10, -0.90] in order that 
    a value of slope in such range with higher resolution is taken.                       
    """
    
    # 7.1 SlopeRange == "closest": Taking, in slope2 list the closest value of slope to -1. 
    if(SlopeRange == "closest"):
    
        u = -1 
        for i in slope: 
            closest_u = [abs(x-u) for x in slope]
    
        index_closest = closest_u.index(min(closest_u))    
        closest_deriv = slope[index_closest] 
    
    
    
    # 7.2 SlopeRange(percentage): After defining a range close to -1 (from -1.10 to -0.90 is the default), it is possible to take the value in such range with higher Resolution. 
    #                            in other words, the rightmost value of slope (the last value) in such range is saved.  
    if(SlopeRange != "closest"):        
    
        ## 7.2.1 Taking the value of slope (without % symbol named 'SlopeRange_1Part') and using 'check_Int_Float_Str' for checking if the value of slope is integer, or float
        SlopeRange_1Part = SlopeRange[:-1]  
        SlopeRange_1Part = check_Int_Float_Str(SlopeRange_1Part)
    
        ## 7.2.2 Computing 'ValueSlope' i.e., by definition of percentage, SlopeRange_1Part/100
        ValueSlope       = SlopeRange_1Part/100
    
        ## 7.2.3 Computing the left and the right extremes of intervals (the center value is -1)
        SlopeLeft        = -1 - ValueSlope
        SlopeRight       = -1 + ValueSlope
    
        ## 7.2.4 Taking the rightmost value of slope in the range defined before: at the end, the last value in such range is saved (since it overwrite previous values)   
        ##       If in such interval no value is found, then an error is returned (7.2.5)
        u = -1 
        for i in slope: 
            if(SlopeLeft < i < SlopeRight):
                index_closest = slope.index(i)   
        
        ## 7.2.5 Checking if 'index_closest' is defined: indeed, if a SlopeRange is defined in terms of percentage, if the latter is a tiny value, 
        #        the interval is tiny too, and the 'slope' list could not contain a value in such interval. 
        #        For instance: SlopeRange = 3%. => Interval = [-1.03, -0.97] => No value of 'slope' list in such interval => An error is returned.                
        try: 
            index_closest 
        except NameError:
            print("\n####################################################################################")
            print("ERROR. '-s/--SlopeRange {}' set. The range is, thus, -1 \u00B1 {}, i.e. [{}, {}]".format(SlopeRange, SlopeRange, SlopeLeft, SlopeRight))
            print("       and the rightmost value of slope in such range should be token.")
            print("       However, in such interval, no value of slope has been found.")
            print("       Please, increase the value of such range in terms of percentage, or write the string 'closest' [-s/--SlopeRange closest],")
            print("       in order that the closest value of slope to -1 is token. In the latter case no range requires to be defined.")
            print("####################################################################################\n\n")
            quit()
        
        ## 7.2.6 Saving 'closest_deriv' that corresponds at the best value slope according the method chosen (range or 'closest' value to -1)            
        closest_deriv = slope[index_closest]     

    
    print("\n● The calculation of the best interval where the slope is close to -1 correctly done... 60% completed")
    
              ##################################################################################################
    ############  8. Finding Resolution, Relevance, & the number of sites in the interval where slope is -1   ############
              ##################################################################################################

    """
    According with the 'density' or 'bin' options, three lists are computed:
        1) 'Hs_u_minus_1': List of Resolution values in the interval found in (7) where slope is -1 (the closest value or the rightmost value in a range)
        2) 'Hk_u_minus_1': List of Relevance values in the interval found in (7) where slope is -1 (the closest value or the rightmost value in a range) 
        3) 'N_u_minus_1' : List of the number of sites corresponding at such Hs & Hk in the interval found in (7) where slope is -1 (the closest value or the rightmost value in a range)  
    """


    # 8.1 The closest slope value to -1 or the righmost value (with higher Relevance) of 'slope' list in the interval selected (default: [-1.10,-0.90]) has been chosen.
    #     The definition of list changes according with the option chosen: "bin" or "density".
     
    if(option.strip() == 'density'):
        Hs_area_best_slope = [x[0] for x in Hs_Hk_N[index_closest*DensityPoints : (index_closest*DensityPoints + DensityPoints)]]
        Hk_area_best_slope = [x[1] for x in Hs_Hk_N[index_closest*DensityPoints : (index_closest*DensityPoints + DensityPoints)]]
        N_area_best_slope  = [x[2] for x in Hs_Hk_N[index_closest*DensityPoints : (index_closest*DensityPoints + DensityPoints)]] 
       
    if(option.strip() == "bin"):
        Hs_area_best_slope = [x[0] for x in Hs_Hk_N if bin*(index_closest + 1/2) <= x[0] < bin*(index_closest + 3/2)] 
        Hk_area_best_slope = [x[1] for x in Hs_Hk_N if bin*(index_closest + 1/2) <= x[0] < bin*(index_closest + 3/2)] 
        N_area_best_slope  = [x[2] for x in Hs_Hk_N if bin*(index_closest + 1/2) <= x[0] < bin*(index_closest + 3/2)]


    # 8.2 Finding the optimal number of sites, corresponding at the number of sites having highest frequency.
    #     Considering that there could be more than one 'Nsites' having same highest frequency, 
    #     in such case, the bigger value of Nsites will be token. 
    #     Example: 
    #              dict_N_area_best_slope = {"1200": 5, '1300': 35, '1400': 25, '1500': 35}
    #
    #              In this case 'Nsites'=1300 and 'Nsites'=1500 have the highest frequency => 'Nsites' = 1500 is token
    #              and it corresponds at the optimal number of sites. Let us this this criterion as convention. 
    
    ## 8.2.1 Trasforming 'N_area_best_slope' list in a dict: 'Counter' function has the scope of counting how many times 'N' is repetaed:
    ##       Example: [1,1,2,3,3,1,4,5] => {'1':3, '2':1, '3':2, '4':1, '5':1}  (1 is repeated 3 times, 2 is repeated once, and so on)
      
    dict_N_area_best_slope = Counter(N_area_best_slope)
    
    ## 8.2.2 Finding the highest frequency: Ex. dict={"1200": 5, '1300': 35, '1400': 25, '1500': 35} => max_frequency = 35 
    max_frequency = max(dict_N_area_best_slope.values())

    ## 8.2.3 Finding the 'keys' i.e. the number of sites having the highest frequency. Ex: max_frequency = 35  =>  keys = [1300,1500]
    list_Nsites_max_frequency = [int(k) for k,v in dict_N_area_best_slope.items() if (v == max_frequency)]
    
    ## 8.2.4 Saving the bigger value in the latter list: such value is the OPTIMAL NUMBER OF SITES 
    optimal_number_of_sites = max(list_Nsites_max_frequency)
    
    ## 8.2.5 Printing on screen what is the optimal number of sites 
    print(f"\n● The optimal number of sites is {optimal_number_of_sites}... 80% completed.")
    
    ## 8.2.6 Writing on a file a summary of the arguments, and the optiam lnumber of sites
    g = open("Opt-number-of-sites.txt", "w") 
    write_file_summary_HsHk_plot(g, DataFile, DensityPoints, NumberWindows, SlopeRange)
    g.write(f"The optimal number of sites is {optimal_number_of_sites}")
    g.close()
        



                        #################################################
    ######################  9. Drawing and saving plots in pdf format  ###############################
                        #################################################
    """
    This last part has the purpose of drawing 4 plots and saving them in pdf format:
        a) 1st plot: Scatter Plot, Average Hs-Hk curve, Zoom of the region of interest (where the slope is -1)    
        b) 2nd plot: Zoom of Resolution and Relevance curve (Hk-Hs) in the windows where the slope is -1      
        c) 3rd plot: Slope of derivative plot and straight line with equation u = -1 
        d) 4th plot: Histogram of frequencies, that is the number of sites with more occourrences             
    """

    # Plot 1: Scatter Plot, Average Hs-Hk curve, Zoom of the region of interest (where the slope is -1) 
    ax = plt.subplot()
    Hs_Hk_scatter_plot(ax, Hs, Hk, N, Hs_avg, Hk_avg, Hs_area_best_slope, Hk_area_best_slope) 
    plt.savefig("Reso.pdf", format="pdf", bbox_inches="tight")
    plt.close('all')

    
    # Plot 2: Zoom of Resolution and Relevance curve (Hk-Hs) in the windows where the slope is -1 
    plt.subplot()   
    Hs_Hk_scatter_plot_zoom(Hs_avg, Hk_avg, Hs_area_best_slope, Hk_area_best_slope, N_area_best_slope)
    plt.savefig("Zoom-Reso.pdf", format="pdf", bbox_inches="tight")
    plt.close('all')

    
    # Plot 3: Slope of derivative plot and straight line with equation u = -1 
    slope_plot_density_opt(slope, SlopeRange)
    plt.savefig("slope.pdf", format="pdf", bbox_inches="tight")
    plt.close()
    
    
    # Plot 4: HISTOGRAM OF FREQUENCIES: Number of sites with more occourrences
    ax = plt.subplot()
    histo_Nsites_plot(ax, N_area_best_slope)
    plt.savefig("histo_Nsites.pdf", format="pdf", bbox_inches="tight")
    plt.close()
    
    
    print("\n● Saving plots in PDF format... 95% completed.")
    print("\n● No Errors! 100% completed.")
    
    



