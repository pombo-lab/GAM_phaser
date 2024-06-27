#!~/.guix-profile/bin/python3

#=====================================================================================================================
# Script written by Alexander Kukalev 
# This script is part of the study: Irastorza-Azcarate I,  Kukalev A,  Kempfer R et al  
# "Extensive folding variability between homologous chromosomes in mammalian cells", bioRxiv, 2024
# The script perform window calling and creates phased GAM segregation tables from the phased and nonphased coverage tables
#=====================================================================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import sys
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def remove_path_from_columns_names (df):
    list_of_columns = df.columns.tolist()
    simplified_columns = []
    for item in list_of_columns:
        file_name = item.split('/')[-1]
        file_name = file_name.split ('.')[0]
        simplified_columns.append (file_name)
    df1 = df.copy()
    df1.columns = simplified_columns
    return df1

def get_lowest_percentile (coverage_df):
    # Defined optimal threshold based on the coverage, sequencing depth and resolution
    # For resolutions below 200Kb the algorithm scan the coverage distribution from 0 to 98 percentile
    # For resolutions above 200Kb the algorithm scan the coverage distribution from 40 to 98 percentile
    import pandas as pd
    import numpy as np
    if len (coverage_df.index.names) == 1: # Check if supplied segregation table has multiindex
        coverage_df.set_index(['chrom','start','stop'], inplace = True) 
    start_values = np.array (coverage_df.index.get_level_values('start'))
    resolution = start_values [1] - start_values [0] # Find the resolution from table 
    if resolution >= 200000: # Set percentile boundaries depending on the resolution
        lowest_percentile = 40 # For resolutions above 250Kb go from 40 to 98 percentile
    else:
        lowest_percentile = 0 # For resolutions below 250Kb go from 0 to 98 percentile
    return lowest_percentile

def get_chromosomes_list (coverage_df):
    import pandas as pd
    import numpy as np
    chromosomes_list = []
    if len (coverage_df.index.names) == 1: # Check if supplied segregation table has multiindex
        coverage_df.set_index(['chrom','start','stop'], inplace = True)
    chromosome_values = np.array (coverage_df.index.get_level_values('chrom'))
    chromosomes = set (chromosome_values)
    for chromosome in chromosomes:
        if (chromosome.find ('random') == -1) and (chromosome.find ('chrUn') == -1):
            if chromosome not in ['chrX', 'chrY', 'chrM']:
                chromosomes_list.append (chromosome)
    chromosomes_list.sort()
    return chromosomes_list

def get_resolution (segregation_table):
    start_values = np.array (segregation_table.index.get_level_values('start'))
    resolution = start_values [1] - start_values [0] # Find the resolution from table 
    return resolution 

def get_segregation (coverage_per_bin, threshold):
    segregation_per_bin = coverage_per_bin > threshold
    segregation_per_bin = segregation_per_bin.astype(int)
    return segregation_per_bin
             
def get_orphan_windows (segregation_per_bin, chromosomes):
    import numpy as np
    import pandas as pd
    orphan_windows = []
    aux, counter = 0, 0
    # count orphan windows
    for chromosome in chromosomes:
        value_list = segregation_per_bin.loc[chromosome].values
        for i,window in enumerate(value_list):
            if window == 1:
                counter += 1
                if i == 0:
                    if sum(value_list[i:i+2]) == 1:
                        aux += 1
                elif i == len(value_list):
                    if sum(value_list[i-1:i+1]) == 1:
                        aux += 1
                else:
                    if sum(value_list[i-1:i+2]) == 1:
                        aux += 1     
    try:
        orphan_percentage = (aux/float(counter))*100
    except:
        orphan_percentage = 0
    return orphan_percentage

def calculate_threshold (coverage_per_bin, coverage_per_bin_log10, lowest_percentile, chromosomes):
    biggest_percentile = 98
    list_of_orphans = []
    for percentile in reversed(range(lowest_percentile, biggest_percentile + 1)):
        threshold = np.percentile (coverage_per_bin_log10, percentile)
        threshold_in_nucleotides = pow (10, threshold)
        sample_segregation = get_segregation (coverage_per_bin, threshold_in_nucleotides)
        percent_of_orphan_windows = get_orphan_windows (sample_segregation, chromosomes)
        list_of_orphans.append (percent_of_orphan_windows)
    for pos in range (len(list_of_orphans)):
        if list_of_orphans[pos]==0:
            list_of_orphans[pos]=100
    percentile_with_lowest_orphans = np.argmin(list_of_orphans)
    optimal_threshold = np.percentile (coverage_per_bin_log10, (biggest_percentile - percentile_with_lowest_orphans))
    optimal_threshold_estimate = pow (10, optimal_threshold)
    if optimal_threshold_estimate > 78:
        optimal_threshold_in_nucleotides = optimal_threshold_estimate
    else:
        optimal_threshold_in_nucleotides = 78
    return optimal_threshold_in_nucleotides

def coverage_df (non_phased_coverage, genome1_coverage, genome2_coverage, sample):
    # Build a df with phased and non-phased coverage
    coverages_df = pd.DataFrame (index=non_phased_coverage.index)
    coverages_df.loc [:, 'Non_phased_coverage'] = non_phased_coverage [sample].values
    coverages_df.loc [:, 'genome1_coverage'] = genome1_coverage [sample].values
    coverages_df.loc [:, 'genome2_coverage'] = genome2_coverage [sample].values
    return coverages_df

def get_phased_segregation_tables (path_to_the_nonphased_coverage_file, path_to_the_genome1_coverage_file, path_to_the_genome2_coverage_file, dataset_name):
    from datetime import date
    import pandas as pd
    import numpy as np
    # Load coverage tables in memory
    print ('Reading coverage tables...')
    print ('Nonphased coverage -> ' + path_to_the_nonphased_coverage_file)
    print ('Genome1 coverage -> ' + path_to_the_genome1_coverage_file)    
    print ('Genome2 coverage -> ' + path_to_the_genome2_coverage_file)    
    non_phased_coverage_long_names = pd.read_csv(path_to_the_nonphased_coverage_file, index_col=[0, 1, 2], sep='\t')
    genome1_coverage_long_names = pd.read_csv(path_to_the_genome1_coverage_file, index_col=[0, 1, 2], sep='\t')
    genome2_coverage_long_names = pd.read_csv(path_to_the_genome2_coverage_file, index_col=[0, 1, 2], sep='\t')
    # Remove path from the column names and make sure that the sample names are consistent between nonphased and phased tables
    non_phased_coverage = remove_path_from_columns_names (non_phased_coverage_long_names)
    genome1_coverage = remove_path_from_columns_names (genome1_coverage_long_names)
    genome2_coverage = remove_path_from_columns_names (genome2_coverage_long_names)
    # window calling and phasing
    chromosomes = get_chromosomes_list (non_phased_coverage) # Get the list of chromosomes
    lowest_percentile = get_lowest_percentile (non_phased_coverage) # Get the lowest percentile depending on the resolution
    non_phased_positive_windows, genome1_positive_windows, genome2_positive_windows, dual_positive_windows = 0,0,0,0
    nonphased_segregation_df = pd.DataFrame(index=non_phased_coverage.index)    
    genome1_segregation_df = pd.DataFrame(index=non_phased_coverage.index)
    genome2_segregation_df = pd.DataFrame(index=non_phased_coverage.index)
    list_of_samples = non_phased_coverage.columns.tolist()
    for sample in list_of_samples:
        import sys
        sys.stdout.write(".")
        sys.stdout.flush()
        coverages_df = coverage_df (non_phased_coverage, genome1_coverage, genome2_coverage, sample)
        # Calculate percentile from non-phased data
        coverage_per_bin = coverages_df ['Non_phased_coverage']
        coverage_above0 = coverage_per_bin[coverage_per_bin > 0].tolist()
        coverage_per_bin_log10 = np.log10(coverage_above0)
        if coverage_per_bin_log10.size != 0: # Calculate optimal threshoold for samples that have coverage 
            optimal_threshold_in_nucleotides = calculate_threshold (coverage_per_bin, coverage_per_bin_log10, lowest_percentile, chromosomes) 
            nonphased_segregation_per_sample = get_segregation (coverages_df ['Non_phased_coverage'], optimal_threshold_in_nucleotides)
            genome1_segregation_per_sample = get_segregation (coverages_df ['genome1_coverage'], optimal_threshold_in_nucleotides)
            genome2_segregation_per_sample = get_segregation (coverages_df ['genome2_coverage'], optimal_threshold_in_nucleotides)
        else: # Generate empty segregation values if sample have no coverage
            zeros_array = np.zeros (len(coverage_df.index), dtype=int)
            nonphased_segregation_per_sample = pd.Series(zeros_array, index=coverage_df.index)
            genome1_segregation_per_sample = pd.Series(zeros_array, index=coverage_df.index)
            genome2_segregation_per_sample = pd.Series(zeros_array, index=coverage_df.index)
        nonphased_segregation_df.loc [:, sample] = nonphased_segregation_per_sample.values
        genome1_segregation_df.loc [:, sample] = genome1_segregation_per_sample.values
        genome2_segregation_df.loc [:, sample] = genome2_segregation_per_sample.values
        # Calculate percent of phased windows
        non_phased_positive_windows = non_phased_positive_windows + len (coverages_df.query('Non_phased_coverage > @optimal_threshold_in_nucleotides'))
        genome1_positive_windows = genome1_positive_windows + len (coverages_df.query('genome1_coverage > @optimal_threshold_in_nucleotides and genome2_coverage <= @optimal_threshold_in_nucleotides'))
        genome2_positive_windows = genome2_positive_windows + len (coverages_df.query('genome2_coverage > @optimal_threshold_in_nucleotides and genome1_coverage <= @optimal_threshold_in_nucleotides'))
        dual_positive_windows = dual_positive_windows + len (coverages_df.query('genome1_coverage > @optimal_threshold_in_nucleotides and genome2_coverage > @optimal_threshold_in_nucleotides'))
    # Print stats
    print ('\nTotal number of non-phased positive windows ' + str (non_phased_positive_windows))
    print ('Number of genome1-only phased positive windows ' + str (genome1_positive_windows))
    print ('Number of genome2-only phased positive windows ' + str (genome2_positive_windows))
    print ('Number of dual-phased positive windows ' + str (dual_positive_windows))
    print ('Percent of genome1-only phased windows ' + str((genome1_positive_windows/non_phased_positive_windows)*100))
    print ('Percent of genome2-only phased windows ' + str((genome2_positive_windows/non_phased_positive_windows)*100))
    print ('Percent of dual phased windows ' + str((dual_positive_windows/non_phased_positive_windows)*100))
    # Save segregation tables
    resolution = get_resolution (non_phased_coverage)
    nonphased_name_of_segregation_table = dataset_name + '.nonphased.segregation_at' + str(resolution) + '.' + str(date.today()) + '.table'
    nonphased_segregation_df.to_csv (nonphased_name_of_segregation_table, sep='\t', index=True, header=True)
    genome1_name_of_segregation_table = dataset_name + '.genome1.segregation_at' + str(resolution) + '.' + str(date.today()) + '.table'
    genome1_segregation_df.to_csv (genome1_name_of_segregation_table, sep='\t', index=True, header=True)
    genome2_name_of_segregation_table = dataset_name + '.genome2.segregation_at' + str(resolution) + '.' + str(date.today()) + '.table'
    genome2_segregation_df.to_csv (genome2_name_of_segregation_table, sep='\t', index=True, header=True)
    print ('Nonphased segregation table was saved as ' + nonphased_name_of_segregation_table)
    print ('Genome1 segregation table was saved as ' + genome1_name_of_segregation_table)
    print ('Genome2 segregation table was saved as ' + genome2_name_of_segregation_table)


#========================================================================================

if len(sys.argv)==5:
    (non_phased_coverage, genome1_coverage, genome2_coverage, dataset_name) = sys.argv[1:]
else:
    print ('You need to provide the path to the nonphased coverage table and two phased coverage tables')
    sys.exit()

get_phased_segregation_tables (non_phased_coverage, genome1_coverage, genome2_coverage, dataset_name)
