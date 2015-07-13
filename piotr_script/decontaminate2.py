#!/usr/bin/env python
from __future__ import division, print_function

import sys

print('\n ----- decontaminate.py v0.1 -----\nby Piotr Lukasik', end = '\n\n', file = sys.stderr)

Exit_message = """ERROR! CHECK YOUR INPUT PARAMETERS!
This scripts takes a mothur-generated count_table file as well as a list of blank (negative control) library names, and conducts a contaminant removal procedure as detailed in [reference]. Briefly, it removes from the table any unique genotype UNLESS (a) the ratio of the maximum relative abundance it attains in any experimental sample to the maximum relative abundance in any blank is above a threshold; and (b) in at least one experimental library it attains a threshold relative abundance. Next (c), the script removes any library which lost more than a threshold proportion of reads during steps (a) and (b), as well as all blank libraries. Then (d) the script repeats step (a) using the maximum relative abundance values calculated for the remaining experimental libraries. Finally, it outputs a .decontaminated.count_table file, suitable for downstream analysis using mothur.

Usage: decontaminate.py <count_table> <list_of_blanks> <ThresholdA; recommended value 10> <ThresholdB; recommended value 0.001> <ThresholdC; recommended value 0.3>

Parameters:
<count_table>       Tab-delimited table produced by mothur's count.seqs() command. It must include experimental as well as blank libraries.
<list_of_blanks>    Text file with names of blank (negative control) libraries, one name per line.
<ThresholdA>        Value representing the threshold parameter for filtation step (a) above; a unique genotype will be removed UNLESS the maximum relative abundance it attains in at least one experimental library is more than ThresholdA * of the maximum relative abundance it attains in any blank library. Recommended value: 10
<ThresholdB>        Value representing the threshold parameter for filtation step (b) above; a unique genotype will be removed UNLESS the maximum relative abundance it attains in at least one experimental library is more than ThresholdB. Recommended value: 0.001
<ThresholdC>        Value representing the threshold parameter for filtation step (c) above; an experimental library will be removed if it lost more than ThresholdC of the starting number of reads. Recommended value: 0.7

"""

if len(sys.argv) != 6:
   print(Exit_message)
   sys.exit()

Script, Input_count_table, List_of_blank_names, ThresholdA, ThresholdB, ThresholdC = sys.argv

print('Processing count table %s ' % (Input_count_table), end = '\n', file = sys.stderr)
print('... using blank list   %s ' % (List_of_blank_names), end = '\n', file = sys.stderr)
print('... and the following contaminant removal thresholds: A - %s; B - %s; C - %s' % (ThresholdA, ThresholdB, ThresholdC), end = '\n', file = sys.stderr)
print('(run the script without parameters for a description of these thresholds).\n', end = '\n', file = sys.stderr)





###### Reading input count_table as a list of lists. Lists correspond to columns, and indexes to rows.
###### Row #0 - headings; Column #0 - Unique genotype names; Column #1 - Total counts for given genotype
print('Reading input count_table.........................................................', end = '', file = sys.stderr)

TABLE = open(Input_count_table, 'r')
headings = TABLE.readline().strip('\t\n').split('\t')
Count_table = []

for i in range(0, len(headings)):
   Count_table.append([headings[i]])


for line in TABLE: 
    line = line.strip('\t\n').split('\t')
    Count_table[0].append(line[0])
    for i in range(1, len(headings)):
       Count_table[i].append(int(line[i]))

TABLE.close()
print('   DONE!\nReading list of blanks............................................................', end = '', file = sys.stderr)




########################################################################
###### Listing numbers of columns representing blank libraries and experimental libraries
BLANKS = open(List_of_blank_names, 'r')
Blank_list = []
for line in BLANKS:
   line = line.strip('\n')
   Blank_list.append(line)
   
   
Blank_columns = []
for column_no in range(2, len(headings)):
   if headings[column_no] in Blank_list:
      Blank_columns.append(column_no)


Experimental_columns = list(range(2, len(headings)))
for column_no in Blank_columns:
   Experimental_columns.remove(column_no)


BLANKS.close()
print('   DONE!\nConverting counts into relative abundance values..................................', end = '', file = sys.stderr)




########################################################################
###### For each column except 0, calculating sum of counts from previous rows, and saving it in a newly created last row
###### Also, creating a separate spreadsheet with library summaries

Summary_table = ['NA','NA'] # First two entries

Count_table[0].append('Sum of counts')
for column_no in range(2, len(Count_table)):
   Sum = 0
   Unique_gen_count = 0
   for row_no in range(1, len(Count_table[column_no])):
      Sum += Count_table[column_no][row_no]
      if Count_table[column_no][row_no] > 0:
         Unique_gen_count += 1
   Count_table[column_no].append(Sum)
   Summary_table.append([Count_table[column_no][0], Sum, Unique_gen_count])

   ### Converting counts into relative abundance
for column_no in range(2, len(Count_table)):
   for row_no in range(1, len(Count_table[column_no])-1):
      Count_table[column_no][row_no] /= Count_table[column_no][-1]

print('   DONE!\nCalculating maximum relative abundance values.....................................', end = '', file = sys.stderr)




########################################################################
###### Finding the maximum relative abundance of each unique sequence
###### in a blank, or in an experimental library
Count_table.append(['Max_relabund_in_blank'])
Count_table.append(['Max_relabund_in_experimental_sample'])
for row_no in range(1, len(Count_table[0]) - 1):
   max_relabund_in_blank = 0
   for column_no in Blank_columns:
      if Count_table[column_no][row_no] > max_relabund_in_blank:
            max_relabund_in_blank = Count_table[column_no][row_no]
   Count_table[len(Count_table) - 2].append(max_relabund_in_blank)
   max_relabund_in_experimental = 0
   for column_no in Experimental_columns:
      if Count_table[column_no][row_no] > max_relabund_in_experimental:
            max_relabund_in_experimental = Count_table[column_no][row_no]
   Count_table[len(Count_table) - 1].append(max_relabund_in_experimental)


print('   DONE!\nContaminant removal steps (a) and (b).............................................', end = '', file = sys.stderr)



   ### At this point, column 0 - labels; 1 - sums; column n-2 - max_relabund_in_blank; column n-1 - max_relabund_in_experimental;
   ###                row 1 - library IDs; row  n-1 - total_no_reads



########################################################################
###### Contaminant removal - two steps, both using information on max relative abundances calculated previously
###### a) Discards unique genotypes whose max_relabund_in_blank < max_relabund_in_blank * ThresholdA
###### b) Discards unique genotypes whose max_relabund_in_blank < ThresholdB
for row_no in range(len(Count_table[-1]) - 1, 0, -1): # Iterates through the list backwards!
   if Count_table[-1][row_no] < Count_table[-2][row_no] * float(ThresholdA):
      for column in Count_table:
         del(column[row_no])
   elif Count_table[-1][row_no] < float(ThresholdB):
      for column in Count_table:
         del(column[row_no])


print('   DONE!\nCalculating proportions of reads remaining in libraries...........................', end = '', file = sys.stderr)




########################################################################
###### Calculating proportions of reads remaining in libraries, saving the values in an added row
###### Appending values to Summary_table
Count_table[0].append('Proportion of reads remaining')
for column_no in range(2, len(Count_table) - 2):
   Sum = 0
   Unique_gen_count = 0
   for row_no in range(1, len(Count_table[column_no]) - 1):
      Sum += Count_table[column_no][row_no]
      if Count_table[column_no][row_no] > 0:
         Unique_gen_count += 1
   Count_table[column_no].append(Sum)
   Summary_table[column_no].append(round(Sum * Count_table[column_no][-2]))
   Summary_table[column_no].append(Unique_gen_count)


print('   DONE!\nContaminant removal step (c) .....................................................', end = '', file = sys.stderr)




########################################################################
###### Contaminant removal step (c):
###### Delete libraries which retain less than ThresholdC of the starting number of reads
###### That should include all blanks! Info on max relative abundances of genotypes in blanks is retained.
for column_no in range(len(Count_table) - 3, 1, -1): # Iterates through the list backwards!
   if Count_table[column_no][-1] < float(ThresholdC):
      del(Count_table[column_no])
   elif column_no in Blank_columns:
      del(Count_table[column_no])


print('   DONE!\nCalculating maximum relative abundance values in remaining libraries .............', end = '', file = sys.stderr)




########################################################################
###### Calculating max_relabund_in_experimental_sample using the remaining libraries,
###### Then overwriting values in the last column
for row_no in range(1, len(Count_table[-1])):
   max_relabund_in_experimental = 0
   for column_no in range(2, len(Count_table) - 2):
      if Count_table[column_no][row_no] > max_relabund_in_experimental:
            max_relabund_in_experimental = Count_table[column_no][row_no]
   Count_table[-1][row_no] = max_relabund_in_experimental


print('   DONE!\nContaminant removal step (d) .....................................................', end = '', file = sys.stderr)




########################################################################
###### Now, the final decontamination step (d):
###### Discarding unique genotypes whose max_relabund_in_blank < max_relabund_in_blank * ThresholdA
for row_no in range(len(Count_table[-1]) - 1, 0, -1): # Iterates through the list backwards!
   if Count_table[-1][row_no] < Count_table[-2][row_no] * float(ThresholdA):
      for column in Count_table:
         del(column[row_no])
         
print('   DONE!\nTranslating relative abundance values back into counts ...........................', end = '', file = sys.stderr)




########################################################################
###### Info on relative abundances not needed any more!
del(Count_table[-1])
del(Count_table[-1])


###### Converting relative abundance back into counts
###### And appending totals to new dictionary Summary_remaining. Later they are added to Summary_table.
Summary_remaining = {}
for column_no in range(2, len(Count_table)):
   Unique_gen_count = 0
   Sum = 0
   for row_no in range(1, len(Count_table[column_no])-2):
      Count_table[column_no][row_no] = int(Count_table[column_no][row_no] * Count_table[column_no][-2])
      Sum += Count_table[column_no][row_no]
      if Count_table[column_no][row_no] > 0:
         Unique_gen_count += 1      
   Summary_remaining[Count_table[column_no][0]] = [Sum, Unique_gen_count]


###### Updating Summary_table using values from Summary_remaining dictionary
for column_no in range(2, len(Summary_table)):
   if Summary_table[column_no][0] in Summary_remaining:
      Summary_table[column_no].append(Summary_remaining[Summary_table[column_no][0]][0])
      Summary_table[column_no].append(Summary_remaining[Summary_table[column_no][0]][1])
   else:
       Summary_table[column_no].append(0)
       Summary_table[column_no].append(0)      

print('   DONE!\nFormatting and exporting final count_table........................................', end = '', file = sys.stderr)




########################################################################
###### Removing the last two rows from all columns
for column_no in ([0] + list(range(2, len(Count_table)))):
   del(Count_table[column_no][-1])
   del(Count_table[column_no][-1])
   
   
###### Update totals in column #1
for row_no in range(1, len(Count_table[0])):
   row_sum = 0
   for column_no in range(2, len(Count_table)):
      row_sum += Count_table[column_no][row_no]
   Count_table[1][row_no] = row_sum


###### Selecting output file name
Input_name_bits = Input_count_table.split('.')
Output_count_table_name = ''
Output_summary_name = ''
for name_bit_no in range(0, len(Input_name_bits) - 1):
   Output_count_table_name += Input_name_bits[name_bit_no] + '.'
   Output_summary_name += Input_name_bits[name_bit_no] + '.'

Output_count_table_name += 'decontaminated.count_table'
Output_summary_name += 'decontamination.summary'
OUTPUT_COUNT_TABLE = open(Output_count_table_name, 'w')
OUTPUT_SUMMARY = open(Output_summary_name, 'w')


### Exporting data as a decontaminated count table!
for row_no in range(0, len(Count_table[0])):
   for column in Count_table:
      print(column[row_no], '\t', end = '', sep = '', file = OUTPUT_COUNT_TABLE)
   print('\n', end = '', sep = '', file = OUTPUT_COUNT_TABLE)

print('   DONE!\n', end = '', file = sys.stderr)




########################################################################
###### Calculating and printing summary statistics

###### Printing summary table to screen
print('\n----------------------------\nLibrary statistics:\n\n', end = '', file = sys.stderr)

print('Library status:            Library name;      no_reads @ start / after steps (a)+(b) / @ end;      no_unique_genotypes @ start / after steps (a)+(b) / @ end', end = '\n', file = sys.stderr)

###### Printing summaries for blanks
for column_no in Blank_columns:
   print('BLANK - REMOVED: %22s;            ' % (Summary_table[column_no][0]), '%6d / %6d (%3.2f%%) / %6d (%3.2f%%);                '  % (Summary_table[column_no][1], Summary_table[column_no][3], Summary_table[column_no][3]/Summary_table[column_no][1]*100, Summary_table[column_no][5], Summary_table[column_no][5]/Summary_table[column_no][1]*100), '%5d / %5d (%3.2f%%) / %5d (%3.2f%%)' % (Summary_table[column_no][2], Summary_table[column_no][4], Summary_table[column_no][4]/Summary_table[column_no][2]*100, Summary_table[column_no][6], Summary_table[column_no][6]/Summary_table[column_no][2]*100), sep = '', end = '\n', file = sys.stderr)

###### Printing summaries for experimental samples
for column_no in Experimental_columns:
   if Summary_table[column_no][0] in Summary_remaining:
      print('LIBRARY RETAINED: %21s;            ' % (Summary_table[column_no][0]), '%6d / %6d (%3.2f%%) / %6d (%3.2f%%);                '  % (Summary_table[column_no][1], Summary_table[column_no][3], Summary_table[column_no][3]/Summary_table[column_no][1]*100, Summary_table[column_no][5], Summary_table[column_no][5]/Summary_table[column_no][1]*100), '%5d / %5d (%3.2f%%) / %5d (%3.2f%%)' % (Summary_table[column_no][2], Summary_table[column_no][4], Summary_table[column_no][4]/Summary_table[column_no][2]*100, Summary_table[column_no][6], Summary_table[column_no][6]/Summary_table[column_no][2]*100), sep = '', end = '\n', file = sys.stderr)
   else:
      print('LIBRARY REMOVED: %22s;            ' % (Summary_table[column_no][0]), '%6d / %6d (%3.2f%%) / %6d (%3.2f%%);                '  % (Summary_table[column_no][1], Summary_table[column_no][3], Summary_table[column_no][3]/Summary_table[column_no][1]*100, Summary_table[column_no][5], Summary_table[column_no][5]/Summary_table[column_no][1]*100), '%5d / %5d (%3.2f%%) / %5d (%3.2f%%)' % (Summary_table[column_no][2], Summary_table[column_no][4], Summary_table[column_no][4]/Summary_table[column_no][2]*100, Summary_table[column_no][6], Summary_table[column_no][6]/Summary_table[column_no][2]*100), sep = '', end = '\n', file = sys.stderr)


###### printing summary table to .summary file (the same set of count values as above, but no proportions, and in tab-delimited format)
print('Library status', 'Library name', 'no_reads @ start', 'no_reads after steps (a)+(b)',  'no_reads @ end', 'no_unique_genotypes @ start', 'no_unique_genotypes after steps (a)+(b)',  'no_unique_genotypes after steps (a)+(b) @ end', sep = '\t', end = '\n', file = OUTPUT_SUMMARY)

for column_no in Blank_columns:
   print('BLANK_REMOVED', (Summary_table[column_no][0]), Summary_table[column_no][1], Summary_table[column_no][3], Summary_table[column_no][5], Summary_table[column_no][2], Summary_table[column_no][4], Summary_table[column_no][6], sep = '\t', end = '\n', file = OUTPUT_SUMMARY)


for column_no in Experimental_columns:
   if Summary_table[column_no][0] in Summary_remaining:
      print('LIBRARY_RETAINED', (Summary_table[column_no][0]), Summary_table[column_no][1], Summary_table[column_no][3], Summary_table[column_no][5], Summary_table[column_no][2], Summary_table[column_no][4], Summary_table[column_no][6], sep = '\t', end = '\n', file = OUTPUT_SUMMARY)
   else:
      print('LIBRARY_REMOVED', (Summary_table[column_no][0]), Summary_table[column_no][1], Summary_table[column_no][3], Summary_table[column_no][5], Summary_table[column_no][2], Summary_table[column_no][4], Summary_table[column_no][6], sep = '\t', end = '\n', file = OUTPUT_SUMMARY)



########################################################################
###### Done :) Printing output file names
print('\nOutput files:\n      ', Output_count_table_name, '\n      ', Output_summary_name, end = '\n\n', file = sys.stderr)