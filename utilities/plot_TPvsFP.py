#!/usr/bin/env python3

import sys, argparse, numpy, math
import scipy.stats as scipystats
import matplotlib.pyplot as pyplot

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',  '--input-tsv-file',  type=str, help='Input TSV file',  required=True, default=None)
parser.add_argument('-label1',  '--label1',  type=str, help='Label 1 for TrueVariants_or_False',  required=False, default='True Positives')
parser.add_argument('-label0',  '--label0',  type=str, help='Label 0 for TrueVariants_or_False',  required=False, default='False Positives')

parser.add_argument('-save',    '--save-figures', action='store_true', help='Save figs', required=False, default=False)
parser.add_argument('-prefix',  '--figure-prefix', type=str, help='fig', required=False, default='fig')
parser.add_argument('-width',   '--pic-width', type=float,  help='pic width', required=False, default=16)
parser.add_argument('-height',  '--pic-height', type=float, help='pic height', required=False, default=9)
parser.add_argument('-text',    '--extra-title-text', type=str, help='2nd line of title text', required=False, default='')

args       = parser.parse_args()

fn         = args.input_tsv_file
label0     = args.label0
label1     = args.label1
save_figs  = args.save_figures
prefix     = args.figure_prefix
w          = args.pic_width
h          = args.pic_height
extra_text = args.extra_title_text



def namestr(obj, namespace=globals() ):
    
    collection = [name for name in namespace if namespace[name] is obj]
    
    new_collection = []
    for i in collection:
        if not i.startswith('_'):
            new_collection = i
    
    return new_collection


def plot_2hist(false_values, true_values, hist_bins=None, labels=None, hist_title=None ):
    
    if not hist_bins:
        try:
            hist_bins = ( min( min(false_values), min(true_values)), max( max(false_values), max(true_values)), 50 )
        except ValueError:
            hist_bins = (0, 1, 2)
        
    bins = numpy.linspace( hist_bins[0], hist_bins[1], hist_bins[2] )
    
    if not labels:
        label0 = namestr(false_values)
        label1 = namestr(true_values)
    else:
        label0 = labels[0]
        label1 = labels[1]
    
    pyplot.hist([ false_values, true_values ], bins, color=('red', 'blue'), label=(label0,label1), density=True )
    
    pyplot.legend(loc='upper right', fontsize=18)
    pyplot.title( hist_title,  fontsize=18)
    
    pyplot.tick_params(axis='x', labelsize=16)
    pyplot.tick_params(axis='y', labelsize=16)

    return 1



def figure_display(fn_string, pre_string=prefix, fig_number=1):
    if save_figs:
        pyplot.savefig('{}.{}.{}.{}'.format(pre_string, str(fig_number), fn_string, 'pdf') )
    else:
        pyplot.show()


# Import data into array
data = numpy.genfromtxt(fn, skip_header=1, delimiter='\t')


with open(fn) as f:
    header = f.readline().rstrip().split('\t')
    judgement_idx = header.index('TrueVariant_or_False')

right_calls = data[:,judgement_idx]==1
wrong_calls = data[:,judgement_idx]==0


# Plot:
# First variable:
var1_index = header.index('ALT')+1


print('Variable', 'NaN FalsePositive', 'NaN CorrectCall', 'Mean FalsePositive', 'Mean CorrectCall', 'STD FalsePositive', 'STD CorrectCall', 'Median FalsePositive', 'Median CorrectCall', 'Min FP', 'Min TP', '10th FP', '10th TP', '20th FP', '20th TP', '30th FP', '30th TP', '40th FP', '40th TP', '50th FP', '50th TP', '60th FP', '60th TP', '70th FP', '70th TP', '80th FP', '80th TP', '90th FP', '90th TP', 'Max FP', 'Max TP', sep='\t')


##### 20-->5
for i in range(var1_index, len(header) ):
    
    print( header[i], end='\t' )
    
    is_nan = numpy.isnan( data[:, i] )
    is_numeric = ~is_nan
    
    # NaN in wrong calls, and right calls:
    print( sum(is_nan[wrong_calls]), sum(is_nan[right_calls]), sep='\t', end='\t' )
    
    # Data for the wrong calls, and right calls:
    data_is_right = data[ right_calls & is_numeric, i]
    data_is_wrong = data[ wrong_calls & is_numeric, i]
    
    vars()[header[i]+'_ALL']  = data[:, i]
    vars()[header[i]+'_True']  = data_is_right
    vars()[header[i]+'_False'] = data_is_wrong

    # ... for wrong calls, and right calls
    #print(data_is_wrong.mean(), data_is_right.mean(), sep='\t', end='\t' )
    #print(data_is_wrong.std(), data_is_right.std(), sep='\t', end='\t' )
    #print(numpy.median(data_is_wrong), numpy.median(data_is_right), sep='\t', end='\t' )
    
    # try:
        # print(data_is_wrong.min(), data_is_right.min(), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 10), numpy.percentile(data_is_right, 10), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 20), numpy.percentile(data_is_right, 20), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 30), numpy.percentile(data_is_right, 30), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 40), numpy.percentile(data_is_right, 40), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 50), numpy.percentile(data_is_right, 50), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 60), numpy.percentile(data_is_right, 60), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 70), numpy.percentile(data_is_right, 70), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 80), numpy.percentile(data_is_right, 80), sep='\t', end='\t' )
        # print(numpy.percentile(data_is_wrong, 90), numpy.percentile(data_is_right, 90), sep='\t', end='\t' )
        # print(data_is_wrong.max(), data_is_right.max(), sep='\t', end='\t' )
        
    # except ValueError:
        # pass
    
    # print('')
        


for fig_n,var_i in enumerate( header[var1_index::] ):
    
    pyplot.figure(fig_n+1)
    
    try:
        plot_2hist( vars()[var_i + '_False'], vars()[var_i + '_True'], labels=(label0, label1), hist_title=var_i )
        figure_display( var_i, prefix, '%03d' % (fig_n+1) )
        
    except ValueError:
        pass
    


i_ref = header.index('REF')
i_alt = header.index('ALT')

GC2CG = [0, 0]
GC2TA = [0, 0]
GC2AT = [0, 0]
TA2AT = [0, 0]
TA2GC = [0, 0]
TA2CG = [0, 0]

# Seperate procedure to plot nucleotide changes:
with open(fn) as f:
    
    line_i = f.readline().rstrip()
    line_i = f.readline().rstrip()
    
    while line_i:
        
        item_i = line_i.split('\t')
        
        ref, alt = item_i[ i_ref ], item_i[ i_alt ]
        status = int( item_i[judgement_idx] )
        
        if   (ref,alt) == ('G','C') or (ref,alt) == ('C','G'):
            GC2CG[status] = GC2CG[status] + 1
            
        elif (ref,alt) == ('G','T') or (ref,alt) == ('C','A'):
            GC2TA[status] = GC2TA[status] + 1

        elif (ref,alt) == ('G','A') or (ref,alt) == ('C','T'):
            GC2AT[status] = GC2AT[status] + 1

        elif (ref,alt) == ('T','A') or (ref,alt) == ('A','T'):
            TA2AT[status] = TA2AT[status] + 1

        elif (ref,alt) == ('T','G') or (ref,alt) == ('A','C'):
            TA2GC[status] = TA2GC[status] + 1

        elif (ref,alt) == ('T','C') or (ref,alt) == ('A','G'):
            TA2CG[status] = TA2CG[status] + 1
        
        line_i = f.readline().rstrip()




# pyplot.figure( fig_n+1 )
# pyplot.bar( (1,2,3,4,5,6), (GC2CG[0], GC2TA[0], GC2AT[0], TA2AT[0], TA2GC[0], TA2CG[0]  ), 0.1, color="red")
# pyplot.bar( (1.1, 2.1, 3.1, 4.1, 5.1, 6.1), (GC2CG[1], GC2TA[1], GC2AT[1], TA2AT[1], TA2GC[1], TA2CG[1]  ), 0.1, color="green")
# pyplot.legend( (label0, label1) )
# pyplot.xticks( (1,2,3,4,5,6), ('G>C', 'G>T', 'G>A', 'T>A', 'T>G', 'T>C') )
# pyplot.savefig('fig.' + '%03d' % (fig_n+1) + '_' + 'ntChange' + '.pdf')



## Plot histogram of p-scores
#plot_2hist(SCORE_False, SCORE_True, hist_bins=(0,1,21), labels=('False Positive', 'True Positive'), hist_title='')
#pyplot.yscale('log', nonposy='clip')
#pyplot.legend(loc='upper right', fontsize=14)
#pyplot.tick_params(axis='x', labelsize=14)
#pyplot.tick_params(axis='y', labelsize=14)
#pyplot.ylabel('Log (N)', fontsize=16)
#pyplot.xlabel('P', fontsize=16)

#input('Press Enter to quit. ')
