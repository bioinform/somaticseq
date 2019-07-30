#!/usr/bin/env python3


import operator as op
from functools import reduce
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vaf',     '--vafs',     type=float, nargs='*', help='VAF', default=0.05)
parser.add_argument('-dp',      '--dp',      type=int, help='WGS depth', default=50)
parser.add_argument('-deep1',   '--deep-dp1', type=int, help='combined WGS depth', default=300)
parser.add_argument('-deep2',   '--deep-dp2', type=int, help='combined WGS depth', default=380)
parser.add_argument('-varWgs',  '--variant-read-WGS', type=int, help='min reads in WGS to get a call', default=3)
parser.add_argument('-varDeep', '--variant-read-Deep', type=int, help='min reads in combined WGS to get a call', default=3)



args     = parser.parse_args()
vafs     = args.vafs
dp       = args.dp
deepCov1 = args.deep_dp1
deepCov2 = args.deep_dp2
vReads   = args.variant_read_WGS
vDeep    = args.variant_read_Deep




# Function for "n Choose r"
def nCr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom

# N = total possibilityes
N = 0
for i in range(dp+1):
    j=nCr(dp,i)
    N+=j


# Sanity check P has to be 1, and it is
P=0
for i in range(dp+1):
    P += nCr(dp,i) * ( vafs[0]**i * (1-vafs[0])**(dp-i) )


print('VAF\tP(≥1/21)\tP(HighConfRightAway)\tP(Rescued)')
for vaf in vafs:
    
    # P_notCalled is the probability that a replicate will have either 0, 1, or just 2 variant reads, i.e., not enough to confidently call a somatic mutation in a particular replicate
    P_notCalled = 0
    for i in range(vReads):
        P_notCalled += nCr(dp,i) * ( vaf**i * (1-vaf)**(dp-i) )
    
    # THe probability that there will be at least 3 variant reads in a replicate, enough to be confidently called:
    P_called = 1 - P_notCalled
    
    
    # Not Called in deeperSeq
    P1_notCalledDeep = 0
    for i in range(vDeep):
        P1_notCalledDeep += nCr(deepCov1,i) * ( vaf**i * (1-vaf)**(deepCov1-i) )
    
    P1_calledDeep = 1 - P1_notCalledDeep


    P2_notCalledDeep = 0
    for i in range(vDeep):
        P2_notCalledDeep += nCr(deepCov2,i) * ( vaf**i * (1-vaf)**(deepCov2-i) )
    
    P2_calledDeep = 1 - P2_notCalledDeep


    # Probability of variant reads found in both 300X Data Sets
    P_Rescued = P1_calledDeep * P2_calledDeep

    
    # There are 21 replicates, so the probability it will be called at least once:
    P = 1 - P_notCalled**21
    
    
    # The probability that in a Data Group with 3 replicates, the variant is called at least twice:
    # There are 4 such groups
    P_atLeast2outta3 = P_called**3 + nCr(3,2) * (P_called**2 * P_notCalled**1)
    
    # The NS group has 9 replicates:
    P_atLeast5outta9 = 0
    for i in 5,6,7,8,9:
        P_atLeast5outta9 += nCr(9,i) * (P_called**i * P_notCalled**(9-i))
    
    
    # There are 5 Data Groups, probability that 3 will have majorities, so a 5% VAF is classified as HighConf directly:
    # 1) The NS Group get it, so the other 4 can be 4, 3, or 2:
    # 2) The NS Group does not get it, so the other 4 can be 4 or 3:
    P_HighConfDirectly = \
    P_atLeast5outta9     * nCr(4,4) * (P_atLeast2outta3**4 * (1-P_atLeast2outta3)**0) + \
    P_atLeast5outta9     * nCr(4,3) * (P_atLeast2outta3**3 * (1-P_atLeast2outta3)**1) + \
    P_atLeast5outta9     * nCr(4,2) * (P_atLeast2outta3**2 * (1-P_atLeast2outta3)**2) + \
    (1-P_atLeast5outta9) * nCr(4,4) * (P_atLeast2outta3**4 * (1-P_atLeast2outta3)**0) + \
    (1-P_atLeast5outta9) * nCr(4,3) * (P_atLeast2outta3**3 * (1-P_atLeast2outta3)**1)
    
    
    
    P_gotTruth = P_HighConfDirectly + (1-P_HighConfDirectly)*P_Rescued
    
    print(vaf, P, P_HighConfDirectly, P_gotTruth, sep='\t')
    
