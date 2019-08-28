import pandas as pd
import numpy as np
import os
import argparse
import pdb

def makeScoreFiles(chrom, var_beta_file, frq_file_prefix, output_prefix):
    '''
    Make score files stores the values of beta

    var_beta_file is the file name of variance of betas.
    It has to contain the column 'VARBETA'

    frq_file_prefix contains the MAF information.
    The columns are ['CHR','SNP','A1','A2','MAF','NCHROBS']

    output_prefix is the prefix of the output score files

    The score files are chromosome separated, with columns
    ['CHR','SNP','A2','beta_scaled'], contain chromosome information,
    SNP ID, minor allele, the multiplier for genotype respectively
    '''
    chrom = str(chrom)
    print('Generating score files for chromosom {}'.format(chrom))
    frq_df = pd.read_csv(frq_file_prefix+chrom+'.frq',delim_whitespace=True)
    # We assume the genotype associated to the .frq file is not scaled
    # Therefore we need to scale beta by the standard error of each SNP
    # which equals the 2*maf*(1-maf)
    # The centering step is done with computing phenotype
    sigma = np.sqrt(2*frq_df['MAF']*(1-frq_df['MAF']))
    score_df = pd.DataFrame(None,columns=['CHR','SNP','A1','BETA_SCALED'])
    score_df[['CHR','SNP','A1']] = frq_df[['CHR','SNP','A1']]

