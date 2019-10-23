import pandas as pd
import numpy as np
import os
import argparse
import pdb
import warnings
import subprocess

def makeScoreFile(chrom, varbeta_prefix, frq_file_prefix, output_prefix):
    '''
    Make score files that stores the values of beta

    varbeta_prefix is the file prefix of variance of betas, file extension is .varbeta
    It has to contain the column 'VARBETA'

    frq_file_prefix contains the MAF information.
    The columns are ['CHR','SNP','A1','A2','MAF','NCHROBS']

    output_prefix is the prefix of the output score files

    The score files are chromosome separated, with columns
    ['CHR','SNP','A2','BETA_SCALED'], contain chromosome information,
    SNP ID, minor allele, the multiplier for genotype respectively
    '''
    chrom = str(chrom)
    print('Generating score files for chromosom {}'.format(chrom))
    frq_df = pd.read_csv(frq_file_prefix+chrom+'.frq',delim_whitespace=True)
    varbeta_df = pd.read_csv(varbeta_prefix+chrom+'.varbeta',delim_whitespace=True)
    Nsnps_frq = frq_df.shape[0]
    Nsnps_varbeta = varbeta_df.shape[0]
    if Nsnps_frq!=Nsnps_varbeta:
        raise ValueError('Frequent file has {} SNPs, but varbeta file has {} SNPs'.format(Nsnps_frq,Nsnps_varbeta))
    df = pd.merge(frq_df,varbeta_df,how="inner",on='SNP')
    Nsnps_merged = df.shape[0]
    if Nsnps_frq!=Nsnps_merged:
        warnings.warn('Frequent file has {} SNPs, but merged dataframe has {} SNPs'.format(Nsnps_frq,Nsnps_merged))
    # We assume the genotype associated to the .frq file is not scaled
    # Therefore we need to scale beta by the standard error of each SNP
    # which equals the 2*maf*(1-maf)
    # The centering step is done with computing phenotype
    sigma = np.sqrt(2*frq_df['MAF']*(1-frq_df['MAF']))
    score_df = pd.DataFrame(None,columns=['CHR','SNP','A1','BETA_SCALED'])
    score_df[['CHR','SNP','A1']] = frq_df[['CHR','SNP','A1']]
    beta = np.random.normal(0,np.sqrt(merged['VARBETA']))
    scaled_beta = beta/sigma
    score_df['BETA_SCALED'] = scaled_beta
    score_df.to_csv(output_prefix+'.'+chrom+'.score',index=False,sep='\t')
    return

def makeProfileFile(chrom, score_file_prefix, bed_file_prefix, output_prefix):
    '''
    Multiplying the scaled beta with genotype

    score_file_prefix is the prefix to the file that contains scaled beta
    It has columns ['CHR','SNP','A2','BETA_SCALED']

    bed_file_prefix is the prefix to the bed file containing genotype
    '''
    print('Computing profile files for chromosome {}'.format(chrom))
    chrom = str(chrom)
    subprocess.call(['plink','--noweb','--bfile',bed_file_prefix+chrom,
        '--score',score_file_prefix+chrom+'.score','2','sum','center',
        '--allow-no-sex',
        '--out',output_prefix+chrom])
    return

def makePhenoFile(profile_file_prefix, h2g, fam_file_prefix, output_prefix):
    '''
    Make files contain phenotype information

    h2g is a float between 0 and 1

    Output has columns ['FID','IID','PHENO']
    '''
    fam_df = pd.read_csv(fam_file_prefix+'.fam',delim_whitespace=True,header=None)
    phe_df = pd.DataFrame(None,columns=['FID','IID','PHENO'])
    phe_df[['FID','IID']] = fam_df[['FID','IID']]
    profile_dfs = [pd.read_csv(profile_file_prefix+str(chrom)+'.profile',delim_whitespace=True) for chrom in range(1,23)]
    scoresums = np.array([df['SCORESUM'].values for df in profile_dfs])
    phenotype = np.sum(scoresums,axis=1)
    Nfam = fam_df.shape[0]
    Nphe = phenotype.shape[0]
    if Nfam!=Nphe:
        raise ValueError('There are {} samples in .fam file, but {} in the sum'.format(Nfam,Nphe))
    phenotype += np.random.normal(0,np.sqrt(1-h2g),len(phe_df))
    phe_df['PHENO'] = phenotype
    phe_df.to_csv(output_prefix+'.phe',index=False,sep='\t')
    return

def makeSumstatFile(chrom, bed_file_prefix, pheno_file, output_prefix):
    '''
    Compute summary statistics for one chromosome

    Output .qassoc files has columns [CHR,SNP,NMISS,BETA,SE,R2,T,P]
    where NMISS: Number of non-missing individuals for this analysis
    T: t-statistic for regression of phenotype on allele count
    P: Asymptotic significance value for coefficient
    '''
    subprocess.call(['plink','--bfile',bed_file_prefix+chrom,
        '--pheno',pheno_file,
        '--assoc',
        '--allow-no-sex',
        '--out',output_prefix+'.'+chrom])
    return 

def formatSumstat(bed_file_prefix, ss_file, output_prefix):
    '''
    Format qassoc files to sLDSC format

    Output .sumstats file has columns ['CHR','N','Z','A1','A2']
    '''























