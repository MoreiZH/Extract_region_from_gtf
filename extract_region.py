#!usr/bin/python3
# -*- coding: utf-8 -*- 
"""
Project: coding_test
File: Extract_region.py
IDE: PyCharm
Creator: morei
Email: zhangmengleiba361@163.com
Create time: 2021-05-26 15:36
Introduction:
"""
import time
import argparse
import os
import tqdm
import pandas as pd
from pathlib import Path


def merge_gene_utr(df):
    uni_g = list(set(df['gene_id']))
    merge_l = []
    for g in uni_g:
        temp = df[df['gene_id'] == g].reset_index(drop=True)
        seq, start, end, gene, tran = temp['seq_id'][0], temp['start'][0], temp['end'][0], \
                                      temp['gene_id'][0], temp['transcript_id'][0]
        if len(temp) > 1:
            for i in range(len(temp)):
                sq, s, e, g, t = temp['seq_id'][i], temp['start'][i], temp['end'][i], \
                                      temp['gene_id'][i], temp['transcript_id'][i]
                if (s >= start & s <= end) | (start >= s & start <= e):
                    seq, start, end, gene, tran = sq, min([start, s]), min([end, e]), g, 'merged'
                else:
                    merge_l.append([sq, s, e, g, t])
            merge_l.append([seq, start, end, gene, tran])
        else:
            merge_l.append([seq, start, end, gene, tran])
    df = pd.DataFrame(merge_l, columns=['seq_id', 'start', 'end', 'gene_id', 'transcript_id'])
    df = df.sort_values(by=['seq_id', 'start'], axis=0, ascending=[True,True])
    return df


class ExtractRegion:
    def __init__(self, file):
        self.file = file

    def get_region(self):
        out_dir = 'output'
        if not Path(out_dir).exists():
            os.makedirs(out_dir)
        cmd = f"cat {self.file} | cut -f1,3,4,5,9 | cut -f1,2 -d ';' " \
              "| awk '{print $1, $2, $3, $4, $6, $8}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' " \
              f"| sed -e 's/\;//g' | sed '1,5d' > {out_dir}/temp.noheader"
        os.system(cmd)
        info = pd.read_table(f'{out_dir}/temp.noheader', header=None)
        info.columns = ['seq_id', 'type', 'start', 'end', 'gene_id', 'attribute']
        # gene, isoforms, coding exons
        gene = info[info['type'] == 'gene']
        isoforms = info[info['type'] == 'transcript']
        isoforms['isoform_id'] = isoforms['gene_id'] + '_' + isoforms['attribute']
        exon_cds = info[info['type'].isin(['exon', 'CDS'])]
        # UTR
        utr = info[info['type'].isin(['UTR'])]
        utr_trans = utr.merge(isoforms, on=['seq_id', 'gene_id', 'attribute'], how='left')
        utr_trans['start_d'] = utr_trans['start_x'] - utr_trans['start_y']  # utr and transcript start point distance
        utr_trans['end_d'] = utr_trans['end_y'] - utr_trans['end_x']  # utr and transcript end point distance
        # 3'UTR isoform level
        utr3 = utr_trans[utr_trans['start_d'] < utr_trans['end_d']]  # closer to transcript start point is 3'UTR
        utr3['type'] = "3'UTR"
        feat = ['seq_id', 'start_x', 'end_x', 'gene_id', 'attribute']
        utr3 = utr3[feat].rename(columns={'start_x': 'start', 'end_x': 'end', 'attribute': 'transcript_id'})
        # 5'UTR isoform level
        utr5 = utr_trans[utr_trans['start_d'] > utr_trans['end_d']]  # closer to transcript start point is 5'UTR
        utr5['type'] = "5'UTR"
        utr5 = utr5[feat].rename(columns={'start_x': 'start', 'end_x': 'end', 'attribute': 'transcript_id'})
        # UTR gene level
        utr3_gene = merge_gene_utr(utr3)
        utr5_gene = merge_gene_utr(utr5)
        # output data
        item_df_l = [gene, isoforms, exon_cds, utr3, utr5, utr3_gene, utr5_gene]
        item_l = ['gene', 'isoforms', 'exon_cds', 'utr3', 'utr5', 'utr3_gene', 'utr5_gene']
        for i in range(len(item_l)):
            item_df_l[i].to_csv(f'{out_dir}/{item_l[i]}.bed', sep='\t', index=False, header=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This tool helps to extract region from gtf file')
    parser.add_argument('-f', '--file',
                        type=str,
                        default="",
                        help='gtf file to be extracted'
                        )
    args = vars(parser.parse_args())
    # input
    file = args['file']
    # extract region
    ExtractRegion = ExtractRegion(file)
    startTime = time.time()  # analysis start time
    print('Analysis starts, please wait for a moment')
    ExtractRegion.get_region()  # analysis processing
    executionTime = (time.time() - startTime)  # analysis end time
    print('Execution time in minutes: ' + str(executionTime/60))

