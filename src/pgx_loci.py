"""
Generate bed
"""
import re, json
import pandas as pd
import numpy as np
import requests, sys

## Diplotype related genes
gene_list = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "MT-RNR1", "NUDT15", "IFNL3", 
             "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]

exist_rsid = []
pos_without_rsid = []
for g in gene_list:
  gene_definition_file = "./data/pgx_diplotypes/definition/%s_allele_definition_table.txt" % (g)
  df = pd.read_csv(gene_definition_file, sep='\t', header=None)
  # Remove the column which represents 'Structural Variation'
  if df.iloc[0, df.shape[1]-1] == 'Structural Variation':
    df = df.iloc[:, 0:(df.shape[1]-2)]
  rsID = df.iloc[4, ].to_list()
  # Saving the existed rsIDs
  exist_rsid.extend(rsID[1:])
  NC = df.iloc[2, ].to_list()
  chr = re.split('\s', NC[0])[6].strip(',')
  if g == 'MT-RNR1':
    chr = 'M'
  for i in range(1, len(NC)):
    match = re.findall('(\d{1,})\w', NC[i].split('; ')[0])
    if len(match) == 1:
      one_line = ['chr%s' % chr, match[0], match[0], rsID[i]]
    else:
      one_line = ['chr%s' % chr, match[0], match[1], rsID[i]]
    pos_without_rsid.append(one_line)



## Load PharmGKB.Variant
all_vars = pd.read_csv('PharmGKB.Variant', sep="\t", header=None)
all_rsid = list(set(all_vars[all_vars[0].str.contains('^rs', regex=True)][0].to_list()))
print(len(all_rsid))
all_rsid.extend(exist_rsid)
all_rsid = list(set(all_rsid))
print(len(all_rsid))


## Ensemble API
def ensembl_api(rsid):
  server = "https://rest.ensembl.org"
  ext = "/variation/human/%s?" % rsid
  r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()
  
  decoded = r.json()
  # decipher the positions
  output = []
  for item in decoded['mappings']:
    bed = re.split(':|\-', item['location'])
    if bed[0].startswith('CHR') is False:
      bed[0] = 'chr%s' % bed[0]
      output.append([bed[0], bed[2], bed[1], rsid])
  
  return(output)


pos_from_rsid = []
fail_rsid = []
for rsid in all_rsid:
  try:
    pos = ensembl_api(rsid)
  except:
    fail_rsid.append(rsid)
  pos_from_rsid.extend(pos)


## Processing failed rsIDs manually. rs1799735 and rs36056065 have been withdrawn.
print(fail_rsid)
pos_from_rsid.append(['chr10', '94947938', '94947938', 'rs2031531005'])
pos_from_rsid.append(['chr22', '42126956', '42126956', 'rs1931013246'])
pos_from_rsid.append(['chr12', '21178680', '21178680', 'rs1940852753'])
pos_from_rsid.append(['chr13', '48037885', '48037885', 'rs1950545307'])
pos_from_rsid.append(['chr10', '94938803', '94938803', 'rs2031308986'])

merged = pos_without_rsid + pos_from_rsid
# Formatting
# final = []
# for item in merged:
#   if item != []:
#     if type(item[0]) is list:
#       for ele in item:
#         final.append(ele)
#     else:
#       final.append(item)
#
# merged = final
# final = []
# for item in merged:
#   if item[0].startswith('chrCHR') is False:
#     final.append(item)

# Convert Start and End locations
final = []
for item in merged:
  if int(item[1]) > int(item[2]):
    final.append([item[0], item[2], item[1], item[3]])
  else:
    final.append([item[0], item[1], item[2], item[3]])

pos_rsid_df = pd.DataFrame(final, columns=['chrom', 'start', 'end', 'rsid']).drop_duplicates()
pos_rsid_df[['start', 'end']] = pos_rsid_df[['start', 'end']].astype(int)
pos_rsid_df.sort_values(['chrom', 'start'], inplace=True)
pos_rsid_df.to_csv('./data/bed/pgx_loci_raw.bed', sep='\t', header=0, index=0)


########## pgx_diplotypes
import pandas as pd
pgx_loci = pd.read_csv('./data/bed/pgx_loci_raw.bed', sep='\t', header=None)
chr = pgx_loci[pgx_loci[0].str.contains('^chr', regex=True)]
chr[1] = chr[1] - 5
chr[2] = chr[2] + 5
hla = pgx_loci[pgx_loci[0].str.contains('^HLA', regex=True)]
pgx_loci_padding_bed = pd.concat([chr, hla])
pgx_loci_padding_bed.to_csv('./data/bed/pgx_loci.padding.bed', sep='\t', index=0)



########## pgx_diplotypes
gene_list = [#"G6PD", "MT-RNR1", "ABCG2", "CACNA1S", "CFTR", "IFNL3", "VKORC1", "RYR1",
             "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",
             "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "NUDT15",
             "SLCO1B1", "TPMT", "UGT1A1"]
pgx_loci = pd.read_csv('./data/bed/pgx_loci_raw.bed', sep='\t', header=None)
pgx_diplotypes = json.loads(open("./data/panno/pgx_diplotypes.json").read())
res_list = []
for gene in gene_list:
  info = pgx_diplotypes[gene]
  chrom = 'chr%s' % info['chrom']
  pos_rs = info['haplotype_definition'][info['reference_haplotype']].keys()
  for ele in pos_rs:
    pos = ele.split(':')[0].split('-')
    if len(pos) == 1:
      pos = pos * 2
    res = [chrom, pos[0], pos[1], '']
    res_list.append(res)

df = pd.DataFrame(res_list)
merged_df = pd.concat([df, pgx_loci])#.iloc[:,0:3]

# hla = merged_df[merged_df[0].str.contains('^HLA', regex=True)]
# not_hla = merged_df[merged_df[0].str.contains('^HLA', regex=True) == False]
# from pybedtools import BedTool
# x = BedTool.from_dataframe(not_hla)
# merged_bed = x.sort().merge()
# merged_bed_df = merged_bed.to_dataframe()
# merged_bed_df.start = merged_bed_df.start - 10
# merged_bed_df.end = merged_bed_df.end + 10

# hla.columns = merged_bed_df.columns
# final_bed = pd.concat([merged_bed_df, hla])
final_bed = merged_df
final_bed.to_csv('./data/panno/pgx_loci.bed', sep='\t', index=0, header=0)