import pandas as pd
import numpy as np
import re, json, itertools, os

gene_list = ["ABCG2", "CACNA1S", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "CFTR", "DPYD", "G6PD", "HLA-A", "HLA-B", "IFNL3", "MT-RNR1", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1", "CYP3A4", "CYP2C8"]

# gene_list = [#"G6PD", "MT-RNR1", "ABCG2", "CACNA1S", "CFTR", "IFNL3", "VKORC1", "RYR1",
#              "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6",
#              "CYP3A4", "CYP3A5", "CYP4F2", "DPYD", "NUDT15",
#              "SLCO1B1", "TPMT", "UGT1A1"]

# Degenerate base means more than one base possibility at a particular position,
# this is usually the case when a DNA sequence is derived from amino acid sequence with codon based sequence.
dic_degenerate_bases = {
  'R': ['A', 'G'],
  'Y': ['C', 'T'],
  'M': ['A', 'C'],
  'K': ['G', 'T'],
  'S': ['G', 'C'],
  'W': ['A', 'T'],
  'H': ['A', 'T', 'C'],
  'B': ['G', 'T', 'C'],
  'V': ['G', 'A', 'C'],
  'D': ['G', 'A', 'T'],
  'N': ['A', 'T', 'C', 'G']
}

######################
##### definition #####
######################
output = {}
for gene in gene_list:
  print(gene)
  allele_definition_table = "./data/pgx_diplotypes/definition/%s_allele_definition_table.txt" % gene
  if os.path.isfile(allele_definition_table):
    define_df = pd.read_csv(allele_definition_table, sep='\t', index_col=0)
  else:
    continue
  # Replace np.nan with epsilon
  define_df = define_df.replace(np.nan, '')
  define_df_raw = define_df.copy()
  single_gene = {}
  single_gene['reference_genome'] = 'GRCh38'
  chrom = re.findall('chromosome (\w+),', define_df.index[1])[0]
  single_gene['chrom'] = chrom
  chrom_nc = 'chr' + chrom + ':' + re.findall('NC\w+.\d+', define_df.index[1])[0]
  pc = define_df.iloc[0,]
  ng = define_df.iloc[1,]
  rs = define_df.iloc[3,].to_list()
  reference = define_df.iloc[5,].to_list()
  
  # Display positions and alleles
  display = []
  for i in range(0, len(ng)):
    display.append(chrom_nc + ':' + ng[i] + ':' + rs[i] + ':' + pc[i])
  
  # Transform the allele definition to help cpat matching
  pos_rs = []
  for i in range(0, len(ng)):
    ng_item = ng[i]
    matchobj = re.search(r'\w\.(\d+)\_(\d+)(del|ins)(\w*)', ng_item)
    if matchobj:
      pos_start = int(matchobj.group(1)) - 1
      pos_end = int(matchobj.group(2))
      if re.search(r'del|ins', ng_item):
        pos_start = int(matchobj.group(1)) - len(reference[i])
      pos = '%s-%s' % (pos_start, pos_end)
    else:
      matchobj = re.search(r'\w\.(\d+)(\w*)', ng_item)
      if matchobj:
        pos = matchobj.group(1)
        if re.findall('del|ins|dup', ng_item):
          tmp = define_df[define_df.columns[i]].iloc[5:,].to_list()
          max_len = len(re.sub('del|ins', '', tmp[0]))
          min_len = len(re.sub('del|ins', '', tmp[0]))
          brackets = 0
          for tt in tmp:
            if "(" in tt:
              brackets = 1
            for t in tt.split('; '):
              if len(re.sub('del|ins', '', t)) > max_len:
                max_len = len(re.sub('del|ins', '', t))
              if len(re.sub('del|ins', '', t)) < min_len and t != '':
                min_len = len(re.sub('del|ins', '', t))
          ### Dup has some unique features
          if 'dup' in ng_item:
            # Change expression in this section
            pos = '%s-%s' % (int(pos) - min_len, pos)
            ref_dup_base = tmp[0]; ref_dup_count = len(ref_dup_base)
            tmp[0] = 'ref%s' % ref_dup_base
            if brackets == 0:
              for z in range(1, len(tmp)):
                if tmp[z] != '':
                  dup_count = int(len(tmp[z])/ref_dup_count)
                  tmp[z] = 'ins%s' % (ref_dup_base * dup_count)
                else:
                  tmp[z] = 'ref%s' % ref_dup_base
            else:
              for z in range(1, len(tmp)):
                if tmp[z] != '':
                  print(tmp[z])
                  dup_count = int(re.findall('\((\d+)\)', tmp[z])[0])
                  dup_base = tmp[z].split('(')[0]
                  tmp[z] = 'ins%s' % (ref_dup_base * dup_count)
                  print(tmp[z])
                else:
                  tmp[z] = 'ref%s' % ref_dup_base
            define_df[define_df.columns[i]].iloc[5:,] = tmp
            print(ng_item, define_df[define_df.columns[i]].iloc[5:,].to_list(), max_len, pos)
          else:
            pos = '%s-%s' % (int(pos) - max_len - 1, pos)
      else:
        print('Warning!')
        print(ng_item)
    ## Change position
    if re.findall('\((\d+)\)', reference[i]):
      ref_dup_base = reference[i].split('(')[0]
      pos = '%s-%s' % (int(pos)-len(ref_dup_base), pos)
    pos_rs.append('%s:%s' % (str(pos), rs[i]))
  
  ##### Definition of Haplotype
  haplotype_definition = {}
  haplotype_mutated_loci = {}
  single_gene['reference_haplotype'] = define_df.index[5]
  # Because there are some change of reference in the above section
  reference = define_df.iloc[5,].to_list()
  for index, row in define_df.iloc[5:,].iterrows():
    hap_name = index
    defined_base = []; defined_loci = []
    for i in range(0, len(row)):
      # Defined loci
      if row[i] != '' and row[i] != reference[i]:
        defined_loci.append(pos_rs[i])
      elif index == define_df.index[5]:
        defined_loci.append(pos_rs[i])
      # Defined bases
      if '-' in pos_rs[i] and re.findall('\d', reference[i]) != []:
        ### This one is the reference
        # ref haplotype
        ref_dup_base = reference[i].split('(')[0]
        if row[i] == reference[i] or row[i] == '':
          base = 'ref%s' % ref_dup_base
        else:
          if re.findall('\((\d+)\)', reference[i]):
            ref_dup_count = int(re.findall('\((\d+)\)', reference[i])[0])
            # other haplotype
            dup_count = int(re.findall('\((\d+)\)', row[i])[0])
            dup_base = row[i].split('(')[0]
            base = [ref_dup_base] * abs(dup_count - ref_dup_count)
            if dup_count - ref_dup_count <0:
              base = 'del%s' % "".join(base)
            else:
              base = 'ins%s' % "".join(base)
          else:
            print(row[i], reference[i])
        defined_base.append([base])
      else:
        if row[i] == '':
          tmp = [reference[i]]
        else:
          tmp = row[i].split('; ')
        # Degenerate bases
        tmp_res = []
        for t in tmp:
          if t in list(dic_degenerate_bases.keys()):
            tmp_res.extend(dic_degenerate_bases[t])
          else:
            tmp_res.append(t)
        defined_base.append(tmp_res)
      #print(defined_base)
    
    haplotype_mutated_loci[index] = defined_loci
    haplotype_definition[index] = dict(zip(pos_rs, defined_base))
  
  
  # Display alleles
  haplotype_definition_display = {}
  reference_raw = define_df_raw.iloc[5,].to_list()
  for index, row in define_df_raw.iloc[5:,].iterrows():
    defined_base_display = []
    for i in range(0, len(row)):
      if row[i] == '':
        defined_base_display.append(display[i] + ':' + reference_raw[i])
      else:
        defined_base_display.append(display[i] + ':' + row[i])
    haplotype_definition_display[index] = dict(zip(pos_rs, defined_base_display))
  
  single_gene['haplotype_mutated_loci'] = haplotype_mutated_loci
  single_gene['haplotype_definition'] = haplotype_definition
  single_gene['haplotype_definition_display'] = haplotype_definition_display
  
  ##### Frequency of Diplotypes
  diplotype_frequency = {}
  
  # Generate diplotypes based on Hardy–Weinberg equilibrium
  haps = list(haplotype_definition.keys())
  if gene == 'CYP2C19':
    haps = ['*1', '*2', '*3', '*4', '*5', '*6', '*7', '*8', '*9', '*10', '*11', '*12', '*13', '*14', '*15', '*16', '*17', '*18', '*19', '*22', '*23', '*24', '*25', '*26', '*28', '*29', '*30', '*31', '*32', '*33', '*34', '*35', '*38', '*39']
  dips = list(itertools.combinations_with_replacement(haps, 2))
  
  # Frequency of haplotype
  races = ["African American/Afro-Caribbean", "American", "Central/South Asian", "East Asian", "European", "Latino", "Near Eastern", "Oceanian", "Sub-Saharan African"]
  hap_freq_fp = "./data/pgx_diplotypes/population_frequency/%s_frequency_table.txt" % gene
  hap_freq_df = pd.read_csv(hap_freq_fp, sep="\t", index_col=0)
  # If the populations is not complete, insert empty column
  for race in races:
    if race not in hap_freq_df.columns:
      hap_freq_df[race] = np.nan
  # Replace np.nan with epsilon
  hap_freq_df = hap_freq_df[races].replace(np.nan, 1e-5)
  
  # Frequency of diplotype
  for dip in dips:
    diplotype = dip[0] + '/' + dip[1]
    hap1_freq = hap_freq_df.loc[dip[0],]
    hap2_freq = hap_freq_df.loc[dip[1],]
    if dip[0] == dip[1]: # homozygotes
      dip_freq = hap1_freq * hap2_freq
    else: # heterozygotes
      dip_freq = 2 * hap1_freq * hap2_freq
    diplotype_frequency[diplotype] = dip_freq.to_dict()
  
  single_gene['diplotype_frequency'] = diplotype_frequency
  
  output[gene] = single_gene



with open("./data/panno/pgx_diplotypes.json", 'w', encoding='utf-8') as f:
  json.dump(output, f, ensure_ascii=False, indent=2)

