for gene in ['CFTR', 'DPYD', 'G6PD', 'MT-RNR1', 'RYR1']:
  allele_definition_table = "./data/pgx_diplotypes/definition/%s_allele_definition_table.txt" % gene
  define_df = pd.read_csv(allele_definition_table, sep='\t', index_col=0)
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
          print('%s$%s$%s$%s$%s' % (gene, index, pos_rs[i].split(':')[1], reference[i], row[i]))
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
