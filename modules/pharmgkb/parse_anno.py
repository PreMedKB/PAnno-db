import pandas as pd
import numpy as np
import re


def process_uncomplete(uncomplete):
  ###### Uncomplete
  still_manual = pd.DataFrame()
  uncomplete = uncomplete.reset_index(drop=True)
  unknown_regL = ["^There is (?:currently |current )?no (?:available |published )?(?:evidence|studies|study|information|association)", \
    "^There were too few patients .* for a conclusion to be drawn.", \
    "^No patients with .* (?:were )?(?:available|included|seen|present|studied)", \
    "^No information (?:is )?available", "^No information were reported regarding patients with the .. genotype", \
    "(?:were|was) not (?:reported|analyzed) in (?:the|this) (?:studies|study)", "not (?:been |be )?(?:studied|evaluated)", \
    "^Patients with .* have an unknown (?:response|likelihood)", "Expression study did not included heterozygous genotype"]
  unknown_reg = r"|".join(unknown_regL)
  for index, row in uncomplete.iterrows():
    function = np.nan
    text = row['Annotation Text']
    ######### Starting the judging #########
    ######### Common functions which are not distinguishable in categories #########
    # Normal function
    if re.search("normal function (?!allele)|normal ability", text, re.IGNORECASE):
      function = "Normal"
    # # No function
    # elif re.search("no function (?:allele )?by CPIC", text, re.IGNORECASE):
    #   function = "No"
    # Unrelated
    elif re.search("^The current evidence base suggests (?:that there|there that) is (?:no|not) (?:significant )?association|^No significant association|^This genotype has no significant effect on|Current literature evidence finds no significant effect of|did not have a statistically significant different", text, re.IGNORECASE):
      function = "Unrelated"
    # Uncertain function
    elif re.search("The association between the .. genotype and .* is currently unclear|^It is currently unclear|It is not clear*|Patients with .* may (?:not )?respond to(?!.*increased|decreased)|undocumented response|^No significant (?:results|findings)|uncertain significance", text, re.IGNORECASE):
      function = "Uncertain"
    # Unknown function
    elif re.search(unknown_reg, text):
      function = "Unknown"
    
    ######### Specific functions #########
    ### Dosage
    elif row["Phenotype Category"] == "Dosage":
      tmp_regL = ["increase.*genotype(,)? (but |and )?(a )?decrease", "decrease.*genotype(,)? (but |and )?(an )?increase", "average reponse", "less likely to require a reduction in dose", "a greater chance of achieving target concentrations", "a 300 mg/day dose", "clearance and exposure to doxorubicin may be similar", "increase exposure\w*genotype.*decreased exposure\w*genotype", "decreased exposure\w*genotype.*increase exposure\w*genotype"]
      if re.search(r"|".join(tmp_regL), text, re.IGNORECASE):
        function = "Moderate"
      elif re.search("(?:decrease|reduce|lower|lowest).*(?:dose|dosage|consumption|intake|clearance|requirement|inactivation of thiopurines)|dose reduction|fewer number|more likely to require a reduction in dose|decreased metabolism|increased plasma concentrations|consume less|smoke less|(?:increase|higher|larger|highest).*(?:exposure)|increased risk of leukopenia and neutropenia", text, re.IGNORECASE):
        function = "Decreased"
      elif re.search("(?:increase|higher|larger|highest).*(?:dose|dosage|consumption|intake|clearance|requirement|inactivation of thiopurines)|greater number|increased metabolism|(?:decreased|reduced) plasma concentrations|consume more|a lower chance of achieving target concentrations|a dose equivalent other than 300 mg/day|smoke more|(?:decrease|reduce|lower|lowest).*(?:exposure)|decreased risk of leukopenia and neutropenia", text, re.IGNORECASE):
        function = "Increased"
    
    ### Efficacy
    elif row["Phenotype Category"] == "Efficacy":
      tmp_regL = ["increase.*genotype(,)? (but |and )?(a )?decrease", "decrease.*genotype(,)? (but |and )?(an )?increase", \
        "increased likelihood.*as compared to.*genotype.*decreased likelihood.*genotype", "decreased likelihood.*as compared to.*genotype.*increased likelihood.*genotype", "average reponse", \
        "average risk", "respond similarily to treatment with", "similar response", \
        "increased (?:response|risk|clearance|improvement|treatment|platelet|severity|progression-free|survival|chance of positive treatment|pravastatin plasma).*genotype.*decreased (?:response|risk|clearance|improvement|treatment|platelet|severity|progression-free|survival|chance of positive treatment|pravastatin plasma).*genotype", \
        "decreased (?:response|risk|clearance|improvement|treatment|platelet|severity|progression-free|survival|chance of positive treatment|pravastatin plasma).*genotype.*increased (?:response|risk|clearance|improvement|treatment|platelet|severity|progression-free|survival|chance of positive treatment|pravastatin plasma).*genotype"]
      if re.search(r"|".join(tmp_regL), text, re.IGNORECASE):
        function = "Moderate"
      elif re.search("(?:increase|inceased|better|improved|faster|greater|longer|more|larger|higher|good).*(?:response|improvement|overall survival|reduction|abstinence|duration|progression-free survival)?|(more likely to.*response)|(less likely to suffer)|(longer recurrence-free survival times)|(less likely to drop out of treatment)|(may receive benefit)|(may benefit from treatment)|improvement in left ventricular ejection fraction|have an early response|low on-treatment platelet reactivity", text):
        function = "Increased"
      elif re.search("(?:decrease|deceased|poorer|worse|reduced|slower|lower|shorter|lesser|smaller|less|weaker|unfavorable|diminished).*(?:response|improvement|overall survival|reduction|abstinence|duration|progression-free survival)?|less likely to.*response|more likely to suffer|shorter recurrence-free survival times|more likely to drop out of treatment|may not receive benefit|may experience a limited benefit|inhibition of platelet aggregation|higher frequency of asthma exacerbations|lack of early response|increasing creatinine levels|not respond as well|high on-treatment platelet reactivity|prolonged time to progression", text, re.IGNORECASE):
        function = "Decreased"
      elif re.search("a varying degree of G6PD deficient red blood cells|may have response|do not seem to have different response|were not found to have different response|may not respond", text, re.IGNORECASE):
        function = "Uncertain"
    
    ### Toxicity
    elif row["Phenotype Category"] == "Toxicity":
      moderate_regL = ["increase.*genotype(,)? (but |and )?(a )?decrease", "decrease.*genotype(,)? (but |and )?(an )?increase", \
        "increased likelihood.*as compared to.*genotype.*decreased likelihood.*genotype", "decreased likelihood.*as compared to.*genotype.*increased likelihood.*genotype", "average reponse", "have a similar risk", \
        "increased (?:risk|intraocular pressure|weight gain|general side-effect burden|sexual side-effects|chance).*(genotype)?.*decreased (?:risk|intraocular pressure|weight gain|general side-effect burden|sexual side-effects|chance|as compared to).*(genotype)?", \
        "decreased (?:risk|intraocular pressure|weight gain|general side-effect burden|sexual side-effects|chance).*(genotype)?.*increased (?:risk|intraocular pressure|weight gain|general side-effect burden|sexual side-effects|chance).*(genotype)?", \
        "increased risk.*as compared to patients with.*or may have (a )?decreased, but not absent, risk.*as compared to", \
        "increased likelihood.*genotype.*decreased likelihood.*genotype", \
        "increased likelihood of addiction as compared to patients.*or a decreased risk as compared to patients.*"]
      high_regL = ["later onset of bortezomib-induced peripheral neuropathy", "greater elevations of LDL", "higher drug toxicity", "greater prolactin response", \
        "(?:increase|increased|greater|higher|longer|larger|increased \(but not absent\))[, \w\-]{1,}(?:adverse|risk|severity|likelihood of developing|anxiety|cocaine cue-reactivity|major cardiovascular events rate|Diabetes Mellitus|diabetes|new-onset diabetes|syndrome|toxic liver disease|toxicity|side(?:-| )effect|statin-related muscle symptoms|hyponatremia|progression-free survival|vomiting|postural hypotension|weight gain|systolic|fever|Osteonecrosis|QT(c)? interval|not non-existent risk|transplant rejection|cravings|smoking addiction|addiction|obesity|arthralgia|bone (mineral )?density|gastrointestinal intolerance|intraocular pressure|cognitive impairment|myopathy|hospital stay|treatment interruptions|tamoxifen-induced decrease|triglycerides|peptic ulcer|low-density lipoprotein cholesterol|resistance|homocysteine|angioedema|antipsychotic-induced weight|Lymphopenia|Mucositis|constipation|leukopenia|respiratory depression|cough|drug reaction|myelosuppression|creatinine clearance|hypertension|headache|depression|fatigue|cardiotoxicity|suicidality|lumbar bone loss|triglyceride levels|blood pressure|sensory neuropathy|prolactin|dose delay|rhabdomyolysis|sexual dysfunction|response to methotrexate|metabolism|suicide|hypersensitivity|adverse events|blood glucose|pulse rate|simvastatin-induced adverse|neutropenia|drug-induced liver injury|fasting glucose levels|narcolepsy|spontaneous abortion|survival of red blood cells|cholesterol levels|tardive dyskinesia|response|myalgia|hyperbilirubinemia).*compared", \
        "greater decrease in blood pressure", "greater increase in (?:total cholesterol|fasting glucose)", \
        "reduced survival", "higher odds", "longer hospital stay", \
        "(?:increase|more likely to|may) (?:develop|experience|drink|become addicted to alcohol|have relapsed|require medication|respond to treatment|have a child with a craniofacial abnormality|be tetrahydrocannabinol|discontinue treatment|require glucarpidase treatment).*compared", \
        "a greater increase in blood glucose", "increased risk of Stevens-Johnson Syndrome", \
        "a shorter period of time before chemotherapy-induced ovarian failure", \
        "less likely to require treatment", "gain more weight", "less likely to.*(?:sexual dysfunction|quit smoking)", "more severe (?:nicotine dependence|anemia|side effect)", "less likely to resume menses", "less skin", "(?:less of a|smaller) reduction in (?:HDL-C|cardiovascular events)", "high frequency of adverse events", "may develop malignant hyperthermia", "shorter (?:survival time|overall and event-free survival time)", "lower body mass index", "(?:greater|higher) (?:incidence|levels) of (?:toxicity|side effect)", "lesser decline", "greater bone mineral loss", "report more adverse events", "higher aspirin-induced decline in forced expiratory volume in 1 s", "decreased lung function|estimated glomerular filtration rate", "increased risk for developing new-onset diabetes"]
      low_regL = ["decreased function by CPIC", "earlier onset of bortezomib-induced peripheral neuropathy", "smaller elevations of LDL", "lower drug toxicity", "lower odds", \
        "(?:decrease|decreased|deceased|dcreased|reduced|lower|less|shorter|smaller|fewer|decreased \(but not absent\))[, \w\-]{,}(?:adverse|risk|severity|likelihood of developing|anxiety|cocaine cue-reactivity|major cardiovascular events rate|Diabetes Mellitus|diabetes|new-onset diabetes|syndrome|toxic liver disease|toxicity|side(?:-| )effect|statin-related muscle symptoms|hyponatremia|progression-free survival|vomiting|postural hypotension|weight gain|systolic|fever|Osteonecrosis|QT(c)? interval|not non-existent risk|transplant rejection|cravings|smoking addiction|addiction|obesity|arthralgia|bone (mineral )?density|gastrointestinal intolerance|intraocular pressure|cognitive impairment|myopathy|hospital stay|treatment interruptions|tamoxifen-induced decrease|triglycerides|peptic ulcer|low-density lipoprotein cholesterol|resistance|homocysteine|angioedema|antipsychotic-induced weight|Lymphopenia|Mucositis|constipation|leukopenia|respiratory depression|cough|drug reaction|myelosuppression|creatinine clearance|hypertension|headache|depression|fatigue|cardiotoxicity|suicidality|lumbar bone loss|triglyceride levels|blood pressure|sensory neuropathy|prolactin|dose delay|rhabdomyolysis|sexual dysfunction|response to methotrexate|metabolism|suicide|hypersensitivity|adverse events|blood glucose|pulse rate|simvastatin-induced adverse|neutropenia|drug-induced liver injury|fasting glucose levels|narcolepsy|spontaneous abortion|survival of red blood cells|cholesterol levels|tardive dyskinesia|response|myalgia|hyperbilirubinemia).*compared", \
        "smaller decrease in blood pressure", "(?:smaller|lower) increase in (?:total cholesterol|fasting glucose)", \
        "(?:decrease|less likely to|may not) (?:develop|experience|drink|become addicted to alcohol|have relapsed|require medication|respond to treatment|have a child with a craniofacial abnormality|be tetrahydrocannabinol|discontinue treatment|be resistant to treatment|require glucarpidase treatment).*compared", \
        "a lower increase in blood glucose", "decreased risk of Stevens-Johnson Syndrome", \
        "a longer period of time before chemotherapy-induced ovarian failure", \
        "more likely to require treatment", "gain (?:smaller|less) weight", "more likely to.*(?:sexual dysfunction|quit smoking)", "less severe (?:nicotine dependence|anemia|side effect)", "more likely to resume menses", "more skin", "(?:greater likelihood of|greater) reduction in (?:HDL-C|cardiovascular events)", "longer (?:survival time|overall and event-free survival time)", "greater body mass index", "more less to develop hyperbilirubinemia", "lower (?:incidence|levels) of (?:toxicity|side effect)", "greater decline", "less bone mineral loss", "lower aspirin-induced decline in forced expiratory volume in 1 s", "increased lung function|estimated glomerular filtration rate", "Cardiotox was not not found"]
      if re.search(r"|".join(moderate_regL), text, re.IGNORECASE):
        function = "Moderate"
      elif re.search(r"|".join(high_regL), text, re.IGNORECASE):
        function = "Increased"
      elif re.search(r"|".join(low_regL), text, re.IGNORECASE):
        function = "Decreased"
      elif re.search("an unknown (?:risk|rate|effect)", text, re.IGNORECASE):
        function = "Unknown"
      elif re.search("may have a risk for peptic ulcer|may be at risk for developing narcolepsy", text, re.IGNORECASE):
        function = "Uncertain"

    ### Metabolism/PK
    elif row["Phenotype Category"] == "Metabolism/PK":
      moderate_regL = ["increase.*(genotype|genotypes|allele)(,)? (but |and |or |or may have )?(a )?decrease", "decrease.*(genotype|genotypes|allele)(,)? (but |and |or |or may have )?(an )?increase", "similar enzyme activity", "average catalytic activity", "clearance and exposure.*may be similar", "similar levels of clearance", "similar exposure"]
      if re.search(r"|".join(moderate_regL), text, re.IGNORECASE):
        function = "Moderate"
      elif re.search("(?:reduced|decrease|decease|lower|worse|poorer).*(?:exposure|concentration|sulfation|formation|glucuronidation|maximal rate|uptake|plasma|steady-state levels of vitamin E)|metabolize.*more rapidly|(?:increase|higher|improved|elevated|more rapid).*(?:clearance|metabolism|metabolite|enzymatic activity|kinase|area under the curve|absorption|bioavailability)", text, re.IGNORECASE):
        function = "Increased"
      elif re.search("(?:increase|higher|improved|elevated|more rapid).*(?:exposure|concentration|sulfation|formation|glucuronidation|maximal rate|uptake|plasma|steady-state levels of vitamin E)|metabolize.*more slowly|impaired catalytic activity|(?:reduced|decrease|decrease|lower|worse|poorer).*(?:clearance|metabolism|metabolite|enzymatic activity|kinase|area under the curve|absorption|bioavailability)", text, re.IGNORECASE):
        function = "Decreased"
      else:
        function = "Unknown"
    
    ### Other
    elif row["Phenotype Category"] == "Other":
      function = "Uncertain"
    
    ### Add the function
    row["Allele Function"] = function
    uncomplete.iloc[index] = row
  
  need_manual = uncomplete[uncomplete["Allele Function"].isna()]
  rule_res = uncomplete[uncomplete["Allele Function"].isna() == False]
  need_manual['Comment'] = 'Manual'
  rule_res['Comment'] = 'Rule'
  return(need_manual, rule_res)

