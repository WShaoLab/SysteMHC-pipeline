#!/usr/bin/env python3
# -*- coding: utf-8  -*-


import sys
from msproteomicstoolslib.format import pepXMLReader
import csv

infile = sys.argv[1]
outfile = sys.argv[2]
r = pepXMLReader.pepXMLReader(infile)

writer = csv.writer(open(outfile, 'w'))
writer.writerow(
    ["Spectrum", "Scan_number", "Peptide_Sequence", "Modified_Pep_Seq","Ions", "Modification_Position", "Massshift", "Charge", "RT", "PrecursorMass", "Massdiff", "Expect_value", "Fvalue", "Pvalue", "MatchedIons", "TotalIons",  "Protein", "Probability"])
for hit in r.parse_all():
  #exp = float(hit.expect)
  #pval = float(hit.pvalue)
  #exp = str(hit.expect).replace('e', 'E')
  position = "NA"
  mshfit = "NA"
  ions = str(hit.modified_peptide) + '/' + str(hit.spectrum_query.charge)
  if len(hit.modifications)>0 :
    position = hit.modifications[0][0]
    mshfit = hit.modifications[0][1]
  pval = "NA"
  if hasattr(hit, "pvalue"): pval = "%.20f"  % hit.pvalue
  fval = "NA"
  if hasattr(hit, "fval"): fval = "%.20f"  % hit.fval
  exp = "NA"
  if hasattr(hit, "expect"): exp = "%.20f"  % hit.expect
  writer.writerow([hit.spectrum_query.spectrum, hit.scan_number, hit.peptide, hit.modified_peptide, ions, position, mshfit, hit.spectrum_query.charge, hit.spectrum_query.retTime, hit.spectrum_query.precursorMass, hit.massdiff, exp, fval, pval, hit.matched_ions, hit.total_ions, hit.protein, hit.probability])
  
