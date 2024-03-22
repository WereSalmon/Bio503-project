import string
import time
import sys, os
import argparse
import numpy
import shutil
import re
import fnmatch
import gzip

try:
	import pandas as pd
	HAS_PD_MODULE = True;
except ImportError as e:
	HAS_PD_MODULE = False;
	name = 'pandas'
	if str(e) != "No module named '%s'"%(name):
		print("WARNING: An error occurred while importing the module '%s'")
	else:
		print("WARNING: The module '%s' was not found.\n")

from numpy import NaN

## define functions

def process_tax_table(filepath, tax_system="SILVA"):
	'''overwrite a phyloseq taxonomy table so that it is in green_genes format'''
	
	## import the taxonomy table from phyloseq
	tax_table = pd.read_csv(filepath, sep = "\t")
	    
	## initialize list of correct taxonomic classifications
	tax_list = []
	
	if tax_system =="SILVA":
		## create the list of correctly formatted taxonomic classifications
		for row in tax_table.iterrows():
			d = "D_0__" + row[1]["domain"]
			p = "D_1__" + row[1]["phylum"]
			o = "D_2__" + row[1]["order"]
			c = "D_3__" + row[1]["class"]
			f = "D_4__" + row[1]["family"]
			g = "D_5__" + row[1]["genus"]
			s = "D_6__" + row[1]["species"]
			tax = [d,p,o,c,f,g,s]
			tax_list.append((row[1]["OTUID"], tax))
    
	elif tax_system =="Greengenes":
		## create the list of correctly formatted taxonomic classifications
		for row in tax_table.iterrows():
			d = "k__" + row[1]["domain"]
			p = "p__" + row[1]["phylum"]
			o = "o__" + row[1]["order"]
			c = "c__" + row[1]["class"]
			f = "f__" + row[1]["family"]
			g = "g__" + row[1]["genus"]
			s = "s__" + row[1]["species"]
			tax = [d,p,o,c,f,g,s]
			tax_list.append((row[1]["OTUID"], tax))
	else:
		raise ValueError('Unknown taxonomy system used. Use either SILVA or Greengenes taxonomy system.')
	
	## save green_tax as a dataframe
	new_tax_table = pd.DataFrame(tax_list)
    
	## rename columsn of the otu table
	new_tax_table = new_tax_table.rename(columns={0:"# OTUID", 1:"taxonomy"})
    
	## save the green_genes taxonomy table as TSV
	new_tax_table.to_csv(filepath, index = None, sep = "\t")

## parser

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Parser")

	# parse command line arguments
	parser.add_argument('-i','--input_table', action='store', required=True, help="Path to a classical table file (TSV, CSV or similar).");
	parser.add_argument('-v','--verbose', action='store_true', dest="verbose", default=False, help='Show lots of information.');
	parser.add_argument('--tax_system', action='store', dest="tax_system", default="SILVA", choices = ["SILVA", "Greengenes"], help='Tell python which taxonomy system you are using.');

	
	args = parser.parse_args()
	
	# basic file checking
	if(not os.path.isfile(args.input_table)):
		print("ERROR: Input table '%s' does not exist or is not a file" % (args.input_table))
		sys.exit(1)
	if(args.verbose): print("Reading input table..")
	process_tax_table(args.input_table, args.tax_system)
	if(args.verbose): print("Reformatting taxonomy information into '%s' format"%(args.tax_system))
	sys.exit(0)