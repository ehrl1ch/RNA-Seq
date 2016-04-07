#####################################
######  	 ABOUT THIS SCRIPT 	    #
#####################################
# This script divides the geneSequenceResults DataBase into NT rows and TN rows
# and then merges matching sample/gene pairs 
# so the normalized counts of TN and NT pairs can be compared side-by-side


#####################################
######  INSTRUCTIONS FOR USE:		#
#####################################
#In the shell go to this directory and type 'python match_gene_counts.py' to run the script.


#####################################
# Import necessary libraries
import pandas
from pandas.io import sql
import sqlite3 as lite


# Path directories for this database
working_directory = "/Users/ehrlich/Documents/HCC/"
db_prefix = str.split(working_directory, "/")[-2]
tcgi_prefix = 'TCGI'

# Open connection to the geneSequenceResults.db
con = lite.connect(working_directory + db_prefix + '_geneSequenceResults.db')
cur = con.cursor()

# Read the TN (tumor) data in as one dataframe, and read the NT (benign) data in as another dataframe 
TN_table = pandas.read_sql_query("select * from genes where tissuetype like 'TN'", con)
NT_table = pandas.read_sql_query("select * from genes where tissuetype like 'NT'", con)

# Create a column with the sample so the TN data and NT data rows can be matched
NT_table['sample'] = NT_table['barcode'].apply(lambda x: x[:12])
TN_table['sample'] = TN_table['barcode'].apply(lambda x: x[:12])

# Merge the two dataframes where they have the same geneid and sample
pairs = pandas.merge(NT_table, TN_table, on=['geneid','sample'], how='inner', suffixes=('_NT', '_TN'))
pairs = pairs[['geneid','sample','normcount_TN','normcount_NT']]

# Save into the geneSequenceResults DataBase as a new table
cur.execute('DROP TABLE IF EXISTS gene_pairs')
sql.to_sql(pairs, name='gene_pairs', con=con)

# Close connection to geneSequenceResults.db
con.close()

#####################################
# END OF FILE
#####################################
