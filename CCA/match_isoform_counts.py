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

# Open connection to the geneSequenceResults.db
con = lite.connect('geneSequenceResults.db');
cur = con.cursor()

# Read the TN (tumor) data in as one dataframe, and read the NT (benign) data in as another dataframe 
TN_table = pandas.read_sql_query("select * from isoforms where tissuetype like 'TN'", con)
NT_table = pandas.read_sql_query("select * from isoforms where tissuetype like 'NT'", con)
# print("TN: ", TN_table.shape)
# print("NT: ", NT_table.shape)

# Create a column with the sample so the TN data and NT data rows can be matched
NT_table['sample'] = NT_table['barcode'].apply(lambda x: x[:12])
TN_table['sample'] = TN_table['barcode'].apply(lambda x: x[:12])


# # Let's check that sample,geneid,isoformid uniquely identify a row
s,g,i = NT_table[['sample','geneid','isoformid']].iloc[0]
NT_check = NT_table[((NT_table['sample']==s) & (NT_table['geneid']==g) & (NT_table['isoformid']==i))]
TN_check = TN_table[((TN_table['sample']==s) & (TN_table['geneid']==g) & (TN_table['isoformid']==i))]
# print("We want to check that a sample,geneid,isoformid combination will UNIQUELY identify 1 row in each table...")
# print("...but we get two rows back from each table ?!?!")
# print()
# print("NT_table: ")
# print(NT_check)
# print()
# print("TN_table: ")
# print(TN_check)
 




# Merge the two dataframes where they have the same isoformid, geneid, and sample
pairs = pandas.merge(NT_table, TN_table, on=['geneid','isoformid','sample'], how='inner', suffixes=('_NT', '_TN'))
print("pairs: ", pairs.shape)

# Save the relevant columns to a comma-separated .csv file  
pairs = pairs[['isoformid','sample','geneid','normcount_TN','normcount_NT']]
pairs.to_csv('isoform_pairs_table.csv', sep=',', index=False)

# Save the pairs dataframe into the geneSequenceResults DataBase as a new table
cur.execute('DROP TABLE IF EXISTS isoform_pairs')
sql.to_sql(pairs, name='isoform_pairs', con=con)

# Commit to the changes and close connection to geneSequenceResults.db
con.commit()
con.close()

#####################################
# END OF FILE
#####################################
