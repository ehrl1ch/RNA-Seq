#####################################
######  	 ABOUT THIS SCRIPT 	    #
#####################################
# This script does a 2-tailed paired t-test on the samples of each isoform


#####################################
######  INSTRUCTIONS FOR USE:		#
#####################################
#In the shell go to this directory and type 'python paired_ttest_isoforms.py' to run the script.


#####################################
# Import necessary libraries
import pandas
from pandas.io import sql
import sqlite3 as lite
from scipy import stats
import numpy


# Open connection to the geneSequenceResults.db
con = lite.connect('geneSequenceResults.db');
cur = con.cursor()
#cur.execute('DROP TABLE IF EXISTS isoform_paired_ttest')
#cur.execute('CREATE TABLE isoform_paired_ttest (isoformid TEXT, tvalue REAL, pvalue REAL)')


# Get list of distinct isoform id's
ids = pandas.read_sql_query("select distinct isoformid from isoform_pairs", con)
count = pandas.read_sql_query("select count(*) from isoform_paired_ttest", con).as_matrix()[0][0]

# Loop through each isoform id
#for id in ids.index:
for id in ids.index[count:]:
    i = ids.loc[id][0] 
    i_samples = pandas.read_sql_query("select * from isoform_pairs where isoformid like '" + i + "'", con) #get all samples of the isoform
    i_samples[['normcount_TN','normcount_NT']] = i_samples[['normcount_TN','normcount_NT']].applymap(lambda x: float(x)) #convert the count from strings to floats
    t_value, p_value = stats.ttest_rel(i_samples['normcount_TN'], i_samples['normcount_NT']) #calculate the t-value and p-value
    tmpList = []
    tmpList.extend([i, t_value, p_value])
    cur.execute('INSERT INTO isoform_paired_ttest VALUES (?,?,?)', tmpList) #write the isoformid, tvalue, pvalue into the database
    con.commit()
    print(str(id) + " of " + str(len(ids)) )



alpha = 0.001
ttest = pandas.read_sql_query("select * from isoform_paired_ttest", con)
significant_isoforms = ttest[ttest['pvalue']<alpha]
significant_positive = significant_isoforms[significant_isoforms['tvalue']>0].sort(columns=['pvalue'], axis=0, ascending=True)
significant_negative = significant_isoforms[significant_isoforms['tvalue']<0].sort(columns=['pvalue'], axis=0, ascending=True)
significant_positive['isoformid'].map(lambda x: str.split(str(x), '|')[0]).to_csv("significant_positive_paired_ttest_isoforms.csv", sep="\t", index=False)
significant_negative['isoformid'].map(lambda x: str.split(str(x), '|')[0]).to_csv("significant_negative_paired_ttest_isoforms.csv", sep="\t", index=False)



# Commit to the changes and close connection to geneSequenceResults.db
con.commit()
con.close()


#####################################
# END OF FILE
#####################################


