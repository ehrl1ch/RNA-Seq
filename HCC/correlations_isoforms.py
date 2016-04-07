import sys
import pandas
import numpy
import sqlite3 as lite
from pandas.io import sql
import os
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import t

# Read in command line arguments
# eg. MEN1 is our gene_of_interest, looking for other genes that are correlated at the alpha=0.25 significance-level
if len(sys.argv)>1:
    gene_of_interest = sys.argv[1]
    alpha = float(sys.argv[2])
else:
    gene_of_interest = 'MEN1'
    alpha = 0.25


# Open connect to database
con = lite.connect('geneSequenceResults.db');
cur = con.cursor() 

# Get the TN/NT counts of MEN1's 9 samples x 9 isoforms 
men1_isoforms = pandas.read_sql_query("select * from isoform_pairs where geneid like '"+gene_of_interest+"%'", con).drop("index", axis=1)
men1_isoforms["diff_normcount"] = men1_isoforms['normcount_TN'].astype(float) - men1_isoforms['normcount_NT'].astype(float)
men1_diff = men1_isoforms.pivot(index='sample', columns='isoformid', values='diff_normcount')

# Get the Gene IDs and TN/NT counts of all other samples
gene_pairs = pandas.read_sql("select * from gene_pairs", con=con)
gene_pairs["diff_normcount"] = gene_pairs['normcount_TN'].astype(float) - gene_pairs['normcount_NT'].astype(float)
gene_diff = gene_pairs.pivot(index='sample', columns='geneid', values='diff_normcount')


# Calculate the pair-wise correlations between the 20532 genes and 9 isoforms
r = pandas.DataFrame(columns=men1_diff.columns, index=gene_diff.columns)
for i in men1_diff.columns:
    r[i] = gene_diff.corrwith( men1_diff[i] )

# Calculate the t-values and p-values for a 2-tailed significance
sample_size = men1_diff.shape[0] 
degrees_freedom = sample_size-2
standard_error = numpy.sqrt( (1-r*r)/degrees_freedom )
t_values = r / standard_error
p_values = t_values.applymap( lambda x: t.sf(abs(x),degrees_freedom)*2 )


# Holm-Bonferroni adjustment
alpha_string = "%1.2g-significant" %alpha
p_unstacked = pandas.DataFrame( p_values.dropna(how='all',axis=1).dropna(how='all').unstack(level=0), columns=['p-values'])
significance_holm_bonferroni = multipletests( p_unstacked['p-values'], alpha=alpha, method='holm' )
alpha_holm_bonferroni = significance_holm_bonferroni[3]
significance = pandas.DataFrame(significance_holm_bonferroni[1], columns=["p-adjusted"]).join( pandas.Series(significance_holm_bonferroni[0], name=alpha_string) ).set_index(p_unstacked.index)
isoforms_significant_men1 = p_unstacked.join( significance )

# Select geneids which are significant at the alpha-level
list_of_significant_genes = isoforms_significant_men1[isoforms_significant_men1[alpha_string]==True]
print("\nNumber of significantly correlated gene-isoform combos: %i\n" % len(list_of_significant_genes) )
print(list_of_significant_genes)


# Save the tables to the geneSequenceResults DataBase
cur.execute('DROP TABLE IF EXISTS genes_significant_with_men1_isoforms')
sql.to_sql(isoforms_significant_men1, name='genes_significant_with_men1_isoforms', con=con)


# Commit to the changes and close connection to geneSequenceResults.db
con.commit()
con.close()






#####################################
# END OF FILE
#####################################

