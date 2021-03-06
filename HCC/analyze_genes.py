import pandas
import numpy
import sqlite3 as lite
from pandas.io import sql
import os

# Path directories for this database
working_directory = "/Users/ehrlich/Documents/HCC/"
db_prefix = str.split(working_directory, "/")[-2]
tcgi_prefix = 'TCGI'

# Open connect to database
con = lite.connect(working_directory + db_prefix + '_geneSequenceResults.db')
cur = con.cursor() 


# Create isoform table to which we will write in the MAIN LOOP below
cur.execute("DROP TABLE IF EXISTS genes");
cur.execute("CREATE TABLE genes (geneid text, normcount real, barcode text, tissuetype text)")




# Read in the file_manifest.txt
file_manifest = pandas.read_csv(working_directory + tcgi_prefix + '_' + db_prefix + '/file_manifest.txt', sep='\t').shift(-2)
# Only need one entry for each sample....can find all files based on prefix of 'File Name' 
subset_manifest = file_manifest.drop_duplicates(['Sample','Barcode']).dropna().reset_index(drop=True)


## MAIN LOOP
# Loops through the samples
for i in subset_manifest.index:
    # Get the i^th sample and its barcode
    sample, barcode = subset_manifest.loc[i,['Sample','Barcode']]

    # Specify the directory path, and the prefix that is unique to the 6 files associated with this sample
    path1 = working_directory + tcgi_prefix + '_' + db_prefix + '/'
    path2 = subset_manifest.loc[i,'Platform Type'] + "/" + subset_manifest.loc[i,'Center'] + '__' + subset_manifest.loc[i,'Platform'] + "/Level_" + subset_manifest.loc[i,'Level'] + "/"
    path = path1 + path2
    file_prefix = '.'.join( str.split(subset_manifest.loc[i,'File Name'], '.')[:3] )

    # Read in sample's relevant file on isoform.normalized_results
    file = os.popen('ls ' + path + file_prefix + '*genes.normalized_results*').read()[:-1]
    genes_results = pandas.read_csv(file, '\t')

    # Add the barcode as a column to the dataframe
    genes_results['barcode'] = barcode
    # Add tissue_type as a column to the dataframe
    if barcode[13]=='0':
        tissue_type = 'TN'
    else:
        tissue_type = 'NT'
    genes_results['tissuetype'] = tissue_type


    # Append the dataframe tothe isoform table in the geneSequenceResults database
    genes_results.rename(columns={'gene_id':'geneid','normalized_count':'normcount'}, inplace=True)
    sql.to_sql(genes_results, name='genes', con=con, flavor='sqlite', if_exists='append', index=False)
    print("Completed sample: ",i)
    ## END OF LOOP

# Commit to the changes and close connection to geneSequenceResults.db
con.commit()
con.close()



#####################################
# END OF FILE
#####################################
