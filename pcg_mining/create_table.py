"""
2023-02-28
Fig 3 supplementary table

combine results of three correlation values into one table
"""

import pandas as pd

def rank_df(total:pd.DataFrame, cols:list, ascending=False):
    """
    rank dataframe (total) by values of predefined columns (cols)
    """
    tups = total[cols].sort_values(cols, ascending=ascending).apply(tuple, 1)
    f, i = pd.factorize(tups)
    factorized = pd.Series(f + 1, tups.index)
    ranked_total  = total.assign(Rank=factorized)
    ranked_total = ranked_total.sort_values('Rank')
    return ranked_total

if __name__ == '__main__':
    k = pd.read_csv('kl.csv')
    p = pd.read_csv('pearson.csv')
    s = pd.read_csv('spearman.csv')

    # merge 3 tables into 1
    total = s.merge(p, how='outer', on=['gene','pcg'])
    total = total.merge(k, how='outer', on=['gene','pcg'])

    # add gene symbols
    names = pd.read_csv('/dellfsqd2/ST_OCEAN/USER/liyao1/scripts/pcg_ap_pattern/input_data/planarian_gene_id_mapping.txt', delimiter='\t')
    # change SMED 'names' into blanks
    names.loc[names['gene_name'].str.startswith('SMED'), 'gene_name'] = ''
    # remove duplicates
    names = names[names.duplicated(['smesg','gene_name'], keep='first')]

    total = names.merge(total, how='right', left_on='smesg', right_on='gene')
    total.drop(['dd_Smed_v6', 'dd_Smed_v4',  'smesg',   'smed'],axis=1,inplace=True)

    # change columns names 
    total.columns = ['name in text', 'potential pcg', 'known pcg', 'spearman\'s', 'pearson\'s', 'kullback-leibler']

    # check if there's duplicates
    dub = total[total.duplicated()]
    print(f'duplicates are {dub.shape[0]}')

    # save final table
    total.to_csv('table.csv', index=False)
    
    # cannot just rank by pearson's, spearman's and kl value
    # bc for pearson's and spearman's value, the higher the better
    # however, for kl divergence value, the lower the better
    # to align the rest two values, first convert kl values into neg
    # last time I checked, no og kl value is neg (2023-02-28)

    # rank gene-pcg pairs for each gene first (gene base)
    total['kullback-leibler'] = total['kullback-leibler']*(-1)
    cols = ['spearman\'s','pearson\'s', 'kullback-leibler']

    dfs = []
    for gene in set(total['potential pcg']):
        gdf = total[total['potential pcg']==gene]
        rgdf = rank_df(gdf,cols)
        dfs.append(rgdf)
    rt = pd.concat(dfs)
    rt['kullback-leibler'] = rt['kullback-leibler']*(-1)
    rt = rt.drop_duplicates(keep='first')
    rt.to_csv('ranked_table_by_gene.csv', index=False)


