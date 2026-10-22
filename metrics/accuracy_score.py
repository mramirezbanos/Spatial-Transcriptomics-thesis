import pandas as pd

def calculate_and_rank(df, value_column, index_column='settings', columns_column='core_uid'):
    median_table = df.pivot_table(
        values=value_column,
        index=index_column,
        columns=columns_column,
        aggfunc='median'
    )
    ranked_table = median_table.rank(ascending=False, method='min').astype('Int64')
    return ranked_table

# ranked_tables['Accuracy_Score'] = (ranked_tables['ARI'] + ranked_tables['NMI'] + ranked_tables['FMI'] + ranked_tables['Purity'] + ranked_tables['Homogeneity'] + ranked_tables['Completeness']) / 6