import pandas as pd
import sys

fn1 = sys.argv[1]
#fn1 = "Mock-clean_fastuniq.markdup.bed.ciro_same.txt"

diff_df = pd.read_csv(fn1, sep='\t', header=None)
diff_df.columns = diff_df.iloc[0]
diff_df = diff_df.drop(diff_df.index[0])
diff_df = diff_df.reset_index(drop=True)

diff_df['combined'] = diff_df['bin_orf_x'] + '-' + diff_df['start_loc_x'].astype(str) + '-' + diff_df['bin_orf_y'] + '-' + diff_df['start_loc_y'].astype(str)

counts = diff_df['combined'].value_counts()

filtered_combined = counts[counts >= 2].index.tolist()
result_df = diff_df[diff_df['combined'].isin(filtered_combined)]

unique_result_df = result_df.drop_duplicates(subset='combined', keep='first')
removed_duplicates_df = result_df[~result_df.index.isin(unique_result_df.index)]

unique_result_df.to_csv(fn1 + ".unique_combined", sep='\t', index=False)

removed_duplicates_df.to_csv(fn1 + ".removed_duplicates", sep='\t', index=False)

removed_reads = removed_duplicates_df['reads_x']
removed_reads_cleaned = removed_reads.str.split(':').str[0]

removed_reads_cleaned.to_csv(fn1 + ".reads_x_cleaned.txt", index=False, header=False)

