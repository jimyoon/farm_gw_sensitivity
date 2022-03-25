pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

if first:
    df_combined = df
    first = False
else:
    df_combined = pd.concat([df_combined, df])

for index, row in gw_cost_curves.iterrows():
    if cum_gw_sum > row['volume']:
        i = index
    else:
        break