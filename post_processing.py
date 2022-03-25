import pandas as pd
import numpy as np
import seaborn as sns

# df = pd.read_csv('farm_arch_gw_20220321_x200y128_v2.csv')
df = pd.read_csv('farm_arch_gw_20220321_x32y103_v2.csv')
# df = pd.read_csv('farm_arch_gw_20220322_x38y97.csv')

df_crop_subset = df[(df.crop == 'Corn')]

sns.set_style("darkgrid")
# palette = sns.color_palette("mako_r", 7) # mako_r sns.light_palette("seagreen", as_cmap=True)
palette = sns.color_palette("PuBu", 12) # sns.color_palette("OrRd", 7) # sns.color_palette("YlGn", 7) #
# palette.reverse()
sns.lineplot(x="year", y="xs_total", hue="trans", palette=palette, data=df, estimator=np.sum, ci = None)

sns.set_style("darkgrid")
# palette = sns.color_palette("mako_r", 7) # mako_r sns.light_palette("seagreen", as_cmap=True)
palette = sns.color_palette("YlGn", 12) # sns.color_palette("OrRd", 7) # sns.color_palette("YlGn", 7) #
# palette.reverse()
sns.lineplot(x="year", y="xs_total", hue="trans", palette=palette, data=df_crop_subset, estimator=np.sum, ci = None)
