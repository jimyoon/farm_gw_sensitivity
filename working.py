pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

if first:
    df_combined = df
    first = False
else:
    df_combined = pd.concat([df_combined, df])

gw_multipliers = []
for farm in cum_gw_sum:
    i = 0
    for index, row in gw_cost_curves.iterrows():
        if farm > row['volume']:
            i = index
        else:
            break
    if gw_cost_curves.loc[i][gw_run] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
        gw_multiplier = gw_cost_curves.loc[i][gw_run]
    else:
        gw_multiplier = 9999999 # multiply cost by large number
    gw_multipliers.append(gw_multiplier)

import csv
import pandas as pd
df = pd.DataFrame(water_constraints_by_farm, index=[0]).T
df.to_csv('temp.csv')



net_prices_gw_subset_update = net_prices_gw_subset.copy()
data_profit_subset = data_profit[['nldas', 'farm_id']].head(53835)
cumul_gw_sum = cumul_gw_sum.merge(data_profit_subset, how='left', on='nldas')
for farm_temp in farm_ids_subset:
    for value in crop_ids_by_farm[farm_temp]:
        net_prices_gw_subset_update[value] *= cumul_gw_sum.iloc[farm_temp]['gw_multiplier']

import pickle

id_nldas = dict(zip(data_profit_short.index, data_profit_short.nldas))

with open('20221101_results_timedepletion.p', 'wb') as handle:
    pickle.dump(results_time_depletion, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Revise to account for removal of "Fodder_Herb category"
crop_ids_by_farm_new = {}
for i in crop_ids_by_farm:
    crop_ids_by_farm_new[i] = crop_ids_by_farm[i][0:10]
crop_ids_by_farm = crop_ids_by_farm_new

with open('20221013_results_intermediate.p', 'wb') as handle:
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)


test = pickle.load(open("20221013_results_intermediate.p", "rb"))


with open('dict.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in alphas_total.items():
       writer.writerow([key, value])

(pd.DataFrame.from_dict(data=net_prices_land, orient='index')
   .to_csv('dict_file_net_prices.csv', header=False))


#################################################

import numpy as np
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import math

sns.set_style('whitegrid', {'axes_linewidth': 0, 'axes.edgecolor': 'white'})


def is_significant(value, confidence_interval, threshold="conf"):
    if threshold == "conf":
        return value - abs(confidence_interval) > 0
    else:
        return value - abs(float(threshold)) > 0


def grouped_radial(SAresults, parameters, radSc=2.0, scaling=1, widthSc=0.5, STthick=1, varNameMult=1.3, colors=None,
                   groups=None, gpNameMult=1.5, threshold="conf"):
    # Derived from https://github.com/calvinwhealton/SensitivityAnalysisPlots
    fig, ax = plt.subplots(1, 1)
    color_map = {}

    # initialize parameters and colors
    if groups is None:

        if colors is None:
            colors = ["k"]

        for i, parameter in enumerate(parameters):
            color_map[parameter] = colors[i % len(colors)]
    else:
        if colors is None:
            colors = sns.color_palette("deep", max(3, len(groups)))

        for i, key in enumerate(groups.keys()):
            # parameters.extend(groups[key])

            for parameter in groups[key]:
                color_map[parameter] = colors[i % len(colors)]

    n = len(parameters)
    angles = radSc * math.pi * np.arange(0, n) / n
    x = radSc * np.cos(angles)
    y = radSc * np.sin(angles)

    # plot second-order indices
    for i, j in itertools.combinations(range(n), 2):
        # key1 = parameters[i]
        # key2 = parameters[j]

        if is_significant(SAresults["S2"][i][j], SAresults["S2_conf"][i][j], threshold):
            angle = math.atan((y[j] - y[i]) / (x[j] - x[i]))

            if y[j] - y[i] < 0:
                angle += math.pi

            line_hw = scaling * (max(0, SAresults["S2"][i][j]) ** widthSc) / 2

            coords = np.empty((4, 2))
            coords[0, 0] = x[i] - line_hw * math.sin(angle)
            coords[1, 0] = x[i] + line_hw * math.sin(angle)
            coords[2, 0] = x[j] + line_hw * math.sin(angle)
            coords[3, 0] = x[j] - line_hw * math.sin(angle)
            coords[0, 1] = y[i] + line_hw * math.cos(angle)
            coords[1, 1] = y[i] - line_hw * math.cos(angle)
            coords[2, 1] = y[j] - line_hw * math.cos(angle)
            coords[3, 1] = y[j] + line_hw * math.cos(angle)

            ax.add_artist(plt.Polygon(coords, color="0.75"))

    # plot total order indices
    for i, key in enumerate(parameters):
        if is_significant(SAresults["ST"][i], SAresults["ST_conf"][i], threshold):
            ax.add_artist(plt.Circle((x[i], y[i]), scaling * (SAresults["ST"][i] ** widthSc) / 2, color='w'))
            ax.add_artist(
                plt.Circle((x[i], y[i]), scaling * (SAresults["ST"][i] ** widthSc) / 2, lw=STthick, color='0.4',
                           fill=False))

    # plot first-order indices
    for i, key in enumerate(parameters):
        if is_significant(SAresults["S1"][i], SAresults["S1_conf"][i], threshold):
            ax.add_artist(plt.Circle((x[i], y[i]), scaling * (SAresults["S1"][i] ** widthSc) / 2, color='0.4'))

    # add labels
    for i, key in enumerate(parameters):
        ax.text(varNameMult * x[i], varNameMult * y[i], key, ha='center', va='center',
                rotation=angles[i] * 360 / (2 * math.pi) - 90,
                color=color_map[key])

    if groups is not None:
        for i, group in enumerate(groups.keys()):
            print(group)
            group_angle = np.mean([angles[j] for j in range(n) if parameters[j] in groups[group]])

            ax.text(gpNameMult * radSc * math.cos(group_angle), gpNameMult * radSc * math.sin(group_angle), group,
                    ha='center', va='center',
                    rotation=group_angle * 360 / (2 * math.pi) - 90,
                    color=colors[i % len(colors)])

    ax.set_facecolor('white')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.axis('equal')
    plt.axis([-2 * radSc, 2 * radSc, -2 * radSc, 2 * radSc])
    # plt.show()

    return fig

# generate Sobol samples
param_samples = saltelli.sample(problem, 1000)

# extract each parameter for input into the lake problem
b_samples = param_samples[:, 0]
q_samples = param_samples[:, 1]
mean_samples = param_samples[:, 2]
stdev_samples = param_samples[:, 3]
delta_samples = param_samples[:, 4]

# run samples through the lake problem using a constant policy of .02 emissions
pollution_limit = np.ones(100) * 0.02

# initialize arrays to store responses
max_P = np.zeros(len(param_samples))
utility = np.zeros(len(param_samples))
inertia = np.zeros(len(param_samples))
reliability = np.zeros(len(param_samples))

# run model across Sobol samples
for i in range(0, len(param_samples)):
    print("Running sample " + str(i) + ' of ' + str(len(param_samples)))
    max_P[i], utility[i], inertia[i], reliability[i] = lake_problem(pollution_limit,
                                                                    b=b_samples[i],
                                                                    q=q_samples[i],
                                                                    mean=mean_samples[i],
                                                                    stdev=stdev_samples[i],
                                                                    delta=delta_samples[i])

# Get sobol indicies for each response
SA_max_P = sobol.analyze(problem, max_P, print_to_console=False)
SA_reliability = sobol.analyze(problem, reliability, print_to_console=True)
SA_inertia = sobol.analyze(problem, inertia, print_to_console=False)
SA_utility = sobol.analyze(problem, utility, print_to_console=False)

# define groups for parameter uncertainties
groups = {"Water Scenarios": ["sw_mult"],
          "Econ Scenarios": ["price_mult"],
          "Model Coefficients": ["gammas_mult","K"]}

fig = grouped_radial(sobol_indices, ['sw_mult','price_mult','gammas_mult','K'], groups=groups, threshold=0)
plt.show()

###########################################################

alphas_mult = 1  # JY temp to disable alphas in sensitivity analysis

# generate gw cost curves if gw_curve_option is "internal"
if gw_curve_option == 'internal':
    gw_cost_curve_internal = Theis_pumping.Analytical(S, m, K, WL, R, IrrDepth, years)

# apply sensitivity multipliers
for key in gammas_total_subset:
    gammas_total_subset[key] *= gammas_mult

for key in nirs_subset:
    nirs_subset[key] *= nir_mult

for key in sw_calib_constr_subset:
    sw_calib_constr_subset[key] *= sw_mult

for key in net_prices_land_subset:
    net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - alphas_total_subset[key]) * price_mult

for key in net_prices_land_subset:
    net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - (alphas_total_subset[key] * alphas_mult))

# initialize counters/trackers
first = True
cumul_gw_sum = 0
gw_multiplier = 1

# start simulation time loop
for t in range(no_of_years):
    print('year' + str(t))
    print('!!JY!! ' + str(nir_mult))

    # initialize pyomo model
    fwm_s = ConcreteModel()
    fwm_s.ids = Set(initialize=ids_subset_sorted)
    if first:
        fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset, mutable=True)
    else:
        net_prices_gw_subset_update = net_prices_gw_subset.copy()
        net_prices_gw_subset_update.update((x, y * gw_multiplier) for x, y in net_prices_gw_subset_update.items())
        fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset_update, mutable=True)
    fwm_s.farm_ids = Set(initialize=farm_id)
    fwm_s.crop_ids_by_farm = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
    fwm_s.crop_ids_by_farm_and_constraint = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
    fwm_s.net_prices_sw = Param(fwm_s.ids, initialize=net_prices_sw_subset, mutable=True)
    fwm_s.net_prices_total = Param(fwm_s.ids, initialize=net_prices_land_subset, mutable=True)
    fwm_s.gammas_total = Param(fwm_s.ids, initialize=gammas_total_subset, mutable=True)
    fwm_s.land_constraints = Param(fwm_s.farm_ids, initialize=max_land_constr_subset, mutable=True)
    fwm_s.sw_constraints = Param(fwm_s.farm_ids, initialize=sw_calib_constr_subset, mutable=True)
    fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset, mutable=True)
    x_start_values = dict(enumerate([0.0] * 3))
    fwm_s.xs_total = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
    fwm_s.xs_sw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
    fwm_s.xs_gw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
    fwm_s.nirs = Param(fwm_s.ids, initialize=nirs_subset, mutable=True)

    # define pyomo objective function
    def obj_fun(fwm_s):
        return 0.00001 * sum(sum((fwm_s.net_prices_total[h] * fwm_s.xs_total[h] - 0.5 * fwm_s.gammas_total[h] *
                                  fwm_s.xs_total[h] * fwm_s.xs_total[h]) for h in fwm_s.crop_ids_by_farm[f]) +
                             sum((fwm_s.net_prices_sw[i] * fwm_s.xs_sw[i]) for i in fwm_s.crop_ids_by_farm[f]) +
                             sum((fwm_s.net_prices_gw[g] * fwm_s.xs_gw[g]) for g in fwm_s.crop_ids_by_farm[f]) for f
                             in fwm_s.farm_ids)

    fwm_s.obj_f = Objective(rule=obj_fun, sense=maximize)

    # define constraints
    def land_constraint(fwm_s, ff):
        return sum(fwm_s.xs_total[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
               fwm_s.land_constraints[ff]

    fwm_s.c1 = Constraint(fwm_s.farm_ids, rule=land_constraint)

    def obs_lu_constraint_sum(fwm_s, i):
        return fwm_s.xs_sw[i] + fwm_s.xs_gw[i] == fwm_s.xs_total[i]

    fwm_s.c5 = Constraint(fwm_s.ids, rule=obs_lu_constraint_sum)

    def sw_constraint(fwm_s, ff):
        return sum(fwm_s.xs_sw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
               fwm_s.sw_constraints[ff]

    fwm_s.c2 = Constraint(fwm_s.farm_ids, rule=sw_constraint)

    def gw_constraint(fwm_s, ff):
        return sum(fwm_s.xs_gw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
               fwm_s.gw_constraints[ff]

    fwm_s.c6 = Constraint(fwm_s.farm_ids, rule=gw_constraint)

    fwm_s.scaling_factor = Suffix(direction=Suffix.EXPORT)
    fwm_s.scaling_factor[fwm_s.xs_total] = 0.0001
    fwm_s.scaling_factor[fwm_s.xs_sw] = 0.0001
    fwm_s.scaling_factor[fwm_s.land_constraints] = 0.000001
    fwm_s.scaling_factor[fwm_s.sw_constraints] = 0.01
    fwm_s.scaling_factor[fwm_s.gw_constraints] = 0.01
    fwm_s.scaling_factor[fwm_s.c1] = 0.000001
    fwm_s.scaling_factor[fwm_s.c5] = 0.000001
    fwm_s.scaling_factor[fwm_s.c2] = 0.01
    fwm_s.scaling_factor[fwm_s.c6] = 0.01

    # create and run the optimization solver
    # opt = SolverFactory("ipopt", solver_io='nl')
    opt = SolverFactory("gurobi", solver_io='python')
    results = opt.solve(fwm_s, keepfiles=False, tee=True)
    print(results.solver.termination_condition)

    # store main model outputs
    result_xs_sw = dict(fwm_s.xs_sw.get_values())
    result_xs_gw = dict(fwm_s.xs_gw.get_values())
    result_xs_total = dict(fwm_s.xs_total.get_values())

    farms_list = [farms_master.loc[[farm_id[0]]].nldas.values[0]]

    # process groundwater production results and
    results_pd = farms_master[farms_master['nldas'].isin(farms_list)]
    results_pd['xs_gw'] = 0
    results_pd['xs_sw'] = 0
    results_pd['xs_total'] = 0
    results_pd['id'] = results_pd['index']
    results_xs_sw_pd = pd.DataFrame.from_dict(result_xs_sw, orient='index')
    results_xs_sw_pd['id'] = results_xs_sw_pd.index + 1
    results_xs_sw_pd = results_xs_sw_pd.rename(columns={0: "xs_sw_temp"})
    results_pd = results_pd.merge(results_xs_sw_pd[['id', 'xs_sw_temp']], how='left', on=['id'])
    results_pd.loc[results_pd['xs_sw_temp'].notnull(), 'xs_sw'] = results_pd['xs_sw_temp']
    results_xs_gw_pd = pd.DataFrame.from_dict(result_xs_gw, orient='index')
    results_xs_gw_pd['id'] = results_xs_gw_pd.index + 1
    results_xs_gw_pd = results_xs_gw_pd.rename(columns={0: "xs_gw_temp"})
    results_pd = results_pd.merge(results_xs_gw_pd[['id', 'xs_gw_temp']], how='left', on=['id'])
    results_pd.loc[results_pd['xs_gw_temp'].notnull(), 'xs_gw'] = results_pd['xs_gw_temp']
    results_xs_total_pd = pd.DataFrame.from_dict(result_xs_total, orient='index')
    results_xs_total_pd['id'] = results_xs_total_pd.index + 1
    results_xs_total_pd = results_xs_total_pd.rename(columns={0: "xs_total_temp"})
    results_pd = results_pd.merge(results_xs_total_pd[['id', 'xs_total_temp']], how='left', on=['id'])
    results_pd.loc[results_pd['xs_total_temp'].notnull(), 'xs_total'] = results_pd['xs_total_temp']
    results_pd = results_pd.drop(['xs_gw_temp', 'xs_sw_temp', 'xs_total_temp'], axis=1)
    results_pd['year'] = t+1
    results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected']
    try:
        results_pd['gw_mult'] = gw_multiplier
    except NameError:
        results_pd['gw_mult'] = 1
    try:
        results_pd['gw_cumul_vol'] = cumul_gw_sum
    except NameError:
        results_pd['gw_cumul_vol'] = 0
    results_pd['gw_run'] = gw_cost_id
    results_pd['trans'] = gw_cost_lookup[(gw_cost_lookup.row_indx == gw_cost_id)]['Trans'].values[0]
    if first:
        results_combined = results_pd
    else:
        results_combined = pd.concat([results_combined, results_pd])
    aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
    gw_sum = results_pd[['nldas','gw_vol','xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
    gw_sum = gw_sum.reset_index()
    # gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
    # gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
    gw_sum['gw_vol_km3'] = gw_sum['gw_vol'] * 1.23348e-6
    if first:
        cumul_gw_sum = gw_sum['gw_vol_km3'].values[0]
    else:
        cumul_gw_sum += gw_sum['gw_vol_km3'].values[0]
    i = 0

    if gw_curve_option == 'external':
        for index, row in gw_cost_curves.iterrows():
            if cumul_gw_sum > row['volume']:
                i = index
            else:
                break
        if gw_cost_curves.loc[i][gw_cost_id] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
            gw_multiplier = gw_cost_curves.loc[i][gw_cost_id]
        else:
            gw_multiplier = 9999999999999  # Set groundwater cost to extremely high value to reflect groundwater exhaustion

    elif gw_curve_option == 'internal':
        for index in range(gw_cost_curve_internal[0][0].size):
            if cumul_gw_sum > gw_cost_curve_internal[1][0][index]:
                i = index
            else:
                break
        if gw_cost_curve_internal[0][0][i] != 0:
            gw_multiplier = gw_cost_curve_internal[0][0][i] / gw_cost_curve_internal[0][0][0]
        else:
            gw_multiplier = 9999999999999

    print(cumul_gw_sum)
    first = False

#  Calculate summary results (percent change in total irrigated area as an initial result)
aggregation_functions = {'xs_total': 'sum'}
summary_pd = results_combined[['nldas','year','xs_total']].groupby(['nldas','year']).aggregate(aggregation_functions)
summary_pd = summary_pd.reset_index()
year_start = summary_pd.year.min()
year_end = summary_pd.year.max()
end_div_start_area = summary_pd[(summary_pd.year==year_end)].xs_total.values[0] / summary_pd[(summary_pd.year == year_start)].xs_total.values[0]


###########################

import networkx as nx
import numpy as np
import itertools
import matplotlib.pyplot as plt


def drawgraphs(SAresults):
    # Get list of parameters
    parameters = list(SAresults['S1'].keys())
    # Set min index value, for the effects to be considered significant
    # index_significance_value = 0.01
    index_significance_value = 0.000001 # JY TEMP
    '''
    Define some general layout settings.
    '''
    node_size_min = 15  # Max and min node size
    node_size_max = 30
    border_size_min = 1  # Max and min node border thickness
    border_size_max = 8
    edge_width_min = 1  # Max and min edge thickness
    edge_width_max = 10
    edge_distance_min = 0.1  # Max and min distance of the edge from the center of the circle
    edge_distance_max = 0.6  # Only applicable to the curved edges

    '''
    Set up some variables and functions that will facilitate drawing circles and 
    moving items around.
    '''
    # Define circle center and radius
    center = [0.0, 0.0]
    radius = 1.0

    # Function to get distance between two points
    def distance(p1, p2):
        return np.sqrt(((p1 - p2) ** 2).sum())

    # Function to get middle point between two points
    def middle(p1, p2):
        return (p1 + p2) / 2

    # Function to get the vertex of a curve between two points
    def vertex(p1, p2, c):
        m = middle(p1, p2)
        curve_direction = c - m
        return m + curve_direction * (edge_distance_min + edge_distance_max * (1 - distance(m, c) / distance(c, p1)))

    # Function to get the angle of the node from the center of the circle
    def angle(p, c):
        # Get x and y distance of point from center
        [dx, dy] = p - c
        if dx == 0:  # If point on vertical axis (same x as center)
            if dy > 0:  # If point is on positive vertical axis
                return np.pi / 2.
            else:  # If point is on negative vertical axis
                return np.pi * 3. / 2.
        elif dx > 0:  # If point in the right quadrants
            if dy >= 0:  # If point in the top right quadrant
                return np.arctan(dy / dx)
            else:  # If point in the bottom right quadrant
                return 2 * np.pi + np.arctan(dy / dx)
        elif dx < 0:  # If point in the left quadrants
            return np.pi + np.arctan(dy / dx)

    '''
    First, set up graph with all parameters as nodes and draw all second order (S2)
    indices as edges in the network. For every S2 index, we need a Source parameter,
    a Target parameter, and the Weight of the line, given by the S2 index itself. 
    '''
    combs = [list(c) for c in list(itertools.combinations(parameters, 2))]

    Sources = list(list(zip(*combs))[0])
    Targets = list(list(zip(*combs))[1])
    # Sometimes computing errors produce negative Sobol indices. The following reads
    # in all the indices and also ensures they are between 0 and 1.
    Weights = [max(min(x, 1), 0) for x in [SAresults['S2'][Sources[i]][Targets[i]] for i in range(len(Sources))]]
    Weights = [0 if x < index_significance_value else x for x in Weights]

    # Set up graph
    G = nx.Graph()
    # Draw edges with appropriate weight
    for s, t, weight in zip(Sources, Targets, Weights):
        G.add_edges_from([(s, t)], w=weight)

    # Generate dictionary of node postions in a circular layout
    Pos = nx.circular_layout(G)

    '''
    Normalize node size according to first order (S1) index. First, read in S1 indices,
    ensure they're between 0 and 1 and normalize them within the max and min range
    of node sizes.
    Then, normalize edge thickness according to S2. 
    '''
    # Node size
    first_order = [max(min(x, 1), 0) for x in [SAresults['S1'][key] for key in SAresults['S1']]]
    first_order = [0 if x < index_significance_value else x for x in first_order]
    node_size = [node_size_min * (1 + (node_size_max - node_size_min) * k / max(first_order)) for k in first_order]
    # Node border thickness
    total_order = [max(min(x, 1), 0) for x in [SAresults['ST'][key] for key in SAresults['ST']]]
    total_order = [0 if x < index_significance_value else x for x in total_order]
    border_size = [border_size_min * (1 + (border_size_max - border_size_min) * k / max(total_order)) for k in
                   total_order]
    # Edge thickness
    edge_width = [edge_width_min * ((edge_width_max - edge_width_min) * k / max(Weights)) for k in Weights]

    '''
    Draw network. This will draw the graph with straight lines along the edges and 
    across the circle. 
    '''
    nx.draw_networkx_nodes(G, Pos, node_size=node_size, node_color='#98B5E2',
                           edgecolors='#1A3F7A', linewidths=border_size)
    nx.draw_networkx_edges(G, Pos, width=edge_width, edge_color='#2E5591', alpha=0.7)
    names = nx.draw_networkx_labels(G, Pos, font_size=12, font_color='#0B2D61', font_family='sans-serif')
    for node, text in names.items():
        position = (radius * 1.1 * np.cos(angle(Pos[node], center)), radius * 1.1 * np.sin(angle(Pos[node], center)))
        text.set_position(position)
        text.set_clip_on(False)
    plt.gcf().set_size_inches(9, 9)  # Make figure a square
    plt.axis('off')

    '''
     We can now draw the network with curved lines along the edges and across the circle.
     Calculate all distances between 1 node and all the others (all distances are 
     the same since they're in a circle). We'll need this to identify the curves 
     we'll be drawing along the perimeter (i.e. those that are next to each other).
     '''
    min_distance = round(min([distance(Pos[list(G.nodes())[0]], Pos[n]) for n in list(G.nodes())[1:]]), 1)

    # Figure to generate the curved edges between two points
    def xy_edge(p1, p2):  # Point 1, Point 2
        m = middle(p1, p2)  # Get middle point between the two nodes
        # If the middle of the two points falls very close to the center, then the
        # line between the two points is simply straight
        if distance(m, center) < 1e-6:
            xpr = np.linspace(p1[0], p2[0], 10)
            ypr = np.linspace(p1[1], p2[1], 10)
        # If the distance between the two points is the minimum (i.e. they are next
        # to each other), draw the edge along the perimeter
        elif distance(p1, p2) <= min_distance:
            # Get angles of two points
            p1_angle = angle(p1, center)
            p2_angle = angle(p2, center)
            # Check if the points are more than a hemisphere apart
            if max(p1_angle, p2_angle) - min(p1_angle, p2_angle) > np.pi:
                radi = np.linspace(max(p1_angle, p2_angle) - 2 * np.pi, min(p1_angle, p2_angle))
            else:
                radi = np.linspace(min(p1_angle, p2_angle), max(p1_angle, p2_angle))
            xpr = radius * np.cos(radi) + center[0]
            ypr = radius * np.sin(radi) + center[1]
            # Otherwise, draw curve (parabola)
        else:
            edge_vertex = vertex(p1, p2, center)
            a = distance(edge_vertex, m) / ((distance(p1, p2) / 2) ** 2)
            yp = np.linspace(-distance(p1, p2) / 2, distance(p1, p2) / 2, 100)
            xp = a * (yp ** 2)
            xp += distance(center, edge_vertex)
            theta_m = angle(middle(p1, p2), center)
            xpr = np.cos(theta_m) * xp - np.sin(theta_m) * yp
            ypr = np.sin(theta_m) * xp + np.cos(theta_m) * yp
            xpr += center[0]
            ypr += center[1]
        return xpr, ypr

    '''
    Draw network. This will draw the graph with curved lines along the edges and 
    across the circle. 
    '''
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(1, 1, 1)
    for i, e in enumerate(G.edges()):
        x, y = xy_edge(Pos[e[0]], Pos[e[1]])
        ax.plot(x, y, '-', c='#2E5591', lw=edge_width[i], alpha=0.7)
    for i, n in enumerate(G.nodes()):
        ax.plot(Pos[n][0], Pos[n][1], 'o', c='#98B5E2', markersize=node_size[i] / 5, markeredgecolor='#1A3F7A',
                markeredgewidth=border_size[i] * 1.15)

    for i, text in enumerate(G.nodes()):
        if node_size[i] < 100:
            position = (
            radius * 1.05 * np.cos(angle(Pos[text], center)), radius * 1.05 * np.sin(angle(Pos[text], center)))
        else:
            position = (
            radius * 1.01 * np.cos(angle(Pos[text], center)), radius * 1.01 * np.sin(angle(Pos[text], center)))
        plt.annotate(text, position, fontsize=12, color='#0B2D61', family='sans-serif')
    ax.axis('off')
    fig.tight_layout()
    plt.show()

def S2_to_dict(matrix, problem):
    result = {}
    names = list(problem["names"])

    for i in range(problem["num_vars"]):
        for j in range(i + 1, problem["num_vars"]):
            if names[i] not in result:
                result[names[i]] = {}
            if names[j] not in result:
                result[names[j]] = {}

            result[names[i]][names[j]] = result[names[j]][names[i]] = float(matrix[i][j])

    return result

result = {} #create dictionary to store new
result['S1']={k : float(v) for k, v in zip(problem["names"], Si["S1"])}
result['S1_conf']={k : float(v) for k, v in zip(problem["names"], Si["S1_conf"])}
result['S2'] = S2_to_dict(Si['S2'], problem)
result['S2_conf'] = S2_to_dict(Si['S2_conf'], problem)
result['ST']={k : float(v) for k, v in zip(problem["names"], Si["ST"])}
result['ST_conf']={k : float(v) for k, v in zip(problem["names"], Si["ST_conf"])}

drawgraphs(result)

#################################


sw_mult = 1
price_mult = 1
gammas_mult = 1
K = 10
alphas_mult = 1  # JY temp to disable alphas in sensitivity analysis
nir_mult = 1 # JY temp to disable alphas in sensitivity analysis

# generate gw cost curves if gw_curve_option is "internal"
if gw_curve_option == 'internal':
    gw_cost_curve_internal = Theis_pumping_with_deepening.Analytical(S, m, K, WL, R, IrrDepth, years)

# apply sensitivity multipliers
gammas_total_subset = gammas_total_subset_og.copy()
for key in gammas_total_subset:
    gammas_total_subset[key] *= gammas_mult

nirs_subset = nirs_subset_og.copy()
for key in nirs_subset:
    nirs_subset[key] *= nir_mult

sw_calib_constr_subset = sw_calib_constr_subset_og.copy()
for key in sw_calib_constr_subset:
    sw_calib_constr_subset[key] *= sw_mult

net_prices_land_subset = net_prices_land_subset_og.copy()
for key in net_prices_land_subset:
    net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - alphas_total_subset[key]) * price_mult

# for key in net_prices_land_subset:
#     net_prices_land_subset[key] = ((yields_subset[key] * prices_subset[key]) - land_costs_subset[key] - (alphas_total_subset[key] * alphas_mult))

# initialize counters/trackers
first = True
depletion_first = True
cumul_gw_sum = 0
gw_multiplier = 1

# start simulation time loop
for t in range(no_of_years):
    print('year' + str(t))
    print('!!JY!! ' + str(nir_mult))

    # initialize pyomo model
    fwm_s = ConcreteModel()
    fwm_s.ids = Set(initialize=ids_subset_sorted)
    if first:
        fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset, mutable=True)
    else:
        net_prices_gw_subset_update = net_prices_gw_subset.copy()
        net_prices_gw_subset_update.update((x, y * gw_multiplier) for x, y in net_prices_gw_subset_update.items())
        fwm_s.net_prices_gw = Param(fwm_s.ids, initialize=net_prices_gw_subset_update, mutable=True)
    fwm_s.farm_ids = Set(initialize=farm_id)
    fwm_s.crop_ids_by_farm = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
    fwm_s.crop_ids_by_farm_and_constraint = Set(fwm_s.farm_ids, initialize=crop_ids_by_farm_subset)
    fwm_s.net_prices_sw = Param(fwm_s.ids, initialize=net_prices_sw_subset, mutable=True)
    fwm_s.net_prices_total = Param(fwm_s.ids, initialize=net_prices_land_subset, mutable=True)
    fwm_s.gammas_total = Param(fwm_s.ids, initialize=gammas_total_subset, mutable=True)
    fwm_s.land_constraints = Param(fwm_s.farm_ids, initialize=max_land_constr_subset, mutable=True)
    fwm_s.sw_constraints = Param(fwm_s.farm_ids, initialize=sw_calib_constr_subset, mutable=True)
    fwm_s.gw_constraints = Param(fwm_s.farm_ids, initialize=gw_calib_constr_subset, mutable=True)
    x_start_values = dict(enumerate([0.0] * 3))
    fwm_s.xs_total = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
    fwm_s.xs_sw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
    fwm_s.xs_gw = Var(fwm_s.ids, domain=NonNegativeReals, initialize=x_start_values)
    fwm_s.nirs = Param(fwm_s.ids, initialize=nirs_subset, mutable=True)

    # define pyomo objective function
    def obj_fun(fwm_s):
        return 0.00001 * sum(sum((fwm_s.net_prices_total[h] * fwm_s.xs_total[h] - 0.5 * fwm_s.gammas_total[h] *
                                  fwm_s.xs_total[h] * fwm_s.xs_total[h]) for h in fwm_s.crop_ids_by_farm[f]) +
                             sum((fwm_s.net_prices_sw[i] * fwm_s.xs_sw[i]) for i in fwm_s.crop_ids_by_farm[f]) +
                             sum((fwm_s.net_prices_gw[g] * fwm_s.xs_gw[g]) for g in fwm_s.crop_ids_by_farm[f]) for f
                             in fwm_s.farm_ids)

    fwm_s.obj_f = Objective(rule=obj_fun, sense=maximize)

    # define constraints
    def land_constraint(fwm_s, ff):
        return sum(fwm_s.xs_total[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
               fwm_s.land_constraints[ff]

    fwm_s.c1 = Constraint(fwm_s.farm_ids, rule=land_constraint)

    def obs_lu_constraint_sum(fwm_s, i):
        return fwm_s.xs_sw[i] + fwm_s.xs_gw[i] == fwm_s.xs_total[i]

    fwm_s.c5 = Constraint(fwm_s.ids, rule=obs_lu_constraint_sum)

    def sw_constraint(fwm_s, ff):
        return sum(fwm_s.xs_sw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
               fwm_s.sw_constraints[ff]

    fwm_s.c2 = Constraint(fwm_s.farm_ids, rule=sw_constraint)

    def gw_constraint(fwm_s, ff):
        return sum(fwm_s.xs_gw[i] * fwm_s.nirs[i] for i in fwm_s.crop_ids_by_farm_and_constraint[ff]) <= \
               fwm_s.gw_constraints[ff]

    fwm_s.c6 = Constraint(fwm_s.farm_ids, rule=gw_constraint)

    fwm_s.scaling_factor = Suffix(direction=Suffix.EXPORT)
    fwm_s.scaling_factor[fwm_s.xs_total] = 0.0001
    fwm_s.scaling_factor[fwm_s.xs_sw] = 0.0001
    fwm_s.scaling_factor[fwm_s.land_constraints] = 0.000001
    fwm_s.scaling_factor[fwm_s.sw_constraints] = 0.01
    fwm_s.scaling_factor[fwm_s.gw_constraints] = 0.01
    fwm_s.scaling_factor[fwm_s.c1] = 0.000001
    fwm_s.scaling_factor[fwm_s.c5] = 0.000001
    fwm_s.scaling_factor[fwm_s.c2] = 0.01
    fwm_s.scaling_factor[fwm_s.c6] = 0.01

    # create and run the optimization solver
    # opt = SolverFactory("ipopt", solver_io='nl')
    opt = SolverFactory("gurobi", solver_io='python')
    results = opt.solve(fwm_s, keepfiles=False, tee=True)
    print(results.solver.termination_condition)

    # store main model outputs
    result_xs_sw = dict(fwm_s.xs_sw.get_values())
    result_xs_gw = dict(fwm_s.xs_gw.get_values())
    result_xs_total = dict(fwm_s.xs_total.get_values())

    farms_list = [farms_master.loc[[farm_id[0]]].nldas.values[0]]

    # process groundwater production results and
    results_pd = farms_master[farms_master['nldas'].isin(farms_list)]
    results_pd['xs_gw'] = 0
    results_pd['xs_sw'] = 0
    results_pd['xs_total'] = 0
    results_pd['id'] = results_pd['index']
    results_xs_sw_pd = pd.DataFrame.from_dict(result_xs_sw, orient='index')
    results_xs_sw_pd['id'] = results_xs_sw_pd.index + 1
    results_xs_sw_pd = results_xs_sw_pd.rename(columns={0: "xs_sw_temp"})
    results_pd = results_pd.merge(results_xs_sw_pd[['id', 'xs_sw_temp']], how='left', on=['id'])
    results_pd.loc[results_pd['xs_sw_temp'].notnull(), 'xs_sw'] = results_pd['xs_sw_temp']
    results_xs_gw_pd = pd.DataFrame.from_dict(result_xs_gw, orient='index')
    results_xs_gw_pd['id'] = results_xs_gw_pd.index + 1
    results_xs_gw_pd = results_xs_gw_pd.rename(columns={0: "xs_gw_temp"})
    results_pd = results_pd.merge(results_xs_gw_pd[['id', 'xs_gw_temp']], how='left', on=['id'])
    results_pd.loc[results_pd['xs_gw_temp'].notnull(), 'xs_gw'] = results_pd['xs_gw_temp']
    results_xs_total_pd = pd.DataFrame.from_dict(result_xs_total, orient='index')
    results_xs_total_pd['id'] = results_xs_total_pd.index + 1
    results_xs_total_pd = results_xs_total_pd.rename(columns={0: "xs_total_temp"})
    results_pd = results_pd.merge(results_xs_total_pd[['id', 'xs_total_temp']], how='left', on=['id'])
    results_pd.loc[results_pd['xs_total_temp'].notnull(), 'xs_total'] = results_pd['xs_total_temp']
    results_pd = results_pd.drop(['xs_gw_temp', 'xs_sw_temp', 'xs_total_temp'], axis=1)
    results_pd['year'] = t+1
    results_pd['gw_vol'] = results_pd['xs_gw'] * results_pd['nir_corrected']
    try:
        results_pd['gw_mult'] = gw_multiplier
    except NameError:
        results_pd['gw_mult'] = 1
    try:
        results_pd['gw_cumul_vol'] = cumul_gw_sum
    except NameError:
        results_pd['gw_cumul_vol'] = 0
    results_pd['gw_run'] = gw_cost_id
    results_pd['trans'] = gw_cost_lookup[(gw_cost_lookup.row_indx == gw_cost_id)]['Trans'].values[0]
    if first:
        results_combined = results_pd
    else:
        results_combined = pd.concat([results_combined, results_pd])
    aggregation_functions = {'gw_vol': 'sum', 'xs_gw': 'sum'}
    gw_sum = results_pd[['nldas','gw_vol','xs_gw']].groupby(['nldas']).aggregate(aggregation_functions)
    gw_sum = gw_sum.reset_index()
    # gw_sum = gw_sum.merge(farm_gw_sum[['nldas', 'no_of_wells']], how='left', on=['nldas'])
    # gw_sum['gw_vol_well'] = gw_sum['gw_vol'] * farms_per_grid / gw_sum['no_of_wells']
    gw_sum['gw_vol_km3'] = gw_sum['gw_vol'] * 1.23348e-6
    if first:
        cumul_gw_sum = gw_sum['gw_vol_km3'].values[0]
    else:
        cumul_gw_sum += gw_sum['gw_vol_km3'].values[0]
    i = 0

    if gw_curve_option == 'external':
        for index, row in gw_cost_curves.iterrows():
            if cumul_gw_sum > row['volume']:
                i = index
            else:
                break
        if gw_cost_curves.loc[i][gw_cost_id] != 0:  # JY temp to deal with zeros in groundwater cost curves (check in with Stephen)
            gw_multiplier = gw_cost_curves.loc[i][gw_cost_id]
        else:
            gw_multiplier = 9999999999999  # Set groundwater cost to extremely high value to reflect groundwater exhaustion

    elif gw_curve_option == 'internal':
        for index in range(gw_cost_curve_internal[0][0].size):
            if cumul_gw_sum > gw_cost_curve_internal[1][0][index]:
                i = index
            else:
                break
        if gw_cost_curve_internal[0][0][i] != 0:
            gw_multiplier = gw_cost_curve_internal[0][0][i] / gw_cost_curve_internal[0][0][0]
        else:
            gw_multiplier = 9999999999999
            if depletion_first:
                time_depletion = t
                depletion_first = False

    print(cumul_gw_sum)
    first = False

# Calculate
results_combined['profit_calc'] = ((price_mult * (results_combined['price'] * results_combined['yield']) -
                                         results_combined['land_cost'] - results_combined['alphas_land']) * results_combined['xs_total']) \
                                       - (0.5 * results_combined['gammas_total'] * results_combined['xs_total'] * results_combined['xs_total']) \
                                       - (results_combined['net_prices_gw'] * results_combined['xs_gw']) \
                                       - (results_combined['net_prices_sw'] * results_combined['xs_sw'])


#  Calculate summary results (percent change in total irrigated area as an initial result)
aggregation_functions = {'xs_total': 'sum', 'profit_calc': 'sum'}
summary_pd = results_combined[['nldas','year','xs_total','profit_calc']].groupby(['nldas','year']).aggregate(aggregation_functions)
summary_pd = summary_pd.reset_index()
year_start = summary_pd.year.min()
year_end = summary_pd.year.max()
end_div_start_area = summary_pd[(summary_pd.year==year_end)].xs_total.values[0] / summary_pd[(summary_pd.year == year_start)].xs_total.values[0]
total_profit = summary_pd.profit_calc.sum()

# Volume depleted (calculated from cost curves)
perc_vol_depleted = cumul_gw_sum / max(gw_cost_curve_internal[1][0].tolist())
if perc_vol_depleted > 1:
    perc_vol_depleted = 1

