import os
from collections import namedtuple
import datetime
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Seaborn/matplotlib plot settings
sns.set()
sns.set_context('paper', font_scale=2.15)

plt.rcParams["figure.figsize"] = (12,6)
font = {'family' : 'sans-serif',
        'weight' : 'regular',
        'size'   : 18
        }
mpl.rc('font', **font)

GenerationType = namedtuple('GenerationType',
                            [
                           'label',
                           'color',
                            ]
                           )
# Colors from MAGMA (R color names)
GENERATION_TYPES = {
    'Z': GenerationType('Battery', '#42F1F4'),  # no color assigned
    'N': GenerationType('Nuclear', '#b22222'),
    'E': GenerationType('Geothermal', '#CDE7B0'),  # no color assigned
    'B': GenerationType('Biomass', '#A3BFA8'),  # no color assigned
    'C': GenerationType('Coal', '#333333'),
    'G': GenerationType('Gas', '#6e8b3d'),
    'O': GenerationType('Oil', '#eea2ad'),
    'H': GenerationType('Hydro', '#add8e6'),
    'W': GenerationType('Wind', '#4f94cd'),
    'S': GenerationType('Solar', '#ffb90f'),
    'SC': GenerationType('SynchCond', '#fff68f')
}

FUEL_TO_CODE = {
    'Oil': 'O',
    'Coal': 'C',
    'NG': 'G',
    'Nuclear': 'N',
    'Hydro': 'H',
    'Solar': 'S',
    'Wind': 'W',
    'Sync_Cond': 'SC',
               }

def _fuel_type_to_code(x):
    code = FUEL_TO_CODE.get(x, '')
    
    return code


def generate_stack_graph(thermal_data, renewables_data, hourly_summary, gen_summary, bar_width=0.9, title='', 
                         plot_individual_generators=False,
                         show_individual_components=False):
    """Creates a stack graph."""    
    def _compute_total_dispatch_per_hour(gen_type):
        """Computes the total dispatch for gen_type generation for each hour in the data set."""
        total_output_by_hour = []
        
        if gen_type in {'S', 'W', 'H'}:
            gen_details = renewables_data
            
            for h in gen_details['Hour'].unique():
                total_hour_output = gen_details.loc[(gen_details['gen type'] == gen_type) & (gen_details['Hour'] == h)]['Output'].sum()
                total_output_by_hour.append(total_hour_output)
        else:
            gen_details = thermal_data
            
            for h in gen_details['Hour'].unique():
                total_hour_output = gen_details.loc[(gen_details['gen type'] == gen_type) & (gen_details['Hour'] == h)]['Dispatch'].sum()
                total_output_by_hour.append(total_hour_output)
                
        return total_output_by_hour

    def _plot_generation_stack_components():
        """
        Plots the stack graph components for generation. 
        If plot_individual_generators is True, each individual generator's output will be keyed.
        If show_individual_components is True, components of each generation type will be broken out. This option is only relevant when plot_individual_generators is False.
        """
        if plot_individual_generators and show_individual_components:
            raise Exception('plot_individual_generators and show_individual_components cannot be simultaneously True.')
            return
        
        bottom = np.zeros(len(indices))
        
        if plot_individual_generators:      
            for gen_type, gen_type_tuple in GENERATION_TYPES.items():
                if gen_type in {'S', 'W', 'H'}:
                    gen_details = renewables_data
                    output_field_name = 'Output'
                else:
                    gen_details = thermal_data
                    output_field_name = 'Dispatch'

                generator_list = gen_details.loc[gen_details['gen type'] == gen_type]['Generator'].unique()

                for ix, gen in enumerate(generator_list, start=0):
                    gen_hourly_output = gen_details.loc[(gen_details['Generator'] == gen)][output_field_name].values

                    ax.bar(indices, gen_hourly_output, bar_width, bottom=bottom, 
                           edgecolor=None, linewidth=0, 
                           label=gen)

                    # Update "bottom" position for each stack. 
                    bottom += gen_hourly_output        
        elif show_individual_components:
            # Plot by individual generator in groups of generation type.
            for gen_type, gen_type_tuple in GENERATION_TYPES.items():
                component_color = gen_type_tuple.color
                vals = _compute_total_dispatch_per_hour(gen_type)

                if sum(vals) <= 0.0:
                    continue

                if gen_type in {'S', 'W', 'H'}:
                    gen_details = renewables_data
                    output_field_name = 'Output'
                else:
                    gen_details = thermal_data
                    output_field_name = 'Dispatch'

                generator_list = gen_details.loc[gen_details['gen type'] == gen_type]['Generator'].unique()

                for ix, gen in enumerate(generator_list, start=0):
                    gen_hourly_output = gen_details.loc[(gen_details['Generator'] == gen)][output_field_name].values

                    # Plot set of bar stack component.
                    if ix == 0:
                        ax.bar(indices, gen_hourly_output, bar_width, bottom=bottom, color=component_color, label=gen_type_tuple.label,
                              edgecolor='#FFFFFF', linewidth=0.5)
                    else:
                        ax.bar(indices, gen_hourly_output, bar_width, bottom=bottom, color=component_color,
                              edgecolor='#FFFFFF', linewidth=0.5)

                    # Update "bottom" position for each stack. 
                    bottom += gen_hourly_output
        else:    
            for gen_type, gen_type_tuple in GENERATION_TYPES.items():
                component_color = gen_type_tuple.color
                vals = _compute_total_dispatch_per_hour(gen_type)

                # Skip if none of the generation type is present.
                if sum(vals) <= 0.0:
                    continue

                # Plot set of bar stack component.
                ax.bar(indices, vals, bar_width, bottom=bottom, color=component_color, label=gen_type_tuple.label,
                      linewidth=0)

                # Update "bottom" position for each stack. 
                bottom += vals
        
        return bottom

    fig, ax = plt.subplots(figsize=(16, 10))

    indices = np.arange(24)

    def _generate_time_labels():    
        thermal_data['delta'] = thermal_data['Hour'].apply(lambda x: datetime.timedelta(hours=x))
        thermal_data['dt'] = pd.to_datetime(thermal_data['Date']) + thermal_data['delta']
        labels = thermal_data['dt'].unique()
        time_labels = [textwrap.fill(pd.to_datetime(x).strftime('%b %d %H:%M'), 6) for x in labels]

        # labels = [str(ix+1) for ix in indices]

        return time_labels
    
    time_labels = _generate_time_labels()
                
    # Plot generation dispatch/output.
    total_dispatch_levels = _plot_generation_stack_components()
    bottom = total_dispatch_levels
    
    # Plot load shedding, if applicable.    
    load_shed_by_hour = hourly_summary[' LoadShedding '].values
    
    if sum(load_shed_by_hour) > 0.0:
        component_color = '#ffff00'
        ax.bar(indices, load_shed_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0,                
               label='Load Shedding')
        bottom += load_shed_by_hour   

    # Plot demand.
    demand_by_hour = hourly_summary[' Demand '].values
    
    ## This is to make it so the step graph covers the total dispatch levels as expected.
    demand_indices = np.arange(25) - 1/2
    demand_by_hour = np.append(demand_by_hour, demand_by_hour[-1])
    
    ax.step(demand_indices, demand_by_hour, linewidth=3, color='#000000', where='post')
    
    # Add reserve requirements, if applicable.
    # TODO: Don't have this explicitly in the output files
    reserve_requirements_by_hour = np.zeros(len(indices))
    #reserve_requirements_by_hour = [0.15*hour_demand for hour_demand in demand_by_hour]
    
    if sum(reserve_requirements_by_hour) > 0.0:
        component_color = '#00c2ff'
        ax.bar(indices, reserve_requirements_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Required Reserve')
        bottom += reserve_requirements_by_hour
    
    # Add reserve shortfalls, if applicable.
    reserve_shortfall_by_hour = hourly_summary[' ReserveShortfall '].values
    #reserve_shortfall_by_hour = np.random.randint(0, 500, size=len(indices))
    
    if sum(reserve_shortfall_by_hour) > 0.0:
        component_color = '#ff00ff'
        ax.bar(indices, reserve_shortfall_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Reserve Shortfall')
        bottom += reserve_shortfall_by_hour
    
    # Add implicit reserves, if applicable.
    # TODO:
    reserve_available_by_hour = gen_summary['Available reserves'].values
    implicit_reserves_by_hour = [max(0.0, reserve_available_by_hour[ix] - reserve_requirements_by_hour[ix]) for ix in range(len(reserve_available_by_hour))]

    if sum(reserve_requirements_by_hour) > 0.0:
        component_color = '#00ffc7'
        ax.bar(indices, implicit_reserves_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Implicit Reserve')
        bottom += implicit_reserves_by_hour
    
    # Add quick-start capacity, if applicable.
    # TODO:
    quickstart_capacity_by_hour = np.zeros(len(indices))
    #quickstart_capacity_by_hour = np.random.randint(0, 500, size=len(indices))
    
    if sum(quickstart_capacity_by_hour) > 0.0:
        component_color = '#494949'
        ax.bar(indices, quickstart_capacity_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Available Quick Start')
        bottom += quickstart_capacity_by_hour
    
    # TODO:
    quickstart_additional_power_by_hour = np.zeros(len(indices))
    #quickstart_additional_power_by_hour = np.random.randint(0, 400, size=len(indices))
    
    if sum(quickstart_additional_power_by_hour) > 0.0:
        component_color = '#0061ff'
        ax.bar(indices, quickstart_additional_power_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor='#000000', linewidth=0, hatch='xxx', 
               label='Quick-Start Generator Output')
        bottom += quickstart_additional_power_by_hour
        
        # mpl.rcParams['hatch.color'] = 'red'
        # mpl.rcParams['hatch.linewidth'] = 1.0
    
    # Add renewable curtailment.
    renewable_curtailment_by_hour = np.zeros(len(indices))
    renewable_curtailment_by_hour = hourly_summary[' RenewablesCurtailed '].values
    
    if sum(renewable_curtailment_by_hour) > 0.0:
        component_color = '#ff0000'
        ax.bar(indices, renewable_curtailment_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Renewables Curtailed')
        bottom += renewable_curtailment_by_hour
    
    # Labels and such.
    tick_frequency = 4
    plt.xticks(indices[::tick_frequency], time_labels[::tick_frequency], rotation=0)
    ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    
    sns.despine(offset=10, trim=True)

    # Put legend outside on the right.
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_title(title)
    ax.set_ylabel('Power [MW]')
    ax.set_xlabel('Hour')
    ax.yaxis.grid(True)

    plt.savefig(title+'.png', format='png')

    return fig, ax


def generate_stack_graph_egret(egret_model_data, bar_width=0.9, 
                                x_tick_frequency=1,
                                title='', 
                                plot_individual_generators=False,
                                show_individual_components=False):
    """Creates a stack graph."""
    def _indexed_dict_to_array(indexed_dict):
        """Converts dictionary with keys of time indices and values of time series values to a NumPy array."""
        return np.array([p[-1] for p in sorted(indexed_dict.items(), key=lambda x: int(x[-1]))])

    def _compute_total_dispatch_per_hour():
        """Computes the total dispatch for gen_type generation for each hour in the data set."""
        total_output_by_hour = []


        
        if gen_type in {'S', 'W', 'H'}:
            gen_details = renewables_data
            
            for h in gen_details['Hour'].unique():
                total_hour_output = gen_details.loc[(gen_details['gen type'] == gen_type) & (gen_details['Hour'] == h)]['Output'].sum()
                total_output_by_hour.append(total_hour_output)
        else:
            gen_details = thermal_data
            
            for h in gen_details['Hour'].unique():
                total_hour_output = gen_details.loc[(gen_details['gen type'] == gen_type) & (gen_details['Hour'] == h)]['Dispatch'].sum()
                total_output_by_hour.append(total_hour_output)
                
        return total_output_by_hour

    def _plot_generation_stack_components():
        """
        Plots the stack graph components for generation. 
        If plot_individual_generators is True, each individual generator's output will be keyed.
        If show_individual_components is True, components of each generation type will be broken out. This option is only relevant when plot_individual_generators is False.
        """
        if plot_individual_generators and show_individual_components:
            raise Exception('plot_individual_generators and show_individual_components cannot be simultaneously True.')
            return
        
        bottom = np.zeros(len(indices))
        
        if plot_individual_generators:      
            INDIVIDUAL_GEN_PLOT_UPPER_LIMIT = 5

            if len(egret_model_data.data['elements']['generator']) > INDIVIDUAL_GEN_PLOT_UPPER_LIMIT:
                raise Exception('There are too many generators in the system to support plotting output individually. (maximum: {0})'.format(INDIVIDUAL_GEN_PLOT_UPPER_LIMIT))

            for generator, generator_data in egret_model_data.data['elements']['generator'].items():
                is_quickstart = generator_data.get('quickstart_capable', False)

                pg = generator_data['pg']['values']
                pg_array = _indexed_dict_to_array(pg)

                if is_quickstart:
                    quickstart_label = 'quickstart'
                else:
                    quickstart_label = 'not quickstart'
                
                if is_quickstart:
                    ax.bar(indices, pg_array, bar_width, bottom=bottom, label=generator,
                        hatch='//',
                        edgecolor='#FFFFFF', linewidth=0.5)
                else:
                    ax.bar(indices, pg_array, bar_width, bottom=bottom, label=generator,
                        edgecolor='#FFFFFF', linewidth=0.5)

                bottom += pg_array        
        elif show_individual_components:
            # Plot by individual generator in groups of generation type.
            generator_generation_by_fuel_type = {}

            for generator, generator_data in egret_model_data.data['elements']['generator'].items():
                fuel_type = generator_data['fuel']
                is_quickstart = generator_data.get('quickstart_capable', False)

                pg = generator_data['pg']['values']
                pg_array = _indexed_dict_to_array(pg)

                if is_quickstart:
                    quickstart_label = 'quickstart'
                else:
                    quickstart_label = 'not quickstart'
        
                try:
                    generator_generation_by_fuel_type[fuel_type].append((pg_array, quickstart_label))
                except KeyError:
                    generator_generation_by_fuel_type[fuel_type] = []
                    generator_generation_by_fuel_type[fuel_type].append((pg_array, quickstart_label))
            
            for fuel_type, generator_output_levels in generator_generation_by_fuel_type.items():
                try:
                    component_label = GENERATION_TYPES[fuel_type].label
                except KeyError:
                    component_label = GENERATION_TYPES[FUEL_TO_CODE[generation_type]].label
                
                try:
                    component_color = GENERATION_TYPES[fuel_type].color
                except KeyError:
                    component_color = GENERATION_TYPES[FUEL_TO_CODE[generation_type]].color

                if len(generator_output_levels) < 1:
                    continue
                else:
                    sorted_generator_output_levels = sorted(generator_output_levels, key=lambda x: x[-1] == 'quickstart')

                    for ix, output_level_tuple in enumerate(sorted_generator_output_levels, start=0):
                        output_level_array = output_level_tuple[0]
                        is_quickstart = True if output_level_tuple[-1] == 'quickstart' else False

                        if is_quickstart:
                            if ix == 0:
                                ax.bar(indices, output_level_array, bar_width, bottom=bottom, color=component_color, label=component_label,
                                    hatch='//',
                                    edgecolor='#FFFFFF', linewidth=0.5)
                            else:
                                ax.bar(indices, output_level_array, bar_width, bottom=bottom, color=component_color,
                                    hatch='//',
                                    edgecolor='#FFFFFF', linewidth=0.5)
                        else:
                            if ix == 0:
                                ax.bar(indices, output_level_array, bar_width, bottom=bottom, color=component_color, label=component_label,
                                    edgecolor='#FFFFFF', linewidth=0.5)
                            else:
                                ax.bar(indices, output_level_array, bar_width, bottom=bottom, color=component_color,
                                    edgecolor='#FFFFFF', linewidth=0.5)

                        bottom += output_level_array
        else:    
            total_generation_by_fuel_type = {}

            for generator, generator_data in egret_model_data.data['elements']['generator'].items():
                fuel_type = generator_data['fuel']
                is_quickstart = generator_data.get('quickstart_capable', False)

                pg = generator_data['pg']['values']
                pg_array = _indexed_dict_to_array(pg)

                if is_quickstart:
                    quickstart_label = 'quickstart'
                else:
                    quickstart_label = 'not quickstart'
        
                try:
                    total_generation_by_fuel_type[fuel_type][quickstart_label] += pg_array
                except KeyError:
                    total_generation_by_fuel_type[fuel_type] = {}
                    total_generation_by_fuel_type[fuel_type][quickstart_label] = pg_array

            # Plot each bar stack component.
            for generation_type, total_generation_array in total_generation_by_fuel_type.items():
                try:
                    component_label = GENERATION_TYPES[generation_type].label
                except KeyError:
                    component_label = GENERATION_TYPES[FUEL_TO_CODE[generation_type]].label
                
                try:
                    component_color = GENERATION_TYPES[generation_type].color
                except KeyError:
                    component_color = GENERATION_TYPES[FUEL_TO_CODE[generation_type]].color

                # Non-quickstart.
                try:
                    component_values = total_generation_array['not quickstart']
                except KeyError:
                    pass
                else:
                    ax.bar(indices, component_values, bar_width, bottom=bottom, color=component_color, label=component_label,
                        linewidth=0)
                    bottom += component_values

                # Quickstart.
                try:
                    component_values = total_generation_array['quickstart']
                except KeyError:
                    pass
                else:
                    ax.bar(indices, component_values, bar_width, bottom=bottom, color=component_color, label=component_label,
                        hatch='//',
                        linewidth=0)
                    bottom += component_values
        
        return bottom

    fig, ax = plt.subplots(figsize=(16, 10))

    indices = np.arange(24)

    def _generate_time_labels():    
        # thermal_data['delta'] = thermal_data['Hour'].apply(lambda x: datetime.timedelta(hours=x))
        # thermal_data['dt'] = pd.to_datetime(thermal_data['Date']) + thermal_data['delta']
        # labels = thermal_data['dt'].unique()
        # time_labels = [textwrap.fill(pd.to_datetime(x).strftime('%b %d %H:%M'), 6) for x in labels]

        time_labels = [str(ix+1) for ix in indices]

        return time_labels
    
    time_labels = _generate_time_labels()
                
    # Plot generation dispatch/output.
    total_dispatch_levels = _plot_generation_stack_components()
    bottom = total_dispatch_levels
    
    # Plot load shedding, if applicable.    
    # TODO:
    load_shed_by_hour = np.zeros(len(indices))
    
    if sum(load_shed_by_hour) > 0.0:
        component_color = '#ffff00'
        ax.bar(indices, load_shed_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0,                
               label='Load Shedding')
        bottom += load_shed_by_hour   

    # Plot demand.
    demand_by_hour = np.zeros(len(indices))
    load_structure = egret_model_data.data['elements']['load']

    for load_bus in load_structure:
        bus_load = load_structure[load_bus]['p_load']['values']
        load_array = _indexed_dict_to_array(bus_load)
        demand_by_hour += load_array
    
    ## This is to make it so the step graph covers the total dispatch levels as expected.
    demand_indices = np.arange(25) - 1/2
    demand_by_hour = np.append(demand_by_hour, demand_by_hour[-1])
    
    ax.step(demand_indices, demand_by_hour, linewidth=3, color='#000000', where='post')
    
    # Add reserve requirements, if applicable.
    reserve_requirements_by_hour = egret_model_data.data['system']['reserve_requirement']['values']
    reserve_requirements_array = _indexed_dict_to_array(reserve_requirements_by_hour)
    
    if sum(reserve_requirements_array) > 0.0:
        component_color = '#00c2ff'
        ax.bar(indices, reserve_requirements_array, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Required Reserve')
        bottom += reserve_requirements_array
    
    # Add reserve shortfalls, if applicable.
    reserve_shortfall_by_hour = egret_model_data.data['system']['reserve_shortfall']['values']
    reserve_shortfall_array = _indexed_dict_to_array(reserve_shortfall_by_hour)
    
    if sum(reserve_shortfall_array) > 0.0:
        component_color = '#ff00ff'
        ax.bar(indices, reserve_shortfall_array, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Reserve Shortfall')
        bottom += reserve_shortfall_array
    
#     # Add implicit reserves, if applicable.
#     # TODO:
#     reserve_available_by_hour = gen_summary['Available reserves'].values
#     implicit_reserves_by_hour = [max(0.0, reserve_available_by_hour[ix] - reserve_requirements_by_hour[ix]) for ix in range(len(reserve_available_by_hour))]

#     if sum(reserve_requirements_by_hour) > 0.0:
#         component_color = '#00ffc7'
#         ax.bar(indices, implicit_reserves_by_hour, bar_width, bottom=bottom, color=component_color, 
#                edgecolor=None, linewidth=0, 
#                label='Implicit Reserve')
#         bottom += implicit_reserves_by_hour
    
#     # Add quick-start capacity, if applicable.
#     # TODO:
#     quickstart_capacity_by_hour = np.zeros(len(indices))
# #     quickstart_capacity_by_hour = np.random.randint(0, 500, size=len(indices))
    
#     if sum(quickstart_capacity_by_hour) > 0.0:
#         component_color = '#494949'
#         ax.bar(indices, quickstart_capacity_by_hour, bar_width, bottom=bottom, color=component_color, 
#                edgecolor=None, linewidth=0, 
#                label='Available Quick Start')
#         bottom += quickstart_capacity_by_hour
    
#     # TODO:
#     quickstart_additional_power_by_hour = np.zeros(len(indices))
# #     quickstart_additional_power_by_hour = np.random.randint(0, 400, size=len(indices))
    
#     if sum(quickstart_additional_power_by_hour) > 0.0:
#         component_color = '#0061ff'
#         ax.bar(indices, quickstart_additional_power_by_hour, bar_width, bottom=bottom, color=component_color, 
#                edgecolor='#000000', linewidth=0, hatch='xxx', 
#                label='Quick-Start Generator Output')
#         bottom += quickstart_additional_power_by_hour
        
# #         mpl.rcParams['hatch.color'] = 'red'
# #         mpl.rcParams['hatch.linewidth'] = 1.0
    
#     # Add renewable curtailment.
#     renewable_curtailment_by_hour = np.zeros(len(indices))
#     renewable_curtailment_by_hour = hourly_summary[' RenewablesCurtailed '].values
    
#     if sum(renewable_curtailment_by_hour) > 0.0:
#         component_color = '#ff0000'
#         ax.bar(indices, renewable_curtailment_by_hour, bar_width, bottom=bottom, color=component_color, 
#                edgecolor=None, linewidth=0, 
#                label='Renewables Curtailed')
#         bottom += renewable_curtailment_by_hour
    
    # Labels and such.
    plt.xticks(indices[::x_tick_frequency], time_labels[::x_tick_frequency], rotation=0)
    ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    
    sns.despine(offset=10, trim=True)

    # Put legend outside on the right.
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_title(title)
    ax.set_ylabel('Power [MW]')
    ax.set_xlabel('Hour')
    ax.yaxis.grid(True)

    # plt.savefig(title+'.png', format='png')

    return fig, ax


def main():
    import json
    import math
    from egret.models.unit_commitment import solve_unit_commitment, create_tight_unit_commitment_model
    from egret.data.model_data import ModelData

    case_prefix = 'rts-gmlc_deterministic'

    RESULTS_DIR = os.path.join('egret', 'viz', 'data', case_prefix)

    test_cases = [os.path.join('egret', 'models', 'tests', 'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1,2)]
    test_int_objvals = [4201915.017320504, 5454367.7670904165, 5999272.361123627, 5461120.3231092375, 6062406.32677043]

    from egret.parsers.rts_gmlc_parser import create_ModelData

    rts_gmlc_dir = os.path.join('..', 'RTS-GMLC')
    begin_time = "2020-07-01"
    end_time = "2020-07-02"

    md = create_ModelData(
        rts_gmlc_dir, begin_time, end_time, 
        # simulation="DAY_AHEAD", t0_state = None,
        )
    
    solved_md = solve_unit_commitment(md,
                    'cbc',
                    mipgap = 0.001,
                    timelimit = None,
                    solver_tee = True,
                    symbolic_solver_labels = False,
                    options = None,
                    uc_model_generator=create_tight_unit_commitment_model,
                    relaxed=False,
                    return_model=False)

    fig, ax = generate_stack_graph_egret(
        solved_md, 
        title=begin_time,
        show_individual_components=False,
        plot_individual_generators=False,
    )

    print(solved_md)

    # for test_case, ref_objval in zip(test_cases, test_int_objvals):
        
    #     md_dict = json.load(open(test_case,'r'))
    #     md = ModelData(md_dict)

    #     solved_md = solve_unit_commitment(md,
    #                       'cbc',
    #                       mipgap = 0.001,
    #                       timelimit = None,
    #                       solver_tee = True,
    #                       symbolic_solver_labels = False,
    #                       options = None,
    #                       uc_model_generator=create_tight_unit_commitment_model,
    #                       relaxed=False,
    #                       return_model=False)

    #     fig, ax = generate_stack_graph_egret(
    #         solved_md, 
    #         title=repr(test_case),
    #         show_individual_components=False,
    #         plot_individual_generators=False,
    #     )


    # # Set the directory to save figures
    # FIGS_DIR = os.path.join('report', 'figs')

    # RTS_SOURCE_DIR = os.path.join('..', 'RTS-GMLC', 'RTS_DATA', 'SourceData')
    # gens_table = pd.read_csv(os.path.join(RTS_SOURCE_DIR, 'gen.csv'))

    # gens_table['gen type'] = gens_table['Fuel'].apply(_fuel_type_to_code)
    # gens_table.set_index('GEN UID', inplace=True)
    # gens_dict = gens_table.to_dict(orient='index')

    # hourly_summary = pd.read_csv(os.path.join(RESULTS_DIR, 'hourly_summary.csv'))
    # hourly_gen_summary = pd.read_csv(os.path.join(RESULTS_DIR, 'Hourly_gen_summary.csv'))
    # thermal_detail = pd.read_csv(os.path.join(RESULTS_DIR, 'thermal_detail.csv'))
    # renewables_detail = pd.read_csv(os.path.join(RESULTS_DIR, 'renewables_detail.csv'))

    # thermal_detail['gen type'] = thermal_detail['Generator'].apply(lambda x: gens_dict[x]['gen type'])
    # renewables_detail['gen type'] = renewables_detail['Generator'].apply(lambda x: gens_dict[x]['gen type'])

    # for date in thermal_detail['Date'].unique()[-1:]:
    #     # Generate sample data for the day
    #     SAMPLE_THERMAL = thermal_detail.loc[thermal_detail['Date'] == date]
    #     SAMPLE_RENEWABLES = renewables_detail.loc[renewables_detail['Date'] == date]
    #     SAMPLE_HOURLY = hourly_summary.loc[hourly_summary['Date '] == date]
    #     SAMPLE_GENS = hourly_gen_summary.loc[hourly_gen_summary['Date'] == date]
        
    #     generate_stack_graph(SAMPLE_THERMAL, SAMPLE_RENEWABLES, SAMPLE_HOURLY, SAMPLE_GENS, 
    #                         plot_individual_generators=True, show_individual_components=False, 
    #                         title='Power Generation for {0}'.format(date))
    
    plt.show()


if __name__ == '__main__':
    main()