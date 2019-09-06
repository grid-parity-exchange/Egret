import os
from collections import namedtuple, defaultdict
import datetime
import textwrap

import matplotlib as mpl
## catch when we're running linux without X
if os.name == 'posix' and 'DISPLAY' not in os.environ:
        mpl.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Seaborn/matplotlib plot settings
sns.set()
sns.set_context('paper', font_scale=2.00)


font = {'family' : 'sans-serif',
        'weight' : 'regular',
        'size'   : 14
        }
mpl.rc('font', **font)


GenerationType = namedtuple('GenerationType',
                            [
                           'label',
                           'color',
                            ]
                           )


GENERATION_TYPES = {
    'Z': GenerationType('Battery', '#42F1F4'),
    'N': GenerationType('Nuclear', '#b22222'),
    'E': GenerationType('Geothermal', '#CDE7B0'),
    'B': GenerationType('Biomass', '#A3BFA8'),
    'C': GenerationType('Coal', '#333333'),
    'G': GenerationType('Gas', '#6e8b3d'),
    'O': GenerationType('Oil', '#eea2ad'),
    'H': GenerationType('Hydro', '#add8e6'),
    'W': GenerationType('Wind', '#4f94cd'),
    'S': GenerationType('Solar', '#ffb90f'),
    'SC': GenerationType('SynchCond', '#fff68f'),
    'Other': GenerationType('Other', '#886688'),
}


FUEL_TO_CODE = defaultdict(lambda: 'Other')


BUILT_IN_FUEL_CODES = [
    ('Oil', 'O'),
    ('Coal', 'C'),
    ('Gas', 'G'),
    ('NG', 'G'),
    ('Solar', 'S'),
    ('PV', 'S'),
    ('Wind', 'W'),
    ('Nuclear', 'N'),
    ('Hydro', 'H'),
    ('Sync_Cond', 'SC'),
    ('Biomass', 'B'),
    # To accomodate data that uses the codes directly
    ('O', 'O'),
    ('C', 'C'),
    ('G', 'G'),
    ('S', 'S'),
    ('W', 'W'),
    ('N', 'N'),
    ('SC', 'SC'),
    ('B', 'B'),
]


for fuel_name, fuel_code in BUILT_IN_FUEL_CODES:
    FUEL_TO_CODE[fuel_name] = fuel_code


GENERATION_TYPE_SORT_KEY = {
    'Nuclear': 0,
    'Coal': 1,
    'Hydro': 2,
    'Gas': 3,
    'Oil': 4,
    'Sync_Cond': 5,
    'Wind': 6,
    'Solar': 7,
    'Biomass': 8,
    'Geothermal': 9,
    'Battery': 10,
    'Other': 999,
}


def _fuel_type_to_code(x):
    """Converts the fuel type string to its generation type code equivalent."""
    code = FUEL_TO_CODE.get(x, '')
    
    return code

def _build_attribute_to_array_func(time_indices):
    '''returns a function for converting EGRET time-valued objects to np arrays'''
    def attribute_to_array(attr):
        '''returns a numpy array for the time-valued attr'''
        if isinstance(attr, dict):
            return np.array([attr['values'][t] for t in time_indices])
        else:
            return np.array([attr for t in time_indices])
    return attribute_to_array

def generate_stack_graph(egret_model_data, bar_width=0.9, 
                            x_tick_frequency=1,
                            title='', 
                            plot_individual_generators=False,
                            show_individual_components=False):
    '''
    Creates a stack graph using an egret ModelData object solved using the egret.models.unit_commitment.solve_unit_commitment() function.

    Parameters
    ----------
    egret_model_data : egret.data.ModelData
        An egret ModelData object with the appropriate data loaded and the model solved.
    bar_width : float
        The width of each bar stack in the time series in (0, 1]; default is 0.9
    x_tick_frequency : int
        Indicates the frequency of labeling the time axis; default is 1
    title : str (optional)
        Title to put on the resulting graph; default is ''
    plot_individual_generators : bool (optional)
        If True, individual generator output will be plotted and labeled. Raises a ValueError if there are more than 5 generators in the model. Raises a ValueError if show_individual_components is simultaneously True; default is False
    show_individual_components : bool (optional)
        If True, individual generator output within a generation type will be discretely indicated. Raises a ValueError if plot_individual_generators is simultaneously True; default is False 

    '''

    ## functions for interfacing with EGRET time structure
    time_periods = egret_model_data.data['system']['time_indices']
    attribute_to_array = _build_attribute_to_array_func(time_periods)

    def _plot_generation_stack_components():
        if plot_individual_generators and show_individual_components:
            raise ValueError('plot_individual_generators and show_individual_components cannot be simultaneously True.')
            return
        
        bottom = np.zeros(len(indices))
        
        if plot_individual_generators:      
            INDIVIDUAL_GEN_PLOT_UPPER_LIMIT = 5

            if len(egret_model_data.data['elements']['generator']) > INDIVIDUAL_GEN_PLOT_UPPER_LIMIT:
                raise ValueError('There are too many generators in the system to support plotting output individually. (maximum: {0})'.format(INDIVIDUAL_GEN_PLOT_UPPER_LIMIT))

            for generator, generator_data in egret_model_data.data['elements']['generator'].items():
                pg_array = attribute_to_array(generator_data['pg'])

                if not sum(pg_array) > 0.0:
                    continue

                is_quickstart = generator_data.get('quickstart_capable', False)

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
                pg_array = attribute_to_array(generator_data['pg'])

                if not sum(pg_array) > 0.0:
                    continue

                reported_fuel_type = generator_data['fuel']

                # Match fuel type to one of pre-determined types or 'other'
                fuel_type = GENERATION_TYPES[FUEL_TO_CODE[reported_fuel_type]].label

                is_quickstart = generator_data.get('quickstart_capable', False)

                if is_quickstart:
                    quickstart_label = 'quickstart'
                else:
                    quickstart_label = 'not quickstart'
        
                try:
                    generator_generation_by_fuel_type[fuel_type].append((pg_array, quickstart_label))
                except KeyError:
                    generator_generation_by_fuel_type[fuel_type] = []
                    generator_generation_by_fuel_type[fuel_type].append((pg_array, quickstart_label))
            
            sorted_total_generation_by_fuel_type = sorted(generator_generation_by_fuel_type.items(), key=lambda x: GENERATION_TYPE_SORT_KEY.get(x[0], 1e3))

            for generation_type, generator_output_levels in sorted_total_generation_by_fuel_type:
                try:
                    component_label = GENERATION_TYPES[generation_type].label
                except KeyError:
                    component_label = GENERATION_TYPES[FUEL_TO_CODE[generation_type]].label
                
                try:
                    component_color = GENERATION_TYPES[generation_type].color
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
                pg_array = attribute_to_array(generator_data['pg'])

                if not sum(pg_array) > 0.0:
                    continue
    
                reported_fuel_type = generator_data['fuel']

                # Match fuel type to one of pre-determined types or 'other'
                fuel_type = GENERATION_TYPES[FUEL_TO_CODE[reported_fuel_type]].label

                is_quickstart = generator_data.get('quickstart_capable', False)

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
            sorted_total_generation_by_fuel_type = sorted(total_generation_by_fuel_type.items(), key=lambda x: GENERATION_TYPE_SORT_KEY.get(x[0], 1e3))

            for generation_type, total_generation_array in sorted_total_generation_by_fuel_type:
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

    fig, ax = plt.subplots(figsize=(16, 8))

    time_labels = [textwrap.fill(time_index, 10) for time_index in egret_model_data.data['system']['time_indices']]
    indices = np.arange(len(time_labels))
                
    # Plot generation dispatch/output.
    total_dispatch_levels = _plot_generation_stack_components()
    bottom = total_dispatch_levels
    
    # Plot load shedding, if applicable.    
    bus_dict = egret_model_data.data['elements']['bus']
    total_load_shed_by_hour = np.zeros(len(indices))

    for bus, bus_data in bus_dict.items():
        p_balance_violation = attribute_to_array(bus_data['p_balance_violation'])
        load_shed = np.maximum(0, p_balance_violation)

        total_load_shed_by_hour += load_shed
    
    if sum(total_load_shed_by_hour) > 0.0:
        component_color = '#ffff00'
        ax.bar(indices, total_load_shed_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0,                
               label='Load Shed')
        bottom += total_load_shed_by_hour   

    # Plot demand.
    demand_by_hour = np.zeros(len(indices))
    load_structure = egret_model_data.data['elements']['load']

    for load_bus in load_structure:
        bus_load = load_structure[load_bus]['p_load']
        load_array = attribute_to_array(bus_load)
        demand_by_hour += load_array
    
    ## This is to make it so the step graph covers the total dispatch levels as expected.
    demand_indices = np.arange(len(time_labels)+1) - 1/2
    demand_by_hour = np.append(demand_by_hour, demand_by_hour[-1])
    
    ax.step(demand_indices, demand_by_hour, linewidth=3, color='#000000', where='post')
    
    # Add reserve requirements, if applicable.
    reserve_requirements_by_hour = egret_model_data.data['system'].get('reserve_requirement')
    if reserve_requirements_by_hour is not None:
        reserve_requirements_array = attribute_to_array(reserve_requirements_by_hour)

        if sum(reserve_requirements_array) > 0.0:
            component_color = '#00c2ff'
            ax.bar(indices, reserve_requirements_array, bar_width, bottom=bottom, color=component_color,
                   edgecolor=None, linewidth=0,
                   label='Required Reserve')
            bottom += reserve_requirements_array
    
        # Add reserve shortfalls, if applicable.
        reserve_shortfall_by_hour = egret_model_data.data['system']['reserve_shortfall']
        reserve_shortfall_array = attribute_to_array(reserve_shortfall_by_hour)
    
        if sum(reserve_shortfall_array) > 0.0:
            component_color = '#ff00ff'
            ax.bar(indices, reserve_shortfall_array, bar_width, bottom=bottom, color=component_color,
                   edgecolor=None, linewidth=0,
                   label='Reserve Shortfall')
            bottom += reserve_shortfall_array
    
    # Add implicit reserves, if applicable.
    generators_dict = egret_model_data.data['elements']['generator']
    reserves_by_hour = np.zeros(len(indices))

    for _, gen_data in generators_dict.items():
        is_quickstart = gen_data.get('quickstart_capable', False)

        if gen_data['generator_type'] == 'thermal':
            headroom = attribute_to_array(gen_data['headroom'])

            if is_quickstart:
                commitment = attribute_to_array(gen_data['commitment'])
                reserves_available = [headroom[ix] if commit_status > 0.0 else 0 for ix, commit_status in enumerate(commitment)]
            else:
                reserves_available = headroom

            reserves_by_hour += reserves_available
    
    if reserve_requirements_by_hour is not None:
        implicit_reserves_by_hour = [max(0.0, reserves_by_hour[ix] - reserves_by_hour[ix]) for ix in range(len(reserve_requirements_by_hour))]
    else:
        implicit_reserves_by_hour = [max(0.0, reserves_by_hour[ix]) for ix in range(len(reserves_by_hour))]
    
    if sum(implicit_reserves_by_hour) > 0.0:
        component_color = '#00ffc7'
        ax.bar(indices, implicit_reserves_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Implicit Reserve')
        bottom += implicit_reserves_by_hour
    
    # Add quickstart capacity, if applicable.
    total_quickstart_capacity_by_hour = np.zeros(len(indices))

    for _, gen_data in generators_dict.items():
        is_quickstart = gen_data.get('quickstart_capable', False)

        if is_quickstart:
            commitment = attribute_to_array(gen_data['commitment'])
            p_max = gen_data['p_max']
            startup_capacity = gen_data['startup_capacity']
            quickstart_capacity = min(p_max, startup_capacity)

            quickstart_capacity_available = [quickstart_capacity if commit_status > 0.0 else 0 for ix, commit_status in enumerate(commitment)]

            total_quickstart_capacity_by_hour += quickstart_capacity_available
    
    if sum(total_quickstart_capacity_by_hour) > 0.0:
        component_color = '#494949'
        ax.bar(indices, total_quickstart_capacity_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               hatch='//',
               label='Available Quickstart')
        bottom += total_quickstart_capacity_by_hour
    
    # Add renewable curtailment.
    total_renewable_curtailment_by_hour = np.zeros(len(indices))

    for _, gen_data in generators_dict.items():
        generator_type = gen_data['generator_type']

        if generator_type == 'renewable':
            p_max = attribute_to_array(gen_data['p_max'])
            pg = attribute_to_array(gen_data['pg'])

            p_curtailed = np.maximum(p_max - pg, 0)

            total_renewable_curtailment_by_hour += p_curtailed
    
    if sum(total_renewable_curtailment_by_hour) > 0.0:
        component_color = '#ff0000'
        ax.bar(indices, total_renewable_curtailment_by_hour, bar_width, bottom=bottom, color=component_color, 
               edgecolor=None, linewidth=0, 
               label='Renewables Curtailed')
        bottom += total_renewable_curtailment_by_hour
    
    # Labels and such.
    plt.xticks(indices[::x_tick_frequency], time_labels[::x_tick_frequency], rotation=0)
    ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    
    # sns.despine(offset=10, trim=True)

    # Put legend outside on the right.
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_title(title)
    ax.set_ylabel('Power [MW]')
    ax.set_xlabel('Time')
    ax.yaxis.grid(True)

    return fig, ax


def main():
    import json
    from egret.models.unit_commitment import solve_unit_commitment, create_tight_unit_commitment_model

    current_dir = os.path.dirname(os.path.abspath(__file__))

    TEST_WITH_RTS_GMLC = True

    if TEST_WITH_RTS_GMLC:
        # Test using RTS-GMLC case, if available
        from egret.parsers.rts_gmlc_parser import create_ModelData

        rts_gmlc_dir = os.path.join(current_dir, '..', '..', '..', 'RTS-GMLC')
        begin_time = "2020-07-01"
        end_time = "2020-07-02"

        md = create_ModelData(
            rts_gmlc_dir, begin_time, end_time, 
            # simulation="DAY_AHEAD", t0_state = None,
            )
        
        solved_md = solve_unit_commitment(md,
                        'gurobi',
                        mipgap = 0.001,
                        timelimit = None,
                        solver_tee = True,
                        symbolic_solver_labels = False,
                        options = None,
                        uc_model_generator=create_tight_unit_commitment_model,
                        relaxed=False,
                        return_model=False)

        fig, ax = generate_stack_graph(
            solved_md, 
            title=begin_time,
            show_individual_components=False,
            plot_individual_generators=False,
            x_tick_frequency=4,
        )
    else:
        ## Test using built-in unit commitment unit test case(s)
        from egret.data.model_data import ModelData

        test_cases = [os.path.join(current_dir, '..', 'models', 'tests', 'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1, 6)]

        for test_case in test_cases:   
            with open(test_case, 'r') as f:
                md_dict = json.load(f)
            md = ModelData(md_dict)

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

            fig, ax = generate_stack_graph(
                solved_md, 
                title=repr(test_case),
                show_individual_components=False,
                plot_individual_generators=False,
            )
    
    plt.show()


if __name__ == '__main__':
    main()
