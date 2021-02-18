#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
  from typing import Dict, Union

import os.path
import pandas as pd
from datetime import datetime, timedelta
from collections import namedtuple

import egret.data.model_data as md

from .parsed_cache import ParsedCache
from ._reserves import is_valid_reserve_name

def create_ModelData(rts_gmlc_dir:str, 
                     begin_time:Union[datetime,str], end_time:Union[datetime,str], 
                     simulation:str="DAY_AHEAD", t0_state:dict = None):

    """
    Create a ModelData object from the RTS-GMLC data.

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to RTS-GMLC directory
    begin_time : datetime.datetime or str
        Beginning of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    end_time : datetime.datetime or str
        End of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    simulation : str
        Either "DAY_AHEAD" or "REAL_TIME", which specifies which time series the data is taken from, 
        default is "DAY_AHEAD".
    t0_state : dict or Nonetype
        Keys of this dict are thermal generator names, each element of which is another dictionary with
        keys "initial_status", "initial_p_output", and "initial_q_output", which specify whether the
        generator is on at t0, the real power output at t0, and the reactive power output at t0. 
        If this is None, default values are loaded.
    
    Returns
    -------
    egret.model_data.ModelData
        Returns a ModelData object with the timeseries data specified
    """
    return md.ModelData(create_model_data_dict(rts_gmlc_dir, begin_time, end_time, simulation, t0_state))

def create_model_data_dict(rts_gmlc_dir:str, 
                           begin_time:Union[datetime,str], end_time:Union[datetime,str],
                           simulation:str="DAY_AHEAD", t0_state:dict = None):

    """
    Create a model_data dictionary from the RTS-GMLC data.

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to RTS-GMLC directory
    begin_time : datetime.datetime or str
        Beginning of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    end_time : datetime.datetime or str
        End of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    simulation : str
        Either "DAY_AHEAD" or "REAL_TIME", which specifies which time series the data is taken from, 
        default is "DAY_AHEAD".
    t0_state : dict or Nonetype
        Keys of this dict are thermal generator names, each element of which is another dictionary with
        keys "initial_status", "initial_p_output", and "initial_q_output", which specify whether the
        generator is on at t0, the real power output at t0, and the reactive power output at t0. 
        If this is None, default values are loaded.
    
    Returns
    -------
        dict : A dictionary in the format required for the ModelData object.
    """
    cache = parse_to_cache(rts_gmlc_dir, begin_time, end_time)
    model = cache.generate_model(simulation, begin_time, end_time)
    if t0_state is not None:
        for name, gen in model['elements']['generator']:
            if gen['generator_type']=='thermal':
                gen['initial_status'] = t0_state[name]['initial_status']
                gen['initial_p_output'] = t0_state[name]['initial_p_output']
                gen['initial_q_output'] = t0_state[name]['initial_q_output']
    return model


def parse_to_cache(rts_gmlc_dir:str, 
                   begin_time:datetime,
                   end_time:datetime) -> ParsedCache:
    ''' Parse data in RTS-GMLC format, keeping the portions between a start and end time

    rts_gmlc_dir : str
        Path to RTS-GMLC directory
    begin_time : datetime.datetime or str
        Beginning of time horizon. 
    end_time : datetime.datetime or str
        End of time horizon.
    simulation : str
        Either "DAY_AHEAD" or "REAL_TIME", which specifies which time series the data is taken from, 
        default is "DAY_AHEAD".
    '''
    if not os.path.exists(rts_gmlc_dir):
        raise ValueError(f'RTS-GMLC directory "{rts_gmlc_dir}" does not exist')

    # Create the skeleton
    model_data = _create_rtsgmlc_skeleton(rts_gmlc_dir)

    # Save the data frequencies
    metadata_df = _read_metadata(rts_gmlc_dir)
    minutes_per_period = {'DAY_AHEAD':int(metadata_df.loc['Period_Resolution', 'DAY_AHEAD'])//60,
                          'REAL_TIME':int(metadata_df.loc['Period_Resolution', 'REAL_TIME'])//60}

    data_start, data_end = _get_data_date_range(metadata_df)

    # TODO: Validate begin_time and end_time

    timeseries_df = _read_timeseries_data(model_data, rts_gmlc_dir,
                           begin_time, end_time, minutes_per_period)

    load_participation_factors = _compute_bus_load_participation_factors(model_data)

    return ParsedCache(model_data, begin_time, end_time,
                       minutes_per_period['DAY_AHEAD'], minutes_per_period['REAL_TIME'], 
                       timeseries_df, load_participation_factors)
    
    
def _read_metadata(base_dir:str) -> pd.DataFrame:
    metadata_path = os.path.join(base_dir, "simulation_objects.csv")
    if not os.path.exists(metadata_path):
        raise ValueError(f'RTS-GMLC directory "{rts_gmlc_dir}" does not contain expected CSV files.')

    # Read metadata about the data
    metadata_df = pd.read_csv(metadata_path, index_col=0)

    return metadata_df

def _get_data_date_range(metadata_df):
    ''' Get the range of dates for which there is data available
    '''
    import dateutil.parser

    # Data start time
    row = metadata_df.loc['Date_From']
    data_start = max(dateutil.parser.parse(row['DAY_AHEAD']),
                     dateutil.parser.parse(row['REAL_TIME']))

    # Data end time is a little more tricky
    row = metadata_df.loc['Date_To']
    def _extract_end_date(which:str):
        # The actual end date is the metadata's Date_To plus a number of look ahead periods.
        # Each look ahead period is a specified number of seconds
        extra_seconds = int(metadata_df.loc['Look_Ahead_Periods_per_Step'][which]) * \
                        int(metadata_df.loc['Look_Ahead_Resolution'][which])
        end_date = dateutil.parser.parse(row[which])
        return end_date + timedelta(seconds=extra_seconds)
    # Get the end date for each kind of data.  Both kinds of data
    # are available up through the later of the two simulation categories.
    data_end = max(_extract_end_date('DAY_AHEAD'),
                   _extract_end_date('REAL_TIME'))

    return (data_start, data_end)

def _read_timeseries_file(file_name:str, minutes_per_period:int,
                          start_time:datetime, end_time:datetime,
                          series_name:str) -> pd.DataFrame:
    """ 
    Read data from a timeseries file, returning only the data that falls within the requested time range

    Parameters
    ----------
    file_name:str
        Path to a CVS file with timeseries data
    minutes_per_period:int
        The number of minutes between time periods in the data
    start_time:datetime
        The earliest time to include in the returned data
    end_time:datetime
        The first time to NOT include in the returned data
    series_name:str
        The name the column holding the resulting data

    Returns
    -------
        df : A DataFrame with data between the specified dates

    Timeseries files can be in one of two formats, columnar or 2D.  A columnar file can hold more 
    than one timeseries in the same file, with one row per time period and one column per series.
    A 2D timeseries file holds only one timeseries, with one row per day and one column per time
    period within the day.  This function reads from either format.
    
    If the indicated file is columnar, the returned DataFrame will include all columns in the file, 
    one column per time series, with columns named as they appear in the file.

    If the indicated file is 2D, the returned DataFrame will have a single column whose name
    is assigned to be the `series_name` passed into the function.
    """
    # Determine which layout we're working with by checking column headers
    headers = pd.read_csv(file_name, nrows=0).columns.tolist()
    if headers[3] == 'Period':
        return _read_columnar_timeseries_file(file_name, minutes_per_period, start_time, end_time)
    else:
        df = _read_2D_timeseries_file(file_name, minutes_per_period, start_time, end_time, series_name)
        return df

def _read_columnar_timeseries_file(file_name:str, minutes_per_period:int,
                                   start_time:datetime, end_time:datetime) -> pd.DataFrame:
    """ 
    Read data from a timeseries file, returning only the data that falls within the requested time range

    Parameters
    ----------
    file_name:str
        Path to a CVS file with timeseries data
    minutes_per_period:int
        The number of minutes between time periods in the data
    start_time:datetime
        The earliest time to include in the returned data
    end_time:datetime
        The first time to NOT include in the returned data

    Returns
    -------
        df : A DataFrame with data between the specified dates

    The returned DataFrame converts the first 4 columns into a datetime which is used as
    the DataFrame's index.  All other CSV columns are included as columns in the DataFrame.
    """
    _date_parser = lambda *columns: datetime(*map(int,columns[0:3])) + \
                                    timedelta(minutes = minutes_per_period*(int(columns[3])-1))
    df = pd.read_csv(file_name, 
                     header=0, 
                     parse_dates=[[0, 1, 2, 3]],
                     date_parser=_date_parser,
                     index_col=0)
    df.index.names = ['DateTime']

    df.sort_index(inplace=True)

    # Remove data outside requested time period.
    # DataFrame slices include the end of the slice, so we need to reduce it slightly to
    # avoid including an extra value at the end.
    end_time = end_time - timedelta(seconds=1)
    df = df[start_time:end_time]

    # Be sure to return a copy instead of a view into the larger data set
    return df.copy()

def _read_2D_timeseries_file(file_name:str, minutes_per_period:int, 
                             start_time:datetime, end_time:datetime, 
                             column_name:str) -> DataFrame:
    """ 
    Read data from a timeseries file with a 2D layout, returning only the data that falls within the requested time range

    Parameters
    ----------
    file_name:str
        Path to a CVS file with reserve-formatted timeseries data
    start_time:datetime
        The earliest time to include in the returned data
    end_time:datetime
        The first time to NOT include in the returned data
    minutes_per_period:int
        The number of minutes between time periods in the data
    column_name:str
        The name that should be given to the DataFrame's column

    Returns
    -------
        df : A single-column pandas DataFrame with data between the specified dates

    The returned DataFrame has one row per cell in the original CSV, indexed by the cell's 
    datetime.  The first 3 columns indicate the date for all cells in the row, and other 
    column headers indicate the time within the day.  Only data within the requested time 
    period is included in the results.  Like a typical python range, the returned data includes 
    the start_time but does not include the end_time.
    """
    _date_parser = lambda *columns: datetime(*map(int,columns[0:3]))
    df = pd.read_csv(file_name, 
                     header=0, 
                     parse_dates=[[0, 1, 2]],
                     date_parser=_date_parser,
                     index_col=0)
    df.sort_index(inplace=True)

    # Remove data outside requested time period.
    # DataFrame slices include the end of the slice, so we need to reduce it slightly to
    # avoid including an extra value at the end.
    end_time = end_time - timedelta(seconds=1)
    df = df[start_time:end_time]

    # Now divide rows into one row per column, resetting index to appropriate datetime
    s = df.stack()
    s.index = map(lambda i: i[0] + timedelta(minutes=minutes_per_period*(int(i[1])-1)), s.index)
    s.index.names = ['DateTime']

    # Filter one more time, trimming out times of day that fall before or after the requested range
    if s.index[0] < start_time or s.index[-1] >= end_time:
        s = s[start_time:end_time]

    # Create and return a new 1-column DataFrame
    return pd.DataFrame({column_name: s})

def _create_rtsgmlc_skeleton(rts_gmlc_dir):
    """
    Creates a data dictionary from the RTS-GMLC data files, without loading hourly data

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to RTS-GMLC directory
    
    Returns
    -------
    data : dict
        Returns a dict loaded from the RTS-GMLC data
    """

    base_dir = rts_gmlc_dir

    model_data = md.ModelData.empty_model_data_dict()

    elements = model_data["elements"]
    system = model_data["system"]

    system["name"] = "RTS-GMLC"

    # this is the default used in the MATPOWER writer for RTS-GMLC
    system["baseMVA"] = 100.

    elements["bus"] = {}
    elements["load"] = {}
    elements["shunt"] = {}

    # add the buses
    bus_types = {'PQ': 'PQ', 
                 'PV': 'PV', 
                 'Ref': 'ref'}
    bus_id_to_name = {}
    bus_areas = set()
    bus_df = pd.read_csv(os.path.join(base_dir,'bus.csv'))
    for idx,row in bus_df.iterrows():
        BUS_TYPE = row['Bus Type']
        if not BUS_TYPE in bus_types:
            raise ValueError(f'Encountered an unsupported bus type: "{BUS_TYPE}" when parsing RTS-GMLC input file')

        bus_name = str(row['Bus Name'])
        bus_dict = {
            "id": str(row['Bus ID']),
            "base_kv": float(row['BaseKV']),
            "matpower_bustype": bus_types[BUS_TYPE],
            "vm": float(row['V Mag']),
            "va": float(row['V Angle']),
            "v_min": 0.95,
            "v_max": 1.05,
            "area": str(row['Area']),
            "zone": str(int(row['Zone'])),
        }

        if bus_dict["base_kv"] <= 0:
            raise ValueError(f'BaseKV value for bus "{bus_name}" is <= 0. Not supported.')

        PD = float(row['MW Load'])
        QD = float(row['MVAR Load'])
        if PD != 0 or QD != 0:
            load_dict = {
                "bus": bus_name, 
                "in_service":True,
                "p_load": PD,
                "q_load": QD,
                "area": bus_dict['area'],
                "zone": bus_dict['zone']
            }
            elements["load"][bus_name] = load_dict

        GS = float(row['MW Shunt G'])
        BS = float(row['MVAR Shunt B'])        
        if GS != 0 or BS != 0:
            shunt_dict = {
                "shunt_type":"fixed", 
                "bus": bus_name,
                "gs": GS,
                "bs": BS
            }
            elements["shunt"][bus_name] = shunt_dict

        if BUS_TYPE == 'Ref':
            va = bus_dict['va']
            if va != 0:
                if abs(va) >= 1e-16:
                    raise ValueError('EGRET only supports reference buses with an angle of 0 degrees.')
                msg = "\nEgret only supports reference buses with an angle of 0 degrees. \nFound a " \
                              "reference bus with an angle close to 0. \n" \
                              "Value: {va}\nSetting reference bus angle to 0."
                warnings.warn(msg)
                bus_dict['va'] = 0.0
            system["reference_bus"] = bus_name
            system["reference_bus_angle"] = 0

        bus_id_to_name[bus_dict['id']] = bus_name
        bus_areas.add(bus_dict['area'])
        elements['bus'][bus_name] = bus_dict

    # add the areas
    elements['area'] = {name:dict() for name in bus_areas}

    # add the branches
    elements["branch"] = {}
    branch_df = pd.read_csv(os.path.join(base_dir,'branch.csv'))
    for idx,row in branch_df.iterrows():

        branch_dict = {
            "from_bus": bus_id_to_name[str(row['From Bus'])], 
            "to_bus": bus_id_to_name[str(row['To Bus'])],
            "in_service": True,
            "resistance": float(row['R']),
            "reactance": float(row['X']),
            "charging_susceptance": float(row['B']),
            "rating_long_term": float(row['Cont Rating']) or None,
            "rating_short_term": float(row['LTE Rating']) or None,
            "rating_emergency": float(row['STE Rating']) or None,
            "angle_diff_min": -90,
            "angle_diff_max": 90,
            "pf": None,
            "qf": None,
            "pt": None,
            "qt": None
        }

        TAP = float(row['Tr Ratio'])
        if TAP != 0.0:
            branch_dict["branch_type"] = "transformer"
            branch_dict["transformer_tap_ratio"] = TAP
            branch_dict["transformer_phase_shift"] = 0.0
        else:
            branch_dict["branch_type"] = "line"

        name = str(row['UID'])
        elements["branch"][name] = branch_dict

    # add the DC branches
    if os.path.exists(os.path.join(base_dir,'dc_branch.csv')):
        branch_df = pd.read_csv(os.path.join(base_dir,'dc_branch.csv'))
        for idx,row in branch_df.iterrows():

            # TODO: I have no idea what field names Egrets expects or supports for DC branches.
            #       The code below is just a placeholder.
            branch_dict = {
                "from_bus": bus_id_to_name[str(row['From Bus'])], 
                "to_bus": bus_id_to_name[str(row['To Bus'])],
                "in_service": True,
                "branch_type": "dc",
                "resistance": float(row['R Line'])
            }

            name = str(row['UID'])
            elements["branch"][name] = branch_dict

    # add the generators
    elements["generator"] = {}
    RENEWABLE_TYPES = {'WIND', 'HYDRO', 'RTPV', 'PV'}
    gen_df = pd.read_csv(os.path.join(base_dir,'gen.csv'))
    for idx,row in gen_df.iterrows():
        # if this is storage we need to handle it differently
        if row['Fuel'] == 'Storage':
            continue

        # NOTE: for now, Egret doesn't handle CSP -- not clear how to model
        if row['Unit Type'] == 'CSP':
            continue

        name = str(row['GEN UID'])
        bus_name = bus_id_to_name[str(row['Bus ID'])]
        gen_dict = {
            "bus": bus_name,
            "in_service": True,
            "mbase": 100.0,
            "pg": float(row['MW Inj']),
            "qg": float(row['MVAR Inj']),
            "vg": float(row['V Setpoint p.u.']),
            "p_min": float(row['PMin MW']),
            "p_max": float(row['PMax MW']),
            "q_min": float(row['QMin MVAR']),
            "q_max": float(row['QMax MVAR']),
            "ramp_q": float(row['Ramp Rate MW/Min']),
            "fuel": str(row['Fuel']),
            "unit_type": str(row['Unit Type']),
            "area": elements['bus'][bus_name]['area'],
            "zone": elements['bus'][bus_name]['zone']
        }

        UNIT_TYPE = str(row['Unit Type'])
        if UNIT_TYPE in RENEWABLE_TYPES:
            gen_dict["generator_type"] = "renewable"
        elif UNIT_TYPE == 'SYNC_COND':
            ## TODO: should we have a flag for these?
            gen_dict["generator_type"] = "thermal" 
        else:
            gen_dict["generator_type"] = "thermal"

        elements["generator"][name] = gen_dict

        # after this is only really needed for thermal units
        if UNIT_TYPE in RENEWABLE_TYPES:
            continue

        # Gen cost
        ## round as in RTS-GMLC Prescient/topysp.py
        pmax = float(row['PMax MW'])
        x = {i: round(float(row[f'Output_pct_{i}'])*pmax, 1)
             for i in range(4) 
        }

        ## /1000. from the RTS-GMLC MATPOWER writer -- 
        ## heat rates are in BTU/kWh, 1BTU == 10^-6 MMBTU, 1kWh == 10^-3 MWh, so MMBTU/MWh == 10^3/10^6 * BTU/kWh
        f = {}
        f[0] = (float(row['HR_avg_0'])*1000./ 1000000.)*x[0]
        for i in range(1,4):
            f[i] = (((x[i]-x[i-1])*(float(row[f'HR_incr_{i}'])*1000. / 1000000.))) + f[i-1]

        fuel_price = float(row['Fuel Price $/MMBTU'])
        y = {i: fuel_price*f[i] for i in range(4)}

        # only include the cost coeffecients that matter
        P_COEFF = [ (x[i], round(y[i],2)) for i in range(4) if (((i == 0) or (x[i-1],y[i-1]) != (x[i], y[i])) and (x[i], y[i]) != (0.,0.)) ]
        if P_COEFF == []:
            P_COEFF = [(pmax, 0.0)] 

        F_COEFF = [ (x[i], round(f[i],2)) for i in range(4) if (((i == 0) or (x[i-1],f[i-1]) != (x[i], f[i])) and (x[i], f[i]) != (0.,0.)) ]
        if F_COEFF == []:
            F_COEFF = [(pmax, 0.0)]
                
        # UC Data
        MIN_DN_TIME = float(row['Min Down Time Hr'])

        # Startup types and costs 
        COLD_HEAT = float(row['Start Heat Cold MBTU'])
        WARM_HEAT = float(row['Start Heat Warm MBTU'])
        HOT_HEAT = float(row['Start Heat Hot MBTU'])

        COLD_TIME = float(row['Start Time Cold Hr'])
        WARM_TIME = float(row['Start Time Warm Hr'])
        HOT_TIME = float(row['Start Time Hot Hr'])

        FIXED_START_COST = float(row['Non Fuel Start Cost $'])

        if (COLD_TIME <= MIN_DN_TIME) or (COLD_TIME == WARM_TIME == HOT_TIME):
            STARTUP_COSTS = [(MIN_DN_TIME, round(COLD_HEAT*fuel_price + FIXED_START_COST, 2))]
            STARTUP_FUEL = [(MIN_DN_TIME, COLD_HEAT)]

        elif WARM_TIME <= MIN_DN_TIME:
            STARTUP_COSTS = [(MIN_DN_TIME, round(WARM_HEAT*fuel_price + FIXED_START_COST, 2)),\
                                (COLD_TIME, round(COLD_HEAT*fuel_price + FIXED_START_COST, 2))]
            STARTUP_FUEL = [(MIN_DN_TIME, WARM_HEAT),\
                                (COLD_TIME, COLD_HEAT)]

        else:
            STARTUP_COSTS = [(MIN_DN_TIME, round(HOT_HEAT*fuel_price+FIXED_START_COST,2)),\
                                (WARM_TIME, round(WARM_HEAT*fuel_price+FIXED_START_COST,2)),\
                                (COLD_TIME, round(COLD_HEAT*fuel_price+FIXED_START_COST,2))]
            STARTUP_FUEL = [(MIN_DN_TIME, HOT_HEAT),\
                                (WARM_TIME, WARM_HEAT),\
                                (COLD_TIME, COLD_HEAT)]
        gen_dict["startup_cost"] = STARTUP_COSTS
        gen_dict["startup_fuel"] = STARTUP_FUEL
        gen_dict["shutdown_cost"] = 0.0

        gen_dict["pc1"] = 0.0
        gen_dict["pc2"] = 0.0
        gen_dict["qc1_min"] = 0.0
        gen_dict["qc1_max"] = 0.0
        gen_dict["qc2_min"] = 0.0
        gen_dict["qc2_max"] = 0.0
        gen_dict["agc_capable"] = True
        gen_dict["p_min_agc"] = gen_dict["p_min"]
        gen_dict["p_max_agc"] = gen_dict["p_max"]

        ramp_q = gen_dict['ramp_q']
        gen_dict["ramp_agc"] = ramp_q
        gen_dict["ramp_10"] = 10.*ramp_q
        gen_dict["ramp_30"] = 30.*ramp_q
        gen_dict["ramp_up_60min"] = 60.*ramp_q
        gen_dict["ramp_down_60min"] = 60.*ramp_q
        
        gen_dict["power_factor"] = 0.0
        gen_dict["p_cost"] = {"data_type": "cost_curve", "cost_curve_type":"piecewise", "values": P_COEFF }
        gen_dict["p_fuel"] = {"data_type": "fuel_curve", "values": F_COEFF }
        gen_dict["fuel_cost"] = fuel_price

        # these assumptions are the same as prescient-rtsgmlc
        gen_dict["startup_capacity"] = gen_dict['p_min']
        gen_dict["shutdown_capacity"]  = gen_dict['p_min']
        gen_dict["min_up_time"] = float(row['Min Up Time Hr'])
        gen_dict["min_down_time"] = MIN_DN_TIME
        gen_dict["must_run"] = False

        elements["generator"][name] = gen_dict

    return model_data

def _compute_bus_load_participation_factors(model_data):
    '''
    compute aggregate load per area, and then compute 
    load participation factors from each bus from that data.

    Returns
    =======
    participation_factors:dict[str,float]
        Maps bus name to the fraction of its area load that it carries (0 to 1)
    '''
    elements = model_data['elements']

    # Sum the loads for each area
    area_total_load = {area:0 for area in elements['area']}
    for name, load in elements["load"].items():
        area = elements["bus"][load["bus"]]["area"]
        area_total_load[area] += load["p_load"]

    bus_load_participation_factors = {}
    for name, load in elements["load"].items():
        area = elements["bus"][load["bus"]]["area"]
        bus_load_participation_factors[name] = load["p_load"] / area_total_load[area]

    return bus_load_participation_factors

def _read_timeseries_data(model_data:dict, rts_gmlc_dir:str,
                          start_time:datetime, end_time:datetime,
                          minutes_per_period:Dict[str,int]):
    """
    Parse all relevant timeseries files

    Returns
    =======
    all_timeseries: DataFrame
        A DataFrame with the following columns:
        [Simulation, Category, Object, Parameter, Scaling Factor, Series]

    The Series column holds the data as a pandas series, indexed by the datetime
    of the value.

    """
    # Where we'll keep our results
    timeseries_data = {'DAY_AHEAD':{}, 'REAL_TIME':{}}

    # All timeseries data that has already been read (map[filename] -> DataFrame)
    timeseries_file_map = {}

    timeseries_pointer_df = pd.read_csv(os.path.join(rts_gmlc_dir, "timeseries_pointers.csv"), header=0)

    elements = model_data['elements']
    params_of_interest = {
        'Generator': { 'PMin MW', 'PMax MW'},
        'Reserve':   {'Requirement'},
        'Area':      {'MW Load'}
        }

    # Add a column to timeseries DF to reference the parsed data 
    # instead of the file name
    timeseries_pointer_df['Series'] = None

    # Store the timeseries data in the timeseries DF
    for idx,row in timeseries_pointer_df.iterrows():
        # Skip rows we don't ingest
        if not row['Category'] in params_of_interest:
            continue
        if not row['Parameter'] in params_of_interest[row['Category']]:
            continue

        # Skip generators not in skeleton
        if row['Category'] == 'Generator' and not row['Object'] in elements['generator']:
            continue
        # Skip areas not in skeleton
        if row['Category'] == 'Area' and not row['Object'] in elements['area']:
            continue

        is_reserve = (row['Category'] == 'Reserve')
        if is_reserve:
            # Skip unrecognized reserve names
            name = str(row['Object'])
            if not is_valid_reserve_name(name, model_data):
                continue

        # Read the timeseries file if we haven't already, using the
        # canonical file path as a key into previously read filenames.
        fname = os.path.abspath(os.path.join(rts_gmlc_dir, row['Data File']))
        if not fname in timeseries_file_map:
            sim = row['Simulation']
            data = _read_timeseries_file(fname, minutes_per_period[sim], 
                                         start_time, end_time, row['Object'])
            timeseries_file_map[fname] = data

        # Save a reference to the relevant data as a Series
        scaling_factor = float(row['Scaling Factor'])
        timeseries_pointer_df.at[idx,'Series'] = timeseries_file_map[fname][row['Object']]*scaling_factor

    # Remove the file path from the DF
    timeseries_pointer_df.pop('Data File')

    # Remove irrelevant rows
    timeseries_pointer_df.dropna(subset=['Series'], inplace=True)

    # Sort by simulation
    timeseries_pointer_df.sort_values(by='Simulation', inplace=True)

    return timeseries_pointer_df

def _get_datetimes_from_strings(begin_time:str, end_time:str):

    datetime_format = "%Y-%m-%d %H:%M:%S"

    datestr = "YYYY-DD-MM"
    midnight = " 00:00:00"

    if isinstance(begin_time,datetime):
        pass
    elif isinstance(begin_time,str):
        if len(begin_time) == len(datestr):
            begin_time += midnight
        begin_time = datetime.strptime(begin_time,datetime_format)
    else:
        raise ValueError("Unable to parse begin_time")

    if isinstance(end_time,datetime):
        pass
    elif isinstance(end_time,str):
        if len(end_time) == len(datestr):
            end_time += midnight
        end_time = datetime.strptime(end_time,datetime_format)
    else:
        raise ValueError("Unable to parse end_time")

    return begin_time, end_time
