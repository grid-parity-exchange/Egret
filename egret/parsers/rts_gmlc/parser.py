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
  from typing import Dict, Union, Optional, Tuple, Set

import sys
import os.path
import pandas as pd
from datetime import datetime, timedelta
import dateutil.parser
from collections import namedtuple

from egret.common.log import logger
import egret.data.model_data as md

from .parsed_cache import ParsedCache
from ._reserves import is_valid_reserve_name, reserve_name_map, ScalarReserveData, ScalarReserveValue

def create_ModelData(rts_gmlc_dir:str, 
                     begin_time:Union[datetime,str], end_time:Union[datetime,str], 
                     simulation:str="DAY_AHEAD", t0_state:Optional[dict] = None) -> md.ModelData:

    """
    Create a ModelData object from the RTS-GMLC data.

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to directory holding csv files in RTS-GMLC format (bus.csv, gen.csv, etc).
    begin_time : datetime.datetime or str
        Beginning of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    end_time : datetime.datetime or str
        End of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    simulation : str
        Either "DAY_AHEAD" or "REAL_TIME", which specifies which time series the data is taken from, 
        default is "DAY_AHEAD".
    t0_state : dict or None
        Keys of this dict are thermal generator names, each element of which is another dictionary with
        keys "initial_status", "initial_p_output", and "initial_q_output", which specify whether the
        generator is on at t0, the real power output at t0, and the reactive power output at t0. 
        If t0_state is None, values are read from initial_status.csv in the rts_gmlc_dir.
        If that file does not exist, no initial state data is set in the model.
    
    Returns
    -------
    egret.model_data.ModelData
        Returns a ModelData object with the timeseries data specified
    """
    return md.ModelData(create_model_data_dict(rts_gmlc_dir, begin_time, end_time, simulation, t0_state))

def create_model_data_dict(rts_gmlc_dir:str, 
                           begin_time:Union[datetime,str], end_time:Union[datetime,str],
                           simulation:str="DAY_AHEAD", t0_state:Optional[dict]=None) -> dict:

    """
    Create a model_data dictionary from the RTS-GMLC data.

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to directory holding csv files in RTS-GMLC format (bus.csv, gen.csv, etc).
    begin_time : datetime.datetime or str
        Beginning of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    end_time : datetime.datetime or str
        End of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    simulation : str
        Either "DAY_AHEAD" or "REAL_TIME", which specifies which time series the data is taken from, 
        default is "DAY_AHEAD".
    t0_state : dict or None
        Keys of this dict are thermal generator names, each element of which is another dictionary with
        keys "initial_status", "initial_p_output", and "initial_q_output", which specify whether the
        generator is on at t0, the real power output at t0, and the reactive power output at t0. 
        If this is None, no initial state data is included in the dict.
    
    Returns
    -------
        dict : A dictionary in the format required for the ModelData object.
    """

    # Convert date string to datetimes, if necessary
    begin_time, end_time = _parse_datetimes_if_strings(begin_time, end_time)

    cache = parse_to_cache(rts_gmlc_dir, begin_time, end_time, t0_state)
    model = cache.generate_model(simulation, begin_time, end_time)
    return model.data


def parse_to_cache(rts_gmlc_dir:str, 
                   begin_time:Union[datetime,str],
                   end_time:Union[datetime,str],
                   t0_state:Optional[dict]=None) -> ParsedCache:
    ''' Parse data in RTS-GMLC format, keeping the portions between a start and end time

    rts_gmlc_dir : str
        Path to directory holding csv files in RTS-GMLC format (bus.csv, gen.csv, etc).
    begin_time : datetime.datetime or str
        Beginning of time horizon. 
    end_time : datetime.datetime or str
        End of time horizon.
    simulation : str
        Either "DAY_AHEAD" or "REAL_TIME", which specifies which time series the data is taken from, 
        default is "DAY_AHEAD".
    t0_state : dict or None
        Keys of this dict are thermal generator names, each element of which is another dictionary with
        keys "initial_status", "initial_p_output", and "initial_q_output", which specify whether the
        generator is on at t0, the real power output at t0, and the reactive power output at t0. 
        If this is None, initial state data is not included in the cache.
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

    constant_reserve_data = _get_scalar_reserve_data(rts_gmlc_dir, metadata_df, model_data)

    begin_time, end_time = _parse_datetimes_if_strings(begin_time, end_time)
    # TODO: Validate begin_time and end_time.
    #       Do we want to enforce that they fall within the data date range?

    timeseries_df = _read_timeseries_data(model_data, rts_gmlc_dir,
                           begin_time, end_time, minutes_per_period)

    load_participation_factors = _compute_bus_load_participation_factors(model_data)

    set_t0_data(model_data, rts_gmlc_dir, t0_state)

    return ParsedCache(model_data, begin_time, end_time,
                       minutes_per_period['DAY_AHEAD'], minutes_per_period['REAL_TIME'], 
                       timeseries_df, load_participation_factors, constant_reserve_data)
    
    
def _read_metadata(base_dir:str) -> pd.DataFrame:
    metadata_path = os.path.join(base_dir, "simulation_objects.csv")
    if not os.path.exists(metadata_path):
        raise ValueError(f'RTS-GMLC directory "{rts_gmlc_dir}" does not contain expected CSV files.')

    # Read metadata about the data
    metadata_df = pd.read_csv(metadata_path, index_col=0)

    return metadata_df

def _get_data_date_range(metadata_df) -> Tuple[datetime, datetime]:
    ''' Get the range of dates for which there is data available
    '''

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

def _create_rtsgmlc_skeleton(rts_gmlc_dir:str) -> dict:
    """
    Creates a data dictionary from the RTS-GMLC data files, without loading hourly data

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to directory holding csv files in RTS-GMLC format (bus.csv, gen.csv, etc).
    
    Returns
    -------
    data : dict
        Returns a dict loaded from the RTS-GMLC data
    """
    from math import isnan

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
    has_shunt_cols = 'MW Shunt G' in bus_df and 'MVAR Shunt B' in bus_df
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
            "zone": str(row['Zone']),
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

        if has_shunt_cols:
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
    branch_df = None

    # add the DC branches
    # TODO: see issue #229
    #if os.path.exists(os.path.join(base_dir,'dc_branch.csv')):
    #    elements["dc_branch"] = {}
    #    branch_df = pd.read_csv(os.path.join(base_dir,'dc_branch.csv'))
    #    for idx,row in branch_df.iterrows():

    #        # TODO: The fields below don't match what Egrets expects or supports for DC branches.
    #        #       The code below is just a placeholder.
    #        branch_dict = {
    #            "from_bus": bus_id_to_name[str(row['From Bus'])], 
    #            "to_bus": bus_id_to_name[str(row['To Bus'])],
    #            "in_service": True,
    #            "branch_type": "dc",
    #            "resistance": float(row['R Line'])
    #        }

    #        name = str(row['UID'])
    #        elements["dc_branch"][name] = branch_dict
    #    branch_df = None

    # add the generators
    elements["generator"] = {}
    RENEWABLE_TYPES = {'WIND', 'HYDRO', 'RTPV', 'PV', 'ROR'}
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

        # Remove optional values if not present
        for key in ('p_min', 'p_max', 'q_min', 'q_max', 'ramp_q'):
            if isnan(gen_dict[key]):
                del gen_dict[key]

        UNIT_TYPE = str(row['Unit Type'])
        if UNIT_TYPE in RENEWABLE_TYPES:
            gen_dict["generator_type"] = "renewable"
            # ROR is treated as HYDRO by Egret
            if UNIT_TYPE == 'ROR':
                gen_dict["unit_type"] = "HYDRO"
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
        # There can be any number of 'Output_pct_<i>' columns.
        # Stop at the first one that doesn't exist or doesn't hold a number
        def valid_output_pcts():
            for i in range(50):
                try:
                    val = float(row[f'Output_pct_{i}'])
                    if isnan(val):
                        return
                    yield (i, val)
                except:
                    return
        x = {i: round(val*pmax, 1)
             for i,val in valid_output_pcts()
            }
        fuel_field_count = len(x)

        if fuel_field_count > 0:
            ## /1000. from the RTS-GMLC MATPOWER writer -- 
            ## heat rates are in BTU/kWh, 1BTU == 10^-6 MMBTU, 1kWh == 10^-3 MWh, so MMBTU/MWh == 10^3/10^6 * BTU/kWh
            f = {}
            f[0] = (float(row['HR_avg_0'])*1000./ 1000000.)*x[0]
            for i in range(1,fuel_field_count):
                f[i] = (((x[i]-x[i-1])*(float(row[f'HR_incr_{i}'])*1000. / 1000000.))) + f[i-1]

            F_COEFF = [ (x[i], round(f[i],2)) for i in range(fuel_field_count) if (((i == 0) or (x[i-1],f[i-1]) != (x[i], f[i])) and (x[i], f[i]) != (0.,0.)) ]
            if F_COEFF == []:
                F_COEFF = [(pmax, 0.0)]
            gen_dict["p_fuel"] = {"data_type": "fuel_curve", "values": F_COEFF }
                
        # UC Data
        MIN_DN_TIME = float(row['Min Down Time Hr'])

        # Startup types and costs, from hot to cold
        startup_heat = (float(row['Start Heat Hot MBTU']),
                        float(row['Start Heat Warm MBTU']),
                        float(row['Start Heat Cold MBTU']))
        startup_time = (float(row['Start Time Hot Hr']),
                        float(row['Start Time Warm Hr']),
                        float(row['Start Time Cold Hr']))

        # Arrange fuel requirements from hottest to coldest, ignoring missing values.
        startup_fuel = []
        for i in range(3):
            # Skip blank values
            if isnan(startup_time[i]) or isnan(startup_heat[i]):
                continue

            t = max(startup_time[i], MIN_DN_TIME)
            f = startup_heat[i]

            # For entries with matching times, use to the colder data
            if len(startup_fuel) > 0 and startup_fuel[-1][0] == t:
                startup_fuel[-1] = (t,f)
            else:
                startup_fuel.append((t,f))

        # If the warmest fuel requirement has a time longer than the minimum
        # down time, extend that warmest requirement down to minimum down time.
        if len(startup_fuel) > 0 and startup_fuel[0][0] > MIN_DN_TIME:
            startup_fuel[0] = (MIN_DN_TIME, startup_fuel[0][1])

        gen_dict["startup_fuel"] = startup_fuel
        fixed_startup_cost = float(row['Non Fuel Start Cost $'])
        if not isnan(fixed_startup_cost):
            gen_dict["non_fuel_startup_cost"] = fixed_startup_cost
        gen_dict["shutdown_cost"] = 0.0

        gen_dict["agc_capable"] = True
        gen_dict["p_min_agc"] = gen_dict["p_min"]
        gen_dict["p_max_agc"] = gen_dict["p_max"]

        ramp_q = gen_dict['ramp_q']
        gen_dict["ramp_agc"] = ramp_q
        gen_dict["ramp_up_60min"] = 60.*ramp_q
        gen_dict["ramp_down_60min"] = 60.*ramp_q
        
        gen_dict["fuel_cost"] = float(row['Fuel Price $/MMBTU'])

        # these assumptions are the same as prescient-rtsgmlc
        gen_dict["startup_capacity"] = gen_dict['p_min']
        gen_dict["shutdown_capacity"]  = gen_dict['p_min']
        gen_dict["min_up_time"] = float(row['Min Up Time Hr'])
        gen_dict["min_down_time"] = MIN_DN_TIME

        elements["generator"][name] = gen_dict
    gen_df = None

    return model_data

def _get_scalar_reserve_data(base_dir:str, metadata_df:df, model_dict:dict) -> ScalarReserveData:
    # Store scalar reserve values as stored in the input
    # 
    # Scalar reserve values that apply to both simulation types are stored in the
    # passed in model dict. Scalar values that vary depending on model type are stored
    # in the returned ScalarReserveData.

    da_scalar_reserves, rt_scalar_reserves = _identify_allowed_scalar_reserve_types(metadata_df)
    shared_reserves = da_scalar_reserves.intersection(rt_scalar_reserves)

    # Collect constant scalar reserves
    da_scalars = []
    rt_scalars = []
    reserve_df = pd.read_csv(os.path.join(base_dir,'reserves.csv'))
    system = model_dict['system']
    areas = model_dict['elements']['area']
    for idx,row in reserve_df.iterrows():
        res_name = row['Reserve Product']
        req = float(row['Requirement (MW)'])

        if res_name in reserve_name_map:
            target_dict = system
            area_name = None
        else:
            # reserve name must be <type>_R<area>.
            # split into type and area
            res_name, area_name = res_name.split("_R", 1)
            if res_name not in reserve_name_map:
                logger.warning(f"Skipping reserve for unrecognized reserve type '{res_name}'")
                continue
            if area_name not in areas:
                logger.warning(f"Skipping reserve for unrecognized area '{area_name}'")
                continue
            target_dict = areas[area_name]

        if res_name in shared_reserves:
            # If it applies to both types, save it in the skeleton
            target_dict[reserve_name_map[res_name]] = req
        elif res_name in da_scalar_reserves:
            # If it applies to just day-ahead, save to DA cache
            da_scalars.append(ScalarReserveValue(res_name, area_name, req))
        elif res_name in rt_scalar_reserves:
            # If it applies to just real-time, save to RT cache
            rt_scalars.append(ScalarReserveValue(res_name, area_name, req))

    return ScalarReserveData(da_scalars, rt_scalars)

def _identify_allowed_scalar_reserve_types(metadata_df:df) -> Tuple[Set[str], Set[str]]:
    ''' Return a list of reserve types that apply to each type of model (DA and RT).

    Arguments
    ---------
    metadata_df:df
        The contents of simulation_objects.csv in a DataFrame

    Returns
    -------
    Returns a tuple with two lists of strings, one list for day-ahead models, 
    and another list for real-time models. Each list holds the names of reserve 
    categories whose scalar reserve values (if specified in reserves.csv) should
    be applied to that type of model.

    (day_ahead_reserves, rt_reserves)
        day_ahead_reserves: Sequence[str]
            The names of reserve categories whose scalar values apply to day-ahead models
        rt_reserves: Sequence[str]
            The names of reserve categories whose scalar values apply to real-time models
    '''
    if not 'Reserve_Products' in metadata_df.index:
        # By default, accept all reserve types in both types of model
        all_reserves = set(reserve_name_map.keys())
        return (all_reserves, all_reserves)

    row = metadata_df.loc['Reserve_Products']
    def parse_reserves(which:str) -> Sequence[str]:
        all = row[which]
        if type(all) is not str:
            return {}
        # strip off parentheses, if present
        all = all.strip('()')
        return set(s.strip() for s in all.split(','))

    return (parse_reserves('DAY_AHEAD'), parse_reserves('REAL_TIME'))
    

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
        [Simulation, Category, Object, Parameter, Series]

    The Series column holds the data as a pandas series, indexed by the datetime
    of the value.

    """

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
        timeseries_pointer_df.at[idx,'Series'] = timeseries_file_map[fname][row['Object']]

    # Remove columns that we don't want to preserve
    keepers= {'Simulation', 'Category', 'Object', 'Parameter', 'Series'}
    for c in timeseries_pointer_df.columns:
        if not c in keepers:
            timeseries_pointer_df.pop(c)

    # Remove irrelevant rows
    timeseries_pointer_df.dropna(subset=['Series'], inplace=True)

    # Sort by simulation
    timeseries_pointer_df.sort_values(by='Simulation', inplace=True)

    return timeseries_pointer_df

def _convert_to_datetime(when:Union[datetime,str]):
    '''
    Convert an object to a datetime, if it is not already one.

    Parameters
    ----------
    when: datetime or str
        The date and time to be returned or parsed
    
    Returns
    -------
    datetime
        The passed in object as a datetime, parsing it if necessary

    If `when` is a datetime it is simply returned. 
    If `when` is a string it is parsed, inferring the format
    If `when` is any other type, a TypeError is raised

    '''
    if isinstance(when, datetime):
        return when
    elif isinstance(when, str):
        return dateutil.parser.parse(when)
    else:
        raise TypeError(f'Invalid argument, expected a datetime or str, got a {type(when)}')

def _parse_datetimes_if_strings(begin_time:Union[datetime,str], end_time:Union[datetime,str]):
    '''
    Ensure both dates are datetimes, parsing date strings if necessary.

    Returns
    -------
    begin_time:datetime
        The begin_time as a datetime, parsing it if necessary
    end_time:datetime
        The end_time as a datetime, parsing it if necessary
    '''
                
    begin_time = _convert_to_datetime(begin_time)
    end_time = _convert_to_datetime(end_time)

    return begin_time, end_time

def set_t0_data(md:dict, base_dir:str="", t0_state:Optional[dict]=None):
    """ Put t0 information into the passed in model dict

    Only t0 data for thermal generators is populated.

    Data comes from:
    * t0_state, if provided
    * otherwise, a file called initial_status.csv, if present
    * otherwise, t0 data is left blank

    If t0_state is provided, it should be organized as t0_state[name][value],
    where `name` is the name of a generator, and `value` is 'initial_status',
    'initial_p_output', and 'initial_q_output'.  For any generator included in
    t0_state, all three values must be present.

    If initial_status.csv is used, it must have a header row and may have 
    from 1 to 3 data rows.  Row 1 is 'initial_status'.  Row 2 is 
    'initial_p_output'.  Row 3 is 'initial_q_output'.  Column headers are 
    the generator names. Default values are used for any missing rows.

    Any generators not mentioned in the data source are left untouched.
    """
    if t0_state is not None:
        for name, gen in md['elements']['generator'].items():
            if gen['generator_type']=='thermal' and name in t0_state:
                gen['initial_status'] = t0_state[name]['initial_status']
                gen['initial_p_output'] = t0_state[name]['initial_p_output']
                gen['initial_q_output'] = t0_state[name]['initial_q_output']
        return

    state_fname = os.path.join(base_dir, 'initial_status.csv')
    if os.path.exists(state_fname):
        import csv
        with open(state_fname, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        # We now have a list of rows, from 1 to 3 rows long.
        # Row 1 is 'initial_status', row 2 is 'initial_p_output', and row 3 is 'initial_q_output'.
        # Any missing row uses defaults
        row_count = len(rows)
        for name, gen in md['elements']['generator'].items():
            if gen['generator_type'] != 'thermal':
                continue
            if name not in reader.fieldnames:
                continue
            gen['initial_status'] = float(rows[0][name])
            if gen['initial_status'] < 0:
                gen['initial_p_output'] = 0.0
                gen['initial_q_output'] = 0.0
            else:
                if row_count >= 2:
                    gen['initial_p_output'] = float(rows[1][name])
                else:
                    gen["initial_p_output"] = gen["p_min"]
                if row_count >= 3:
                    gen['initial_q_output'] = float(rows[2][name])
                else:
                    gen["initial_q_output"] = max(0., gen["q_min"])
    else:
        logger.warning("Setting default t0 state in RTS-GMLC parser")
        for name, gen in md['elements']['generator'].items():
            if gen['generator_type']=='thermal':
                gen['initial_status'] = gen['min_up_time']+1
                gen['initial_p_output'] = gen['p_min']
                gen['initial_q_output'] = 0.
