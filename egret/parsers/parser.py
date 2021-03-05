#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

""" 
This module provides supporting functions for interacting with standard format input data

It includes methods to parse the data and load them into a TemporalGridNetwork object

"""

import os.path
import egret.data.model_data as md
import pandas as pd
import math
from datetime import datetime, timedelta
from collections import namedtuple


def convert_load_by_area_to_source(data_dir, begin_time, end_time, t0_state=None):
    """
    Create a ModelData object from the input data. Assumes data is formatted like the RTS-GMLC repository's 'RTS_Data' directory.

    Parameters
    ----------
    data_dir : str
        Path to data directory
    begin_time : datetime.datetime or str
        Beginning of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    end_time : datetime.datetime or str
        End of time horizon. If str, date/time in "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD" format,
        the later of which assumes a midnight start.
    t0_state : dict or Nonetype
        Keys of this dict are thermal generator names, each element of which is another dictionary with
        keys "initial_status", "initial_p_output", and "initial_q_output", which specify whether the
        generator is on at t0, the real power output at t0, and the reactive power output at t0. 
        If this is None, default values are loaded.
    """
    for simulation in ['DAY_AHEAD', 'REAL_TIME']:
        simulation = simulation.upper()

        base_dir = os.path.join(data_dir, 'SourceData')

        begin_time, end_time = _get_datetimes(begin_time, end_time, base_dir, simulation)

        TimeSeriesPointer = namedtuple('TimeSeriesPointer',
                                       ['Object',
                                        'Simulation',
                                        'Parameter',
                                        'DataFile'])

        DateTimeValue = namedtuple('DateTimeValue',
                                   ['DateTime', 'Value'])

        areas = _get_eligible_areas(rts_gmlc_dir)
        area_names = _get_eligible_area_names(areas)

        Load = namedtuple('Load', ['DateTime'] + area_names)

        timeseries_pointer_df = pd.read_csv(os.path.join(base_dir, "timeseries_pointers.csv"), header=0, sep=',')

        time_delta = end_time - begin_time

        hours = 24 * time_delta.days + math.ceil(time_delta.seconds / 3600.)

        model_data = _create_rtsgmlc_skeleton(rts_gmlc_dir)

        ## create an object for easy iterating
        md_obj = md.ModelData(model_data)

        system = md_obj.data["system"]
        elements = md_obj.data["elements"]

        if simulation == "DAY_AHEAD":
            system["time_period_length_minutes"] = 60
        else:
            system["time_period_length_minutes"] = 5

        # compute aggregate load per area, and then compute
        # load participation factors from each bus from that data.
        region_total_load = {}
        areas = ["Area" + str(i) for i in range(1, 4)]
        for this_region in areas:
            this_region_total_load = 0.0
            ## loads have exactly one bus
            for name, load in md_obj.elements("load"):
                bus = elements["bus"][load["bus"]]
                if bus["area"] == this_region:
                    this_region_total_load += load["p_load"]
            region_total_load[this_region] = this_region_total_load

        bus_load_participation_factor_dict = {}
        bus_Ql_over_Pl_dict = {}
        for name, load in md_obj.elements("load"):
            bus = elements["bus"][load["bus"]]
            bus_load_participation_factor_dict[name] = load["p_load"] / region_total_load[bus["area"]]
            bus_Ql_over_Pl_dict[name] = load["q_load"] / load["p_load"]

        timeseries_pointer_dict = {}
        for timeseries_pointer_index in timeseries_pointer_df.index.tolist():
            this_timeseries_pointer_dict = timeseries_pointer_df.loc[timeseries_pointer_index].to_dict()
            new_timeseries_pointer = TimeSeriesPointer(this_timeseries_pointer_dict["Object"],
                                                       this_timeseries_pointer_dict["Simulation"],
                                                       this_timeseries_pointer_dict["Parameter"],
                                                       os.path.join(base_dir,
                                                                    this_timeseries_pointer_dict["Data File"]))

            timeseries_pointer_dict[
                (new_timeseries_pointer.Object, new_timeseries_pointer.Simulation)] = new_timeseries_pointer

        load_timeseries_spec = timeseries_pointer_dict[("1", simulation)]
        load_timeseries_df = _read_rts_gmlc_table(load_timeseries_spec.DataFile, simulation)
        load_timeseries_df = load_timeseries_df.rename(columns={"Year_Month_Day_Period": "DateTime"})
        start_mask = load_timeseries_df["DateTime"] >= begin_time
        end_mask = load_timeseries_df["DateTime"] < end_time
        masked_load_timeseries_df = load_timeseries_df[start_mask & end_mask]
        load_dict = masked_load_timeseries_df.to_dict(orient='split')
        load_timeseries = []
        for load_row in load_dict["data"]:
            load_timeseries.append(Load(load_row[0],
                                        float(load_row[1]),
                                        float(load_row[2]),
                                        float(load_row[3])))

        times = []
        for load in load_timeseries:
            times.append(str(load.DateTime))

        system["time_keys"] = times

        ## load into grid_network object
        ## First, load Pl, Ql
        for name, load in md_obj.elements("load"):
            pl_dict, ql_dict = dict(), dict()
            bus = elements["bus"][load["bus"]]
            for load_time in load_timeseries:
                area_load = getattr(load_time, bus["area"])
                pl_dict[str(load_time.DateTime)] = round(bus_load_participation_factor_dict[name] * area_load, 2)
                ql_dict[str(load_time.DateTime)] = pl_dict[str(load_time.DateTime)] * bus_Ql_over_Pl_dict[name]
            load["p_load"] = _make_time_series_dict(list(pl_dict.values()))
            load["q_load"] = _make_time_series_dict(list(ql_dict.values()))

        new_load_time_series = []

        day_ahead_load_file = '../timeseries_data_files/Load/new_load_time_series_DA.csv'
        real_time_load_file = '../timeseries_data_files/Load/new_load_time_series_RT.csv'

        for ix, load_time in enumerate(load_timeseries, start=0):
            load_time_series_record = {}
            load_time_series_record['Year'] = load_time.DateTime.year
            load_time_series_record['Month'] = load_time.DateTime.month
            load_time_series_record['Day'] = load_time.DateTime.day

            if simulation == 'DAY_AHEAD':
                load_time_series_record['Period'] = (ix % 24) + 1
            else:
                load_time_series_record['Period'] = (ix % (24 * 12)) + 1

            for name, load in md_obj.elements('load'):
                bus = elements['bus'][load['bus']]
                area_load = getattr(load_time, bus['area'])

                load_time_series_record[name] = round(bus_load_participation_factor_dict[name] * area_load, 2)

            new_load_time_series.append(load_time_series_record)

        new_load_time_series_df = pd.DataFrame(new_load_time_series)
        new_load_time_series_df = new_load_time_series_df[
            ['Year', 'Month', 'Day', 'Period'] + new_load_time_series_df.columns[4:].tolist()]
        new_load_time_series_fname = 'new_load_time_series_{0}.csv'.format('DA' if simulation == "DAY_AHEAD" else 'RT')
        new_load_time_series_df.to_csv(
            os.path.join(data_dir, 'timeseries_data_files', 'Load', new_load_time_series_fname), index=False)

        # Augment time series pointer dataframe.
        for name, load in md_obj.elements('load'):
            new_load_timeseries_spec = {}
            new_load_timeseries_spec['Object'] = name
            new_load_timeseries_spec['Parameter'] = 'Requirement'
            new_load_timeseries_spec['Simulation'] = 'DAY_AHEAD'
            new_load_timeseries_spec['Data File'] = day_ahead_load_file
            timeseries_pointer_df = timeseries_pointer_df.append(new_load_timeseries_spec, ignore_index=True)

            new_load_timeseries_spec = {}
            new_load_timeseries_spec['Object'] = name
            new_load_timeseries_spec['Parameter'] = 'Requirement'
            new_load_timeseries_spec['Simulation'] = 'REAL_TIME'
            new_load_timeseries_spec['Data File'] = real_time_load_file
            timeseries_pointer_df = timeseries_pointer_df.append(new_load_timeseries_spec, ignore_index=True)

        timeseries_pointer_df.loc[timeseries_pointer_df['Object'] != 'Load'].to_csv(
            os.path.join(data_dir, 'SourceData', 'timeseries_pointers.csv'), index=False)


def create_ModelData(data_dir, begin_time, end_time, simulation="DAY_AHEAD", t0_state=None):
    """
    Create a ModelData object from the input data.

    Parameters
    ----------
    data_dir : str
        Path to data directory
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
    return md.ModelData(create_model_data_dict(data_dir, begin_time, end_time, simulation, t0_state))


def create_model_data_dict(rts_gmlc_dir, begin_time, end_time, simulation="DAY_AHEAD", t0_state=None):
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

    simulation = simulation.upper()
    if simulation not in ["DAY_AHEAD", "REAL_TIME"]:
        raise ValueError('simulation must be "DAY_AHEAD" or "REAL_TIME"')

    base_dir = os.path.join(rts_gmlc_dir, 'SourceData')

    begin_time, end_time = _get_datetimes(begin_time, end_time, base_dir, simulation)

    TimeSeriesPointer = namedtuple('TimeSeriesPointer',
                                   ['Object',
                                    'Simulation',
                                    'Parameter',
                                    'DataFile'])

    DateTimeValue = namedtuple('DateTimeValue',
                               ['DateTime', 'Value'])

    areas = _get_eligible_areas(rts_gmlc_dir)
    area_names = _get_eligible_area_names(areas)

    Load = namedtuple('Load', ['DateTime'] + area_names)

    timeseries_pointer_df = pd.read_csv(os.path.join(base_dir, "timeseries_pointers.csv"), header=0, sep=',')

    time_delta = end_time - begin_time

    hours = 24 * time_delta.days + math.ceil(time_delta.seconds / 3600.)

    model_data = _create_rtsgmlc_skeleton(rts_gmlc_dir)

    ## create an object for easy iterating
    md_obj = md.ModelData(model_data)

    system = md_obj.data["system"]
    elements = md_obj.data["elements"]

    if simulation == "DAY_AHEAD":
        system["time_period_length_minutes"] = 60
    else:
        system["time_period_length_minutes"] = 5

    # compute aggregate load per area, and then compute
    # load participation factors from each bus from that data.
    region_total_load = {}
    for this_region in area_names:
        this_region_total_load = 0.0
        ## loads have exactly one bus
        for name, load in md_obj.elements("load"):
            bus = elements["bus"][load["bus"]]
            if bus["area"] == this_region:
                this_region_total_load += load["p_load"]
        region_total_load[this_region] = this_region_total_load

    bus_load_participation_factor_dict = {}
    bus_Ql_over_Pl_dict = {}
    for name, load in md_obj.elements("load"):
        bus = elements["bus"][load["bus"]]
        bus_load_participation_factor_dict[name] = load["p_load"] / region_total_load[bus["area"]]
        bus_Ql_over_Pl_dict[name] = load["q_load"] / load["p_load"]

    timeseries_pointer_dict = {}
    for timeseries_pointer_index in timeseries_pointer_df.index.tolist():
        this_timeseries_pointer_dict = timeseries_pointer_df.loc[timeseries_pointer_index].to_dict()
        new_timeseries_pointer = TimeSeriesPointer(this_timeseries_pointer_dict["Object"],
                                                   this_timeseries_pointer_dict["Simulation"],
                                                   this_timeseries_pointer_dict["Parameter"],
                                                   os.path.join(base_dir, this_timeseries_pointer_dict["Data File"]))

        timeseries_pointer_dict[
            (new_timeseries_pointer.Object, new_timeseries_pointer.Simulation)] = new_timeseries_pointer

    filtered_timeseries = {}
    for name, gen in md_obj.elements("generator", generator_type="renewable"):
        if gen["fuel"] in ["Solar", "Wind", "Hydro"]:
            if (name, simulation) not in timeseries_pointer_dict:
                print("***WARNING - No timeseries pointer entry found for generator=%s" % name)

            else:
                # print("Time series for generator=%s will be loaded from file=%s" % (name, timeseries_pointer_dict[(name,"DAY_AHEAD")].DataFile))
                renewables_timeseries_df = _read_rts_gmlc_table(timeseries_pointer_dict[(name, simulation)].DataFile,
                                                                simulation)
                this_source_timeseries_df = renewables_timeseries_df.loc[:, ["Year_Month_Day_Period", name]]
                this_source_timeseries_df = this_source_timeseries_df.rename(
                    columns={"Year_Month_Day_Period": "DateTime"})

                start_mask = this_source_timeseries_df["DateTime"] >= begin_time
                end_mask = this_source_timeseries_df["DateTime"] < end_time
                this_source_masked_timeseries_df = this_source_timeseries_df[start_mask & end_mask]

                renewables_timeseries_dict = this_source_masked_timeseries_df.to_dict(orient='split')
                renewables_timeseries = []
                for this_row in renewables_timeseries_dict["data"]:
                    renewables_timeseries.append(DateTimeValue(this_row[0],
                                                               float(this_row[1])))
                filtered_timeseries[name] = renewables_timeseries

    for name, load in md_obj.elements("load"):
        load_timeseries_spec = timeseries_pointer_dict[(name, simulation)]
        load_timeseries_df = _read_rts_gmlc_table(load_timeseries_spec.DataFile, simulation)
        load_timeseries_df = load_timeseries_df.rename(columns={"Year_Month_Day_Period": "DateTime"})
        start_mask = load_timeseries_df["DateTime"] >= begin_time
        end_mask = load_timeseries_df["DateTime"] < end_time
        masked_load_timeseries_df = load_timeseries_df[start_mask & end_mask]
        load_dict = masked_load_timeseries_df.to_dict(orient='records')

    reserves_dfs = {}
    spin_reserve_categories = ["Spin_Up_R1", "Spin_Up_R2", "Spin_Up_R3"]

    other_reserve_categories = ["Reg_Down", "Reg_Up", ]
    ## flexiramp products only in day-ahead simulation
    if simulation == "DAY_AHEAD":
        other_reserve_categories += ["Flex_Down", "Flex_Up", ]

    for reserve in spin_reserve_categories:
        reserves_dfs[reserve] = _read_rts_gmlc_table(timeseries_pointer_dict[(reserve, simulation)].DataFile,
                                                     simulation)

    reserves_dict = {}
    for name, reserve_df in reserves_dfs.items():
        reserve_df = reserve_df.rename(columns={"Year_Month_Day_Period": "DateTime"})
        start_mask = reserve_df["DateTime"] >= begin_time
        end_mask = reserve_df["DateTime"] < end_time
        reserve_df = reserve_df[start_mask & end_mask]
        reserve_timeseries = []
        for this_row in reserve_df.to_dict(orient='split')["data"]:
            reserve_timeseries.append(DateTimeValue(this_row[0], float(this_row[1])))
        reserves_dict[name] = reserve_timeseries

    for reserve in other_reserve_categories:
        reserves_dict[reserve] = _read_rts_gmlc_reserve_table(
            timeseries_pointer_dict[(reserve, simulation)].DataFile,
            begin_time,
            end_time,
            simulation,
        )

    times = []
    for load in load_dict:
        times.append(str(load['DateTime']))

    system["time_keys"] = times

    ## load into grid_network object
    ## First, load Pl, Ql
    for name, load in md_obj.elements("load"):
        pl_dict, ql_dict = dict(), dict()
        bus = elements["bus"][load["bus"]]
        for load_row in load_dict:
            pl_dict[str(load_row['DateTime'])] = round(load_row[name], 2)
            ql_dict[str(load_row['DateTime'])] = pl_dict[str(load_row['DateTime'])] * bus_Ql_over_Pl_dict[name]
        load["p_load"] = _make_time_series_dict(list(pl_dict.values()))
        load["q_load"] = _make_time_series_dict(list(ql_dict.values()))

    ## load in area reserve factors
    area_spin_map = _create_rts_gmlc_area_spin_map(rts_gmlc_dir)
    for name, area in md_obj.elements("area"):
        spin_reserve_dict = dict()
        for datetimevalue in reserves_dict[area_spin_map[name]]:
            spin_reserve_dict[str(datetimevalue.DateTime)] = round(datetimevalue.Value, 2)
        area["spinning_reserve_requirement"] = _make_time_series_dict(list(spin_reserve_dict.values()))

    ## load in global reserve factors
    rts_to_egret_reserve_map = {
        "Flex_Down": "flexible_ramp_down_requirement",
        "Flex_Up": "flexible_ramp_up_requirement",
        "Reg_Down": "regulation_down_requirement",
        "Reg_Up": "regulation_up_requirement",
    }
    for reserve in other_reserve_categories:
        system[rts_to_egret_reserve_map[reserve]] = _make_time_series_dict(list(reserves_dict[reserve].values()))

    ## now load renewable generator stuff
    for name, gen in md_obj.elements("generator", generator_type="renewable"):
        if gen["fuel"] not in ["Solar", "Wind", "Hydro"]:
            continue
        renewables_timeseries = filtered_timeseries[name]
        ## for safety, curtailable renewables can go down to 0
        gen["p_min"] = 0.
        output_dict = dict()
        for datetimevalue in renewables_timeseries:
            output_dict[str(datetimevalue.DateTime)] = round(datetimevalue.Value, 2)

        gen["p_max"] = _make_time_series_dict(list(output_dict.values()))
        # set must-take for Hydro and RTPV
        if gen["unit_type"] in ["HYDRO", "RTPV"]:
            ## copy is for safety when overwriting
            gen["p_min"] = _make_time_series_dict(list(output_dict.copy().values()))

    ## get this from the same place the prescient reader does
    if t0_state is None:
        unit_on_time_df = pd.read_csv(os.path.join(base_dir,
                                                   "../FormattedData/PLEXOS/PLEXOS_Solution/DAY_AHEAD Solution Files/noTX/on_time_7.12.csv"),
                                      header=0,
                                      sep=",")
        unit_on_time_df_as_dict = unit_on_time_df.to_dict(orient="split")
        unit_on_t0_state_dict = {}
        for i in range(0, len(unit_on_time_df_as_dict["columns"])):
            gen_id = unit_on_time_df_as_dict["columns"][i]
            unit_on_t0_state_dict[gen_id] = int(unit_on_time_df_as_dict["data"][0][i])

        for name, gen in md_obj.elements("generator", generator_type="thermal"):
            gen["initial_status"] = unit_on_t0_state_dict[name]
            if gen["initial_status"] < 0:
                gen["initial_p_output"] = 0.
                gen["initial_q_output"] = 0.
            else:
                gen["initial_p_output"] = gen["p_min"]
                gen["initial_q_output"] = max(0., gen["q_min"])

    else:
        for name, gen in md_obj.elements("generator", generator_type="thermal"):
            gen["initial_status"] = t0_state[name]["initial_status"]
            gen["initial_p_output"] = t0_state[name]["initial_p_output"]
            gen["initial_q_output"] = t0_state[name]["initial_q_output"]

    return md_obj.data

def _create_rts_gmlc_area_spin_map(rts_gmlc_dir):
    base_dir = os.path.join(rts_gmlc_dir, 'SourceData')
    reserves = pd.read_csv(os.path.join(base_dir, 'reserves.csv'))
    area_spin_map = {}
    areas = _get_eligible_areas(rts_gmlc_dir)
    area_names = _get_eligible_area_names(areas)
    #assuming we have areas that correspond to the "Eligible Regions" category, starting at 1, 2, 3...
    for area, name in zip(areas, area_names):
        spin_name = reserves.loc[reserves['Eligible Regions'] == str(area)]['Reserve Product'].values[0]
        area_spin_map[name] = spin_name
    return area_spin_map

def _get_rts_gmlc_start_end_dates(base_dir, simulation):
    simulation_objects = pd.read_csv(os.path.join(base_dir, 'simulation_objects.csv'))
    date_from = simulation_objects.loc[simulation_objects['Simulation_Parameters'] == 'Date_From']
    date_to = simulation_objects.loc[simulation_objects['Simulation_Parameters'] == 'Date_To']
    from_date_string = ''
    to_date_string = ''
    if simulation == 'DAY_AHEAD':
        from_date_string = date_from.iloc[0]['DAY_AHEAD']
        to_date_string = date_to.iloc[0]['DAY_AHEAD']
    else:
        from_date_string = date_from.iloc[0]['REAL_TIME']
        to_date_string = date_to.iloc[0]['REAL_TIME']
    start_date = datetime.strptime(from_date_string, '%m/%d/%y %H:%M')
    end_date = datetime.strptime(to_date_string, '%m/%d/%y %H:%M')
    return start_date, end_date

def _get_eligible_areas(rts_gmlc_dir):
    base_dir = os.path.join(rts_gmlc_dir, 'SourceData')
    bus = pd.read_csv(os.path.join(base_dir, 'bus.csv'))
    return bus['Area'].drop_duplicates().values.tolist()

def _get_eligible_area_names(areas):
    area_names = list(map(lambda x: 'Area' + str(x), areas))
    return area_names


def _create_rtsgmlc_skeleton(rts_gmlc_dir):
    """
    Creates a grid_data dictionary from the RTS-GMLC data,
    but does not load hourly data

    Parameters
    ----------
    rts_gmlc_dir : str
        Path to RTS-GMLC directory

    Returns
    -------
    grid_data : dict
        Returns a dict loaded from the RTS-GMLC data
    """

    base_dir = os.path.join(rts_gmlc_dir, 'SourceData')

    case_name = "RTS-GMLC"

    model_data = md.ModelData.empty_model_data_dict()

    elements = model_data["elements"]
    system = model_data["system"]

    system["name"] = case_name

    # this is the default used in the MATPOWER writer for RTS-GMLC
    system["baseMVA"] = 100.

    elements["bus"] = {}
    elements["load"] = {}
    elements["shunt"] = {}

    # add the buses
    bus_df = pd.read_csv(os.path.join(base_dir, 'bus.csv'))
    for idx, row in bus_df.iterrows():
        BUS_I = str(row['Bus ID'])
        if row['Bus Type'] == 'PQ':
            BUS_TYPE = 1
        elif row['Bus Type'] == 'PV':
            BUS_TYPE = 2
        elif row['Bus Type'] == 'Ref':
            BUS_TYPE = 3
        else:
            BUS_TYPE = 4

        PD = float(row['MW Load'])
        QD = float(row['MVAR Load'])
        GS = float(row['MW Shunt G'])
        BS = float(row['MVAR Shunt B'])
        BUS_AREA = str(row['Area'])
        VM = float(row['V Mag'])
        VA = float(row['V Angle'])
        BASE_KV = float(row['BaseKV'])
        ZONE = str(int(row['Zone']))
        VMAX = 1.05  # default used in RTS-GMLC MATPOWER writer
        VMIN = 0.95  # default used in RTS-GMLC MATPOWER writer

        bus_dict = dict()

        if BUS_TYPE < 1 or BUS_TYPE > 3:
            raise ValueError(
                "Encountered an unsupported bus type: {} when parsing MATPOWER input file".format(BUS_TYPE))

        bus_types = {1: "PQ", 2: "PV", 3: "ref", 4: "isolated"}
        bus_dict["matpower_bustype"] = bus_types[BUS_TYPE]

        if BUS_TYPE == 3:
            if VA != 0:
                if abs(VA) >= 1e-16:
                    raise ValueError('EGRET only supports reference buses with an angle of 0 degrees.')
                msg = "\nEgret only supports reference buses with an angle of 0 degrees. \nFound a " \
                      "reference bus with an angle close to 0. \n" \
                      "Value: {0}".format(VA) + "\nSetting reference bus angle to 0."
                warnings.warn(msg)
            system["reference_bus"] = BUS_I
            system["reference_bus_angle"] = VA

        if PD != 0 or QD != 0:
            load_dict = {"bus": BUS_I, "in_service": True}
            load_dict["p_load"] = PD
            load_dict["q_load"] = QD
            load_dict["area"] = "Area" + BUS_AREA
            load_dict["zone"] = ZONE
            elements["load"]['load_' + BUS_I] = load_dict

        if GS != 0 or BS != 0:
            shunt_dict = {"shunt_type": "fixed", "bus": BUS_I}
            shunt_dict["gs"] = GS
            shunt_dict["bs"] = BS
            elements["shunt"]['shunt_' + BUS_I] = shunt_dict

        bus_dict["vm"] = VM
        bus_dict["va"] = VA
        if BASE_KV > 0:
            bus_dict["base_kv"] = BASE_KV
        else:
            raise ValueError('BASE_KV value found that is <= 0. Not supported at this time.')

        bus_dict["area"] = "Area" + BUS_AREA
        bus_dict["zone"] = ZONE
        bus_dict["v_min"] = VMIN
        bus_dict["v_max"] = VMAX
        bus_dict["id"] = row['Bus Name']

        elements["bus"][BUS_I] = bus_dict

    # add the areas

    elements["area"] = {}
    areas = _get_eligible_areas(rts_gmlc_dir)
    area_names = _get_eligible_area_names(areas)
    for name in area_names:
        ## TODO: what else should be in here?
        elements["area"][name] = dict()

    elements["branch"] = {}
    # add the branches
    branch_df = pd.read_csv(os.path.join(base_dir, 'branch.csv'))
    for idx, row in branch_df.iterrows():
        name = str(row['UID'])
        F_BUS = str(row['From Bus'])
        T_BUS = str(row['To Bus'])
        BR_R = float(row['R'])
        BR_X = float(row['X'])
        BR_B = float(row['B'])
        RATE_A = float(row['Cont Rating'])
        RATE_B = float(row['Cont Rating'])
        RATE_C = float(row['Cont Rating'])
        if RATE_A == 0:
            RATE_A = None
        if RATE_B == 0:
            RATE_B = None
        if RATE_C == 0:
            RATE_C = None
        TAP = float(row['Tr Ratio'])
        SHIFT = 0.0  # these hard-coded values are the defaults
        BR_STATUS = 1  # from the RTS-GMLC MATPOWER writer
        ANGMIN = -90.
        ANGMAX = 90.
        PF = None  # these values are not given
        QF = None
        PT = None
        QT = None

        branch_dict = {"from_bus": F_BUS, "to_bus": T_BUS}
        branch_dict["resistance"] = BR_R
        branch_dict["reactance"] = BR_X
        branch_dict["charging_susceptance"] = BR_B

        if TAP != 0.0:
            branch_dict["transformer_tap_ratio"] = TAP
            branch_dict["transformer_phase_shift"] = SHIFT
            branch_dict["branch_type"] = "transformer"
        else:
            branch_dict["branch_type"] = "line"

        branch_dict["rating_long_term"] = RATE_A
        branch_dict["rating_short_term"] = RATE_B
        branch_dict["rating_emergency"] = RATE_C
        branch_dict["angle_diff_min"] = ANGMIN
        branch_dict["angle_diff_max"] = ANGMAX
        assert (BR_STATUS == 0 or BR_STATUS == 1)
        if BR_STATUS == 1:
            branch_dict["in_service"] = True
        else:
            branch_dict["in_service"] = False
        branch_dict["pf"] = PF
        branch_dict["qf"] = QF
        branch_dict["pt"] = PT
        branch_dict["qt"] = QT

        elements["branch"][name] = branch_dict

    # add the generators
    elements["generator"] = {}
    RENEWABLE_TYPES = ['WIND', 'HYDRO', 'RTPV', 'PV']
    gen_df = pd.read_csv(os.path.join(base_dir, 'gen.csv'))
    for idx, row in gen_df.iterrows():
        name = str(row['GEN UID'])
        GEN_BUS = str(row['Bus ID'])
        gen_dict = {"bus": GEN_BUS}

        # if this is a renewable, hydro, or storage need to handle differently
        # (hydro schedules in RTS-GMLC are fixed)
        if row['Fuel'] in ['Storage']:
            pass
        else:
            # NOTE: for now, prescient doesn't handle CSP -- not clear how to model
            if row['Unit Type'] == 'CSP':
                continue
            ## (mostly) MATPOWER data
            PG = float(row['MW Inj'])
            QG = float(row['MVAR Inj'])
            QMAX = float(row['QMax MVAR'])
            QMIN = float(row['QMin MVAR'])
            RAMP_Q = 1. * float(row['Ramp Rate MW/Min'])
            VG = float(row['V Setpoint p.u.'])
            MBASE = 100.  # set in RTS-GMLC MATPOWER writer
            GEN_STATUS = 1
            PMAX = float(row['PMax MW'])
            PMIN = float(row['PMin MW'])
            FUEL = str(row['Fuel'])
            UNIT_TYPE = str(row['Unit Type'])

            if UNIT_TYPE in RENEWABLE_TYPES:
                gen_dict["generator_type"] = "renewable"
            elif UNIT_TYPE == 'SYNC_COND':
                ## TODO: should we have a flag for these?
                gen_dict["generator_type"] = "thermal"
            else:
                gen_dict["generator_type"] = "thermal"
            gen_dict["bus"] = GEN_BUS
            gen_dict["mbase"] = MBASE
            gen_dict["in_service"] = True
            gen_dict["pg"] = PG
            gen_dict["qg"] = QG
            gen_dict["vg"] = VG
            gen_dict["p_min"] = PMIN
            gen_dict["p_max"] = PMAX
            gen_dict["q_min"] = QMIN
            gen_dict["q_max"] = QMAX
            gen_dict["ramp_q"] = RAMP_Q
            gen_dict["fuel"] = FUEL
            gen_dict["unit_type"] = UNIT_TYPE
            gen_dict["area"] = elements["bus"][gen_dict["bus"]]["area"]
            gen_dict["zone"] = elements["bus"][gen_dict["bus"]]["zone"]

            # after this is only really needed for thermal units
            if UNIT_TYPE in RENEWABLE_TYPES:
                elements["generator"][name] = gen_dict
                continue

            PC1 = 0.0
            PC2 = 0.0
            QC1MIN = 0.0
            QC1MAX = 0.0
            QC2MIN = 0.0
            QC2MAX = 0.0
            RAMP_AGC = 1. * float(row['Ramp Rate MW/Min'])
            RAMP_10 = 10. * float(row['Ramp Rate MW/Min'])
            RAMP_30 = 30. * float(row['Ramp Rate MW/Min'])
            RAMP_UP_60 = 60. * float(row['Ramp Rate MW/Min'])
            RAMP_DN_60 = 60. * float(row['Ramp Rate MW/Min'])
            APF = 0.0  # 0.0 from RTS-GMLC MATPOWER writer

            # Gen cost
            x = {}
            ## round as in RTS-GMLC Prescient/topysp.py
            x[0] = round(float(row['Output_pct_0']) * float(row['PMax MW']), 1)
            x[1] = round(float(row['Output_pct_1']) * float(row['PMax MW']), 1)
            x[2] = round(float(row['Output_pct_2']) * float(row['PMax MW']), 1)
            x[3] = round(float(row['Output_pct_3']) * float(row['PMax MW']), 1)

            y = {}
            y[0] = float(row['Fuel Price $/MMBTU']) * ((float(row['HR_avg_0']) * 1000. / 1000000.) * x[
                0])  ## /1000. from the RTS-GMLC MATPOWER writer,
            y[1] = float(row['Fuel Price $/MMBTU']) * (((x[1] - x[0]) * (float(row['HR_incr_1']) * 1000. / 1000000.))) + \
                   y[0]
            y[2] = float(row['Fuel Price $/MMBTU']) * (((x[2] - x[1]) * (float(row['HR_incr_2']) * 1000. / 1000000.))) + \
                   y[1]
            y[3] = float(row['Fuel Price $/MMBTU']) * (((x[3] - x[2]) * (float(row['HR_incr_3']) * 1000. / 1000000.))) + \
                   y[2]

            # only include the cost coeffecients that matter
            P_COEFF = [(x[i], round(y[i], 2)) for i in range(4) if
                       (((i == 0) or (x[i - 1], y[i - 1]) != (x[i], y[i])) and (x[i], y[i]) != (0., 0.))]
            if P_COEFF == []:
                P_COEFF = [(PMAX, 0.0)]

                # UC Data
            MIN_UP_TIME = float(row['Min Up Time Hr'])
            MIN_DN_TIME = float(row['Min Down Time Hr'])

            # Startup types and costs
            COLD_HEAT = float(row['Start Heat Cold MBTU'])
            WARM_HEAT = float(row['Start Heat Warm MBTU'])
            HOT_HEAT = float(row['Start Heat Hot MBTU'])

            COLD_TIME = float(row['Start Time Cold Hr'])
            WARM_TIME = float(row['Start Time Warm Hr'])
            HOT_TIME = float(row['Start Time Hot Hr'])

            FUEL_PRICE = float(row['Fuel Price $/MMBTU'])
            FIXED_START_COST = float(row['Non Fuel Start Cost $'])

            if (COLD_TIME <= MIN_DN_TIME) or (COLD_TIME == WARM_TIME == HOT_TIME):
                STARTUP_COSTS = [(MIN_DN_TIME, round(COLD_HEAT * FUEL_PRICE + FIXED_START_COST, 2))]
            elif WARM_TIME <= MIN_DN_TIME:
                STARTUP_COSTS = [(MIN_DN_TIME, round(WARM_HEAT * FUEL_PRICE + FIXED_START_COST, 2)), \
                                 (COLD_TIME, round(COLD_HEAT * FUEL_PRICE + FIXED_START_COST, 2))]
            else:
                STARTUP_COSTS = [(MIN_DN_TIME, round(HOT_HEAT * FUEL_PRICE + FIXED_START_COST, 2)), \
                                 (WARM_TIME, round(WARM_HEAT * FUEL_PRICE + FIXED_START_COST, 2)), \
                                 (COLD_TIME, round(COLD_HEAT * FUEL_PRICE + FIXED_START_COST, 2))]

            SHUTDOWN_COST = 0.0

            gen_dict["pc1"] = PC1
            gen_dict["pc2"] = PC2
            gen_dict["qc1_min"] = QC1MIN
            gen_dict["qc1_max"] = QC1MAX
            gen_dict["qc2_min"] = QC2MIN
            gen_dict["qc2_max"] = QC2MAX
            gen_dict["agc_capable"] = True
            gen_dict["p_min_agc"] = gen_dict["p_min"]
            gen_dict["p_max_agc"] = gen_dict["p_max"]
            gen_dict["ramp_agc"] = RAMP_AGC
            gen_dict["ramp_10"] = RAMP_10
            gen_dict["ramp_30"] = RAMP_30
            gen_dict["ramp_up_60min"] = RAMP_UP_60
            gen_dict["ramp_down_60min"] = RAMP_DN_60
            gen_dict["power_factor"] = APF
            gen_dict["p_cost"] = {"data_type": "cost_curve", "cost_curve_type": "piecewise", "values": P_COEFF}

            gen_dict["startup_cost"] = STARTUP_COSTS
            gen_dict["shutdown_cost"] = SHUTDOWN_COST
            # these assumptions are the same as prescient-rtsgmlc
            gen_dict["startup_capacity"] = PMIN
            gen_dict["shutdown_capacity"] = PMIN
            gen_dict["min_up_time"] = MIN_UP_TIME
            gen_dict["min_down_time"] = MIN_DN_TIME
            gen_dict["must_run"] = False

            elements["generator"][name] = gen_dict

    return model_data


def _read_rts_gmlc_table(file_name, simulation):
    if simulation == "DAY_AHEAD":
        _date_parser = lambda *columns: datetime(*map(int, columns[0:3]), int(columns[3]) - 1)
    else:
        minute_mutli = 5
        hour_divisor = 12
        time_periods_in_day = 24 * hour_divisor
        _date_parser = lambda *columns: datetime(*map(int, columns[0:3]), \
                                                 (int(columns[3]) - 1) // hour_divisor,
                                                 minute_mutli * ((int(columns[3]) - 1) % hour_divisor))
    return pd.read_csv(file_name,
                       header=0,
                       sep=',',
                       parse_dates=[[0, 1, 2, 3]],
                       date_parser=_date_parser)


def _read_rts_gmlc_reserve_table(file_name, begin_time, end_time, simulation):
    table_dict = pd.read_csv(file_name, header=0, sep=',').T.to_dict()

    if simulation == "DAY_AHEAD":
        hour_divisor = 1
        minute_mutli = 0
        time_periods_in_day = 24
    else:
        minute_mutli = 5
        hour_divisor = 12
        time_periods_in_day = 24 * hour_divisor

    by_datetime_dict = dict()
    for day_num, day_data in table_dict.items():
        year = day_data['Year']
        month = day_data['Month']
        day = day_data['Day']
        for i in range(1, time_periods_in_day + 1):
            date_time = datetime(year=int(year), month=int(month), day=int(day),
                                 hour=(i - 1) // hour_divisor, minute=minute_mutli * ((i - 1) % hour_divisor))
            if begin_time <= date_time < end_time:
                by_datetime_dict[str(date_time)] = float(day_data[str(i)])
    return by_datetime_dict


def _make_time_series_dict(values):
    return {"data_type": "time_series", "values": values}


def _get_datetimes(begin_time, end_time, base_dir, simulation):
    datetime_format = "%Y-%m-%d %H:%M:%S"

    datestr = "YYYY-DD-MM"
    midnight = " 00:00:00"

    if isinstance(begin_time, datetime):
        pass
    elif isinstance(begin_time, str):
        if len(begin_time) == len(datestr):
            begin_time += midnight
        begin_time = datetime.strptime(begin_time, datetime_format)
    else:
        raise ValueError("Unable to parse begin_time")

    if isinstance(end_time, datetime):
        pass
    elif isinstance(end_time, str):
        if len(end_time) == len(datestr):
            end_time += midnight
        end_time = datetime.strptime(end_time, datetime_format)
    else:
        raise ValueError("Unable to parse end_time")

    # stay in the times provided
    rts_start_date, rts_end_date = _get_rts_gmlc_start_end_dates(base_dir, simulation)
    assert begin_time >= rts_start_date
    assert end_time <= rts_end_date

    # We only take times in whole hours (for now)
    assert (begin_time.minute == 0. and begin_time.second == 0. and begin_time.microsecond == 0.)
    assert (end_time.minute == 0. and end_time.second == 0. and end_time.microsecond == 0.)

    return begin_time, end_time


if __name__ == '__main__':
    from egret.viz.generate_graphs import generate_stack_graph
    from egret.models.unit_commitment import solve_unit_commitment, create_tight_unit_commitment_model
    import matplotlib.pyplot as plt

    current_dir = os.path.dirname(os.path.abspath(__file__))
    rts_gmlc_dir = os.path.join(current_dir, '..', '..', '..', 'RTS-GMLC',
                                'RTS_Data')  # This is just the root of the RTS-GMLC data set.

    # This converts the load data (in RTS-GMLC format) such that individual loads have their own time series explicitly specified (instead of one system-wide time series).
    # It should only need to be run once.
    convert_load_by_area_to_source(
        rts_gmlc_dir, "2020-01-01", "2020-12-31",
        t0_state=None,
    )

    # Test model creation and UC solve for one day using the newly formatted data.
    begin_time = "2020-07-05"
    end_time = "2020-07-06"

    md = create_ModelData(
        rts_gmlc_dir, begin_time, end_time,
        simulation="DAY_AHEAD",
        t0_state=None,
    )

    solved_md = solve_unit_commitment(md,
                                      'gurobi_persistent',
                                      mipgap=0.001,
                                      timelimit=None,
                                      solver_tee=True,
                                      symbolic_solver_labels=False,
                                      options=None,
                                      uc_model_generator=create_tight_unit_commitment_model,
                                      relaxed=False,
                                      return_model=False
                                      )

    fig, ax = generate_stack_graph(
        solved_md,
        title=begin_time,
        show_individual_components=False,
        plot_individual_generators=False,
        x_tick_frequency=4,
    )

    plt.show()
