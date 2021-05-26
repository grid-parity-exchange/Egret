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
  from typing import Dict
  from pandas import DataFrame
  from datetime import datetime

import copy
from datetime import timedelta

from egret.data.model_data import ModelData

from ._reserves import reserve_name_map, ScalarReserveData

class ParsedCache():

    def __init__(self, model_skeleton:dict, 
                 begin_time:datetime, end_time:datetime,
                 minutes_per_day_ahead_period:int, minutes_per_real_time_period:int,
                 timeseries_data:DataFrame,
                 load_participation_factors:Dict[str,float],
                 scalar_reserve_data:ScalarReserveData):
        self.skeleton = model_skeleton
        self.begin_time = begin_time
        self.end_time = end_time
        self.timeseries_df = timeseries_data
        self.minutes_per_period = {
            'DAY_AHEAD': minutes_per_day_ahead_period,
            'REAL_TIME': minutes_per_real_time_period
        }
        self.load_participation_factors = load_participation_factors

        # Find and save the index of the first row of each sim type in timeseries_df
        cur_sim = self.timeseries_df['Simulation'][0]
        self._first_indices = {cur_sim:0}
        for i in range(1,len(self.timeseries_df)):
            if self.timeseries_df['Simulation'].iat[i] != cur_sim:
                cur_sim = self.timeseries_df['Simulation'].iat[i]
                self._first_indices[cur_sim] = i

        self.scalar_reserve_data = scalar_reserve_data


    def generate_model(self, simulation_type:str, begin_time:datetime, end_time:datetime) -> ModelData:
        """ Create a new model populated with requested data 

        Parameters
        ----------
        simulation_type:str
            Either 'DAY_AHEAD' or 'REAL_TIME'
        begin_time:datetime
            The earliest time to include in the returned data
        end_time:datetime
            The earliest time to NOT include in the returned data
        """
        md = self.get_new_skeleton()
        self.populate_skeleton_with_data(md, simulation_type, begin_time, end_time)
        return ModelData(md)

    def get_new_skeleton(self) -> dict:
        """ Get a new model dict with system elements but no time-specific data
        """
        return copy.deepcopy(self.skeleton)

    def populate_skeleton_with_data(self, skeleton_dict:dict, simulation_type:str, 
                                    begin_time:datetime, end_time:datetime) -> None:
        """ Update an existing model dict with requested data

        Parameters
        ----------
        skeleton_dict:dict
            The skeleton model dict to populate with data
        simulation_type:str
            Either 'DAY_AHEAD' or 'REAL_TIME'
        begin_time:datetime
            The earliest time to include in the returned data
        end_time:datetime
            The earliest time to NOT include in the returned data
        """

        #Because pandas includes the end of a range, reduce our end time by one second
        end_time = end_time - timedelta(seconds=1)
        self._insert_scalar_reserve_data(skeleton_dict, simulation_type)
        self._process_timeseries_data(skeleton_dict, simulation_type, begin_time, end_time)
        self._insert_system_data(skeleton_dict, simulation_type, begin_time, end_time)

    def _process_timeseries_data(self, md:dict, simulation_type:str, 
                                 begin_time:datetime, end_time:datetime) -> None:
        df = self.timeseries_df

        # Go through each timeseries value for this simulation type
        for i in range(self._first_indices[simulation_type], len(df)):
            if df.iat[i, df.columns.get_loc('Simulation')] != simulation_type:
                break

            category = df.iat[i, df.columns.get_loc('Category')]

            if category == 'Generator':
                self._process_generator_timeseries(md, begin_time, end_time, i)
            elif category == 'Area':
                self._process_area_timeseries(md, begin_time, end_time, i)
            elif category == 'Reserve':
                self._process_reserve_timeseries(md, begin_time, end_time, i)

    def _process_generator_timeseries(self, md:dict, begin_time:datetime, 
                                      end_time:datetime, df_index:int):
        df = self.timeseries_df
        i = df_index
        gen_name = df.iat[i, df.columns.get_loc('Object')]
        gen_dict = md['elements']['generator'][gen_name]
        param = df.iat[i, df.columns.get_loc('Parameter')]
        data = df.iat[i, df.columns.get_loc('Series')][begin_time:end_time].to_list()

        if param == 'PMin MW':
            gen_dict['p_min'] = { 'data_type': 'time_series',
                                  'values' : data }
        elif param == 'PMax MW':
            gen_dict['p_max'] = { 'data_type': 'time_series',
                                  'values' : data }
        else:
            raise ValueError(f"Unexpected generator timeseries data: {param}")

    def _process_area_timeseries(self, md:dict, begin_time:datetime, 
                                 end_time:datetime, df_index:int):
        df = self.timeseries_df
        i = df_index

        area_name = df.iat[i, df.columns.get_loc('Object')]
        param = df.iat[i, df.columns.get_loc('Parameter')]
        assert(param == "MW Load")
        data = df.iat[i, df.columns.get_loc('Series')][begin_time:end_time]

        skeleton_loads = self.skeleton['elements']['load']
        for bus, load_dict in md['elements']['load'].items():
            # Skip loads from other areas
            if load_dict['area'] != area_name:
                continue

            # Replace skeleton's p_load with the timeseries data, scaled by the load's
            # portion of the area's total load.
            # Also replace q_load, if present, with timeseries
            p_factor = self.load_participation_factors[bus]
            # save skeleton's scalar p_load
            p_load = skeleton_loads[bus]['p_load'] if 'p_load' in skeleton_loads[bus] else None
            # overwrite p_load with timeseries
            load_dict['p_load'] = { 'data_type': 'time_series',
                                    'values' : [v*p_factor for v in data] }
            if p_load is not None and 'q_load' in load_dict:
                q_over_p = load_dict['q_load'] / p_load
                load_dict['q_load'] = { 'data_type': 'time_series',
                                        'values' : [v*q_over_p for v in load_dict['p_load']['values']] }

    def _process_reserve_timeseries(self, md:dict, begin_time:datetime, 
                                    end_time:datetime, df_index:int):
        df = self.timeseries_df
        i = df_index

        res_name = df.iat[i, df.columns.get_loc('Object')]

        if res_name in reserve_name_map:
            target_dict = md['system']
        else:
            # reserve name must be <type>_R<area>, 
            # split into type and area
            res_name, area_name = res_name.split("_R", 1)
            target_dict = md['elements']['area'][area_name]

        data = df.iat[i, df.columns.get_loc('Series')][begin_time:end_time]
        target_dict[reserve_name_map[res_name]] = { 'data_type': 'time_series',
                                                    'values' : data.to_list() }

    def _insert_system_data(self, md:dict, simulation_type:str, 
                            begin_time:datetime, end_time:datetime):
        md['system']['time_period_length_minutes'] = self.minutes_per_period[simulation_type]
        
        df = self.timeseries_df
        sample_df = df.iat[self._first_indices[simulation_type], df.columns.get_loc('Series')]
        dates = sample_df[begin_time:end_time].index
        md['system']['time_keys'] = [dt.strftime('%Y-%m-%d %H:%M') for dt in dates]
        
    def _insert_scalar_reserve_data(self, md:dict, simulation_type:str):
        ''' Insert scalar reserve values into the model dict
        '''
        system = md['system']
        areas = md['elements']['area']

        reserve_list = self.scalar_reserve_data.get_simulation_reserves(simulation_type)
        for res in reserve_list:
            if res.area_name is None:
                target_dict = system
            else:
                target_dict = areas[res.area_name]
            target_dict[reserve_name_map[res.reserve_type]] = res.value
