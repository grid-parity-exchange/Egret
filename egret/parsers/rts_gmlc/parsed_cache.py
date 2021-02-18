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

from egret.data.model_data import ModelData

from ._reserves import reserve_name_map

class ParsedCache():

    def __init__(self, model_skeleton:dict, 
                 begin_time:datetime, end_time:datetime,
                 minutes_per_day_ahead_period:int, minutes_per_real_time_period:int,
                 timeseries_data:DataFrame,
                 load_participation_factors:Dict[str,float]):
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


    def generate_model(self, simulation_type:str, begin_time:datetime, end_time:datetime) -> ModelData:
        md = copy.deepcopy(self.skeleton)
        self._process_timeseries_data(md, simulation_type, begin_time, end_time)
        self._insert_system_data(md, simulation_type, begin_time, end_time)
        return ModelData(md)

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

        for bus, load_dict in md['elements']['load'].items():
            # Skip loads from other areas
            if load_dict['area'] != area_name:
                continue

            # Replace skeleton's p_load with the timeseries data, scaled by the load's
            # portion of the area's total load.
            # Also replace q_load, if present, with timeseries
            p_factor = self.load_participation_factors[bus]
            # save skeleton's scalar p_load
            p_load = load_dict['p_load'] if 'p_load' in load_dict else None
            # overwrite skeleton's p_load with timeseries
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
