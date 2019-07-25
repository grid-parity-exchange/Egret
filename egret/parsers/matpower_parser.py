#  ___________________________________________________________________________
#
#  EGRET: Electrical Grid Research and Engineering Tools
#  Copyright 2019 National Technology & Engineering Solutions of Sandia, LLC
#  (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
#  Government retains certain rights in this software.
#  This software is distributed under the Revised BSD License.
#  ___________________________________________________________________________

"""
This module provides supporting functions for parsing (and writing) MATPOWER input files

.. todo::
    documentation and examples
"""

import os.path
import logging
import warnings
import egret.data.model_data as md
import numpy as np

logger = logging.getLogger('egret.parsers.matpower_parser')

def create_ModelData(matpower_filename):
    """
    Parse a MATPOWER input file into a ModelData object containing
    the model_data dictionary

    Parameters
    ----------
    matpower_filename : str
        Path and filename of the matpower inp file you wish to load

    Returns
    -------
        ModelData
    """
    data = create_model_data_dict(matpower_filename)
    return md.ModelData(data)

def get_create(d, key, default):
    if key not in d:
        d[key] = default
    return d[key]

def create_model_data_dict(matpower_filename):
    """
    Parse a MATPOWER input file into a model_data dictionary

    Parameters
    ----------
    matpower_filename : str
        Path and filename of the MATPOWER inp file you wish to load

    Returns
    -------
        dict : Returns a dictionary in the format required for the ModelData
               object.
    """
    looking_for_keywords = 1
    in_matlab_matrix_defn = 2
    status = looking_for_keywords
    matlab_matrix_name = None

    base_name = os.path.basename(matpower_filename)
    case_name, ext = os.path.splitext(base_name)

    # create the model data object
    model_data = md.ModelData.empty_model_data_dict()
    system = model_data["system"]
    elements = model_data["elements"]
    gen_costs_temp = dict()

    sections_found = dict()
    sections_found['function'] = False
    sections_found['mpc.version'] = False
    sections_found['mpc.baseMVA'] = False
    sections_found['mpc.bus'] = False
    sections_found['mpc.gen'] = False
    sections_found['mpc.branch'] = False
    sections_found['mpc.gen_cost'] = False
    sections_found['mpc.bus_name'] = False
    sections_found['mpc.areas'] = False

    mfile = open(matpower_filename, 'r', encoding="utf-8")
    if not mfile:
        return None

    for line in mfile:
        # remove comments
        comment_idx = line.find('%')
        if comment_idx != -1:
            line = line[:comment_idx]

        # get rid of leading whitespace
        line = line.lstrip()
        line = line.rstrip()

        # if no line left, then read the next
        if len(line) <= 0:
            continue

        # trim any trailing semicolon
        if line[-1] == ';':
            line = line[:-1]
            line.rstrip()

        # tokenize and proceed
        tokens = line.split()
        if len(tokens) == 0:
            # blank line, therefore skip
            continue

        if status == looking_for_keywords:
            if tokens[0] == 'function':
                sections_found['function'] = True
                # expecting "function mpc = <case name>"
                assert(tokens[1] == 'mpc')
                assert(tokens[2] == '=')

                # we expect the case name defined in the file to
                # be the same as the filename (to eliminate
                # confusion and potentially loading the wrong file)
                assert(case_name == tokens[3] and "We expect the case name in the file to match the filename itself.")
                system["model_name"] = case_name

            elif tokens[0] == 'mpc.version':
                sections_found['mpc.version'] = True
                assert(tokens[1] == '=')
                assert(tokens[2] == "'2'")
                system["mpc_version"] = tokens[2]

            elif tokens[0] == 'mpc.baseMVA':
                sections_found['mpc.baseMVA'] = True
                assert(tokens[1] == '=')

                # read the baseMVA value
                system["baseMVA"] = float(tokens[2])

            elif tokens[0] == 'mpc.bus':
                sections_found['mpc.bus'] = True
                # expecting "mpc.bus = ["
                # this is pretty restrictive. It is possible that someone
                # might put data on the first line (e.g. "mpc.bus = [ 1 2 45 ... ;")
                # but for now, this works for all the data files I have
                assert(tokens[1] == '=')
                assert(tokens[2] == '[')
                status = in_matlab_matrix_defn
                matlab_matrix_name = 'bus'

            elif tokens[0] == 'mpc.gen':
                sections_found['mpc.gen'] = True
                # expecting "mpc.gen = ["
                # this is pretty restrictive. It is possible that someone
                # might put data on the first line (e.g. "mpc.gen = [ 1 2 45 ... ;")
                # but for now, this works for all the data files I have
                assert (tokens[1] == '=')
                assert (tokens[2] == '[')
                assert (matlab_matrix_name is None)
                status = in_matlab_matrix_defn
                matlab_matrix_name = 'gen'

            elif tokens[0] == 'mpc.branch':
                sections_found['mpc.branch'] = True
                # expecting "mpc.branch = ["
                # this is pretty restrictive. It is possible that someone
                # might put data on the first line (e.g. "mpc.branch = [ 1 2 45 ... ;")
                # but for now, this works for all the data files I have
                assert (tokens[1] == '=')
                assert (tokens[2] == '[')
                assert (matlab_matrix_name is None)
                status = in_matlab_matrix_defn
                matlab_matrix_name = 'branch'

            elif tokens[0] == 'mpc.gencost':
                sections_found['mpc.gen_cost'] = True
                assert(tokens[1] == '=')
                assert(tokens[2] == '[')
                assert(matlab_matrix_name is None)
                status = in_matlab_matrix_defn
                matlab_matrix_name = 'gencost'

            elif tokens[0].startswith('mpc.'):
                assert matlab_matrix_name is None
                logger.warning('Skipping unknown section: '+str(line))
                warnings.warn('Skipping unknown section: '+str(line))
                status = in_matlab_matrix_defn
                matlab_matrix_name = 'skipping_unknown_section'

            else:
                raise TypeError('Unknown line encountered in matpower file parser: ' + str(line))

        elif status == in_matlab_matrix_defn:
            if tokens[0] == ']' or tokens[0] == '}':
                # found the end of the matrix definition
                matlab_matrix_name = None
                status = looking_for_keywords
                continue

            if matlab_matrix_name == 'skipping_unknown_section':
                continue

            if matlab_matrix_name == 'bus_name':
                logger.warning('"bus_name" encountered in MATPOWER input file, but '
                               'not currently supported by the MATPOWER parser.')
                continue # we don't parse bus names for now

            if matlab_matrix_name == 'areas':
                logger.warning('"areas" section encountered in MATPOWER input file, but '
                               'not currently supported by the MATPOWER parser.')
                continue  # we don't parse areas for now

            data = list()
            n_tokens = len(tokens)
            for i in range(n_tokens):
                data.append(float(tokens[i]))

            if matlab_matrix_name == 'bus':
                BUS_I = tokens[0]
                BUS_TYPE = int(tokens[1])
                PD = float(tokens[2])
                QD = float(tokens[3])
                GS = float(tokens[4])
                BS = float(tokens[5])
                BUS_AREA = float(tokens[6])
                VM = float(tokens[7])
                VA = float(tokens[8])
                BASE_KV = float(tokens[9])
                ZONE = float(tokens[10])
                VMAX = float(tokens[11])
                VMIN = float(tokens[12])

                bus_dict = {}
                load_dict = {}
                shunt_dict = {}
                if BUS_TYPE < 1 or BUS_TYPE > 3:
                    raise ValueError("Encountered an unsupported bus type: {} when parsing MATPOWER input file".format(BUS_TYPE))

                # TODO: decide if these are the names we want to use
                # and document them somewhere
                bus_types = {1: "PQ", 2: "PV", 3: "ref", 4: "isolated"}
                bus_dict['matpower_bustype'] = bus_types[BUS_TYPE]

                if BUS_TYPE == 3: # Reference bus
                    system["reference_bus"] = str(BUS_I)
                    system["reference_bus_angle"] = float(VA)

                if PD != 0 or QD != 0:
                    # MATPOWER only has one load per bus
                    load_dict['in_service'] = True
                    load_dict['p_load'] = PD
                    load_dict['q_load'] = QD
                    load_dict['bus'] = str(BUS_I)

                if GS != 0 or BS != 0:
                    ## MATPOWER shunts are fixed
                    shunt_dict['shunt_type'] = 'fixed'
                    shunt_dict['bs'] = BS
                    shunt_dict['gs'] = GS
                    shunt_dict['bus'] = str(BUS_I)

                bus_dict['area'] = BUS_AREA
                bus_dict['vm'] = VM
                bus_dict['va'] = VA
                if BASE_KV > 0:
                    bus_dict['base_kv'] = BASE_KV

                bus_dict['zone'] = ZONE
                bus_dict['v_min'] = VMIN
                bus_dict['v_max'] = VMAX

                buses = get_create(elements, 'bus', dict())
                buses[str(BUS_I)] = bus_dict
                if load_dict:
                    loads = get_create(elements, 'load', dict())
                    loads['load_'+str(BUS_I)] = load_dict

                if shunt_dict:
                    ## TODO: for MATPOWER should this be created as a "shunt_type":"fixed"?
                    shunts = get_create(elements, 'shunt', dict())
                    shunts['shunt_'+str(BUS_I)] = shunt_dict

            elif matlab_matrix_name == 'gen':

                name = str(len(get_create(elements, 'generator', dict()))+1)
                GEN_BUS = tokens[0]
                PG = float(tokens[1])
                QG = float(tokens[2])
                QMAX = float(tokens[3])
                QMIN = float(tokens[4])
                VG = float(tokens[5])
                MBASE = float(tokens[6])
                GEN_STATUS = int(tokens[7])
                PMAX = float(tokens[8])
                PMIN = float(tokens[9])

                gen_dict = dict()
                gen_dict['bus'] = str(GEN_BUS)
                gen_dict['pg'] = PG  # TODO Name?
                gen_dict['qg'] = QG  # TODO Name?
                gen_dict['vg'] = VG  # TODO Name?
                gen_dict['mbase'] = MBASE  # TODO Name?

                if GEN_STATUS > 0:
                    gen_dict['in_service'] = True
                else:
                    gen_dict['in_service'] = False

                gen_dict['p_min'] = PMIN
                gen_dict['p_max'] = PMAX
                gen_dict['q_min'] = QMIN
                gen_dict['q_max'] = QMAX

                # Assume all generators are of type thermal
                gen_dict['generator_type'] = 'thermal'

                if len(tokens) > 10:
                    PC1 = float(tokens[10])
                    PC2 = float(tokens[11])
                    QC1MIN = float(tokens[12])
                    QC1MAX = float(tokens[13])
                    QC2MIN = float(tokens[14])
                    QC2MAX = float(tokens[15])
                    RAMP_AGC = float(tokens[16])
                    RAMP_10 = float(tokens[17])
                    RAMP_30 = float(tokens[18])
                    RAMP_Q = float(tokens[19])
                    APF = float(tokens[20])

                    gen_dict['pc1'] = PC1
                    gen_dict['pc2'] = PC2
                    gen_dict['qc1_min'] = QC1MIN
                    gen_dict['qc1_max'] = QC1MAX
                    gen_dict['qc2_min'] = QC2MIN
                    gen_dict['qc2_max'] = QC2MAX
                    gen_dict['ramp_agc'] = RAMP_AGC
                    gen_dict['ramp_10'] = RAMP_10
                    gen_dict['ramp_30'] = RAMP_30
                    gen_dict['ramp_q'] = RAMP_Q
                    gen_dict['power_factor'] = APF

                get_create(elements, 'generator', dict())[name] = gen_dict

            elif matlab_matrix_name == 'branch':

                name = str(len(get_create(elements, 'branch', dict()))+1)
                F_BUS = tokens[0]
                T_BUS = tokens[1]
                BR_R = float(tokens[2])
                BR_X = float(tokens[3])
                BR_B = float(tokens[4])
                RATE_A = float(tokens[5])
                RATE_B = float(tokens[6])
                RATE_C = float(tokens[7])
                if RATE_A == 0:
                    RATE_A = None
                if RATE_B == 0:
                    RATE_B = None
                if RATE_C == 0:
                    RATE_C = None
                TAP = float(tokens[8])
                SHIFT = float(tokens[9])
                BR_STATUS = int(tokens[10])
                ANGMIN = float(tokens[11])
                ANGMAX = float(tokens[12])
                PF = None
                QF = None
                PT = None
                QT = None
                if (len(tokens) > 13):
                    PF = float(tokens[13])
                    QF = float(tokens[14])
                    PT = float(tokens[15])
                    QT = float(tokens[16])

                branch_dict = dict()
                branch_dict['from_bus'] = str(F_BUS)
                branch_dict['to_bus'] = str(T_BUS)
                branch_dict['resistance'] = BR_R
                branch_dict['reactance'] = BR_X
                branch_dict['charging_susceptance'] = BR_B

                if TAP != 0.0:
                    branch_dict['transformer_tap_ratio'] = TAP
                    branch_dict['transformer_phase_shift'] = SHIFT
                    branch_dict['branch_type'] = 'transformer'
                else:
                    branch_dict['branch_type'] = 'line'

                branch_dict['rating_long_term'] = None
                branch_dict['rating_short_term'] = None
                branch_dict['rating_emergency'] = None
                if RATE_A != 0.0:
                    branch_dict['rating_long_term'] = RATE_A
                if RATE_B != 0.0:
                    branch_dict['rating_short_term'] = RATE_B
                if RATE_C != 0.0:
                    branch_dict['rating_emergency'] = RATE_C
                assert(BR_STATUS == 0 or BR_STATUS == 1)
                if BR_STATUS == 1:
                    branch_dict['in_service'] = True
                else:
                    branch_dict['in_service'] = False


                if ANGMIN == 0 and ANGMAX == 0:
                    branch_dict['angle_diff_min'] = -360.0
                    branch_dict['angle_diff_max'] = 360.0
                else:
                    branch_dict['angle_diff_min'] = ANGMIN
                    branch_dict['angle_diff_max'] = ANGMAX

                branch_dict['pf'] = PF
                branch_dict['qf'] = QF
                branch_dict['pt'] = PT
                branch_dict['qt'] = QT

                get_create(elements, 'branch', dict())[name] = branch_dict

            elif matlab_matrix_name == 'gencost':

                gen_cost_dict = dict()
                gen_idx = len(gen_costs_temp)+1
                MODEL = int(tokens[0])
                if MODEL == 1:
                    STARTUP = float(tokens[1])
                    SHUTDOWN = float(tokens[2])
                    NCOST = int(tokens[3])

                    points = list()
                    cost = list()
                    for i in range(NCOST):
                        points.append(float(tokens[4+2*i]))
                        cost.append(float(tokens[4+2*i+1]))

                    gen_cost_dict['startup_cost'] = STARTUP
                    gen_cost_dict['shutdown_cost'] = SHUTDOWN
                    gen_cost_dict['cost_curve'] = {'data_type': 'cost_curve', 'cost_curve_type': 'piecewise',
                                                   'values': list(zip(points, cost))}

                    gen_costs_temp[gen_idx] = gen_cost_dict

                elif MODEL == 2:

                    STARTUP = float(tokens[1])
                    SHUTDOWN = float(tokens[2])
                    NCOST = int(tokens[3])
                    coeffs = list()
                    for i in range(NCOST):
                        coeffs.append(float(tokens[4+i]))

                    ## this is for the 'values' logic below, MATPOWER does (n-1),...,0, we
                    ## want 0,..,(n-1)
                    coeffs.reverse()

                    gen_cost_dict['startup_cost'] = STARTUP
                    gen_cost_dict['shutdown_cost'] = SHUTDOWN
                    gen_cost_dict['cost_curve'] = {'data_type': 'cost_curve', 'cost_curve_type': 'polynomial',
                                                   'values': {i: coeffs[i] for i in range(len(coeffs))}}
                    gen_costs_temp[gen_idx] = gen_cost_dict

                else:
                    raise ValueError("Found unknown MATPOWER cost model {}".format(MODEL))

            else:
                # these should be the only types we have set above
                assert(False)

    mfile.close()

    if sections_found['function'] != True:
        raise ValueError('Expected a "function" declaration in MATPOWER file')
    if sections_found['mpc.version'] != True:
        raise ValueError('Expected a "mpc.version" declaration in MATPOWER file')
    if sections_found['mpc.baseMVA'] != True:
        raise ValueError('Expected a "mpc.baseMVA" declaration in MATPOWER file')
    if sections_found['mpc.bus'] != True:
        raise ValueError('Expected a "mpc.bus" declaration in MATPOWER file')
    if sections_found['mpc.gen'] != True:
        raise ValueError('Expected a "mpc.gen" declaration in MATPOWER file')
    if sections_found['mpc.branch'] != True:
        raise ValueError('Expected a "mpc.branch" declaration in MATPOWER file')
    if sections_found['mpc.gen_cost'] != True:
        raise ValueError('Expected a "mpc.gen_cost" declaration in MATPOWER file')
    # not concerned if mpc.bus_name does not exist

    # set the generator costs
    # NOTE: as-is, this assumes the "actual" startup/shutdown costs are a part of the
    #       real power cost curve
    generators = get_create(elements, 'generator', dict())
    for gen_idx in gen_costs_temp:
        if gen_idx <= len(generators):
            generators[str(gen_idx)]['p_cost'] = gen_costs_temp[gen_idx]['cost_curve']
            generators[str(gen_idx)]['startup_cost'] = gen_costs_temp[gen_idx]['startup_cost']
            generators[str(gen_idx)]['shutdown_cost'] = gen_costs_temp[gen_idx]['shutdown_cost']
        else:
            generators[str(gen_idx - len(generators))]['q_cost'] = gen_costs_temp[gen_idx]['cost_curve']

    return model_data

