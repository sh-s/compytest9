"""
Discomfort Severity of Different Comfort Definitions
Python version used: 3.7
Created by: Shide Salimi, Harvard University, February 2021
ssalimi@gsd.harvard.edu,
shide.salimi@gmail.com,
Version 0.0.0
"""

import pandas as pd
import numpy as np
from pythermalcomfort.psychrometrics import *
from compytest9.thermaldefinitions import *
import math
from operator import *
from shapely.geometry import Point, LinearRing, LineString
from shapely.geometry.polygon import Polygon
from datetime import *

BODY_SURFACE_AREA = 1.8258
# ---------------- 0) General Methods ------------------------------------------------------

def day_of_year(c_date):
    """
    This method receives a date and determines day of year by comparng the date with 1st of January.
    So, the output of this method for 1st of January would be 1, for 10th of January would be 10, and 
    so on.

    Parameters
    ----------
    c_date : string
        Date in '%m/%d' format.

    Returns
    -------
    days : integer
        Day of a year.
    """
    parsed = datetime.strptime(c_date, '%m/%d').date().replace(year=date.today().year-1) 
    # Compare to 1st of January
    date_to_compare = datetime.strptime('01/01', '%m/%d').date().replace(year=date.today().year-1)
    
    days = (parsed - date_to_compare).days  
    if days <= 58:
        days += 1
    else:
        days = days
    
    return days


# ---------------- 1) Comfort Zone Determination Helper Methods ----------------------------
 
def ref_comfortzone(df, i):
    """
    This method determines reference comfort zone based on day of the year and time of the day.
    There are four reference zones based on winter and summer comfort zones for awake (9 am - 12 am)
    and sleep (1 am - 8 am) times.

    Parameters
    ----------
    df : dataframe
        Input dataframe.
    df_helper : dataframe
        Input dataframe, which is a copy of the original df with summer and winter clothing
        but not sleep and awake specific clothing levels.
    i : integer
        Time step.

    Returns
    -------
    ref_cz : list of tuples
        Coordinates of four boundary points of the reference zone,
        in the form of tuples (dry bulb temperature, relative humadity).
    """
    # Reference comfort temperatures
    # for winter (awake) clo = 1, met = 1, v = 0.1
    [t_1_w_a, t_2_w_a, t_3_w_a, t_4_w_a] = [22.8, 26.76, 23.7, 20.23]
    # for winter (sleep) clo = 3.73, met = 0.7, v = 0.1
    [t_1_w_s, t_2_w_s, t_3_w_s, t_4_w_s] = [13.29, 18.95, 17.44, 12.21]     
    # for summer (awake) clo = 0.5, met = 1, v = 0.1
    [t_1_s_a, t_2_s_a, t_3_s_a, t_4_s_a] = [25.92, 28.8, 26.21, 23.67]
    # for summer (sleep) clo = 2, met = 0.7, v = 0.1
    [t_1_s_s, t_2_s_s, t_3_s_s, t_4_s_s] = [23.95, 27.03, 25.16, 22.37]

    
    # Determine Day of year using day_of_year method
    day = day_of_year(df['Date'][i])
    # Number of winter and summer days in a year 
    [wdays_1, sdays, wdays_2] = [90, 183, 92]    
       
    if (day <= wdays_1) or ((sdays + wdays_1) < day <= sum([wdays_1, sdays, wdays_2])):
        # for Winter awake time, met = 1, v = 0.1,    
        if int(df['Time'][i][:2]) >= 9:
            ref_cz = [(t_1_w_a, 0), (t_2_w_a, 0), (t_3_w_a, 100), (t_4_w_a, 100)]
        # for Winter sleep time, met = 0.7, v = 0.1
        else:
            ref_cz = [(t_1_w_s, 0), (t_2_w_s, 0), (t_3_w_s, 100), (t_4_w_s, 100)] 
    else:
        # for Summer awake time, met = 1, v = 0.1,    
        if int(df['Time'][i][:2]) >= 9:
            ref_cz = [(t_1_s_a, 0), (t_2_s_a, 0), (t_3_s_a, 100), (t_4_s_a, 100)]
        # for Summer sleep time, met = 0.7, v = 0.1
        else:
            ref_cz = [(t_1_s_s, 0), (t_2_s_s, 0), (t_3_s_s, 100), (t_4_s_s, 100)] 
    return ref_cz


def ref_parameters(df, i):
    """
    This method determines reference comfort parameters (clo_c, met_c, v_c) based on 
    day of the year and time of the day.
    
    Parameters
    ----------
    df : dataframe
        Input dataframe.
    df_helper : dataframe
        Input dataframe, which is a copy of the original df with summer and winter clothing
        but not sleep and awake specific clothing levels.
    i : integer
        Time step.

    Returns
    -------
    ref_param : list
        List of affecting parameters in comfort condition (clo_c, met_c, v_c).
    """      
    # Determine Day of year
    day = day_of_year(df['Date'][i])
    # Number of winter and summer days in a year
    [wdays_1, sdays, wdays_2] = [90, 183, 92]   
       
    if day <= wdays_1 or (sdays + wdays_1) < day <= sum([wdays_1, sdays, wdays_2]):
        # for Winter awake time, met = 1, v = 0.1,    
        if int(df['Time'][i][:2]) >= 9:            
            [clo_c, met_c, v_c] = [1, 1, 0.1]
        # for Winter sleep time, met = ???, v = 0.1
        else:
            [clo_c, met_c, v_c] = [3.73, 0.7, 0.1] 
    else:
        # for Summer awake time, met = 1, v = 0.1,    
        if int(df['Time'][i][:2]) >= 9:
            [clo_c, met_c, v_c] = [0.5, 1, 0.1]
        # for Summer sleep time, met = ???, v = 0.1
        else:
            [clo_c, met_c, v_c] = [2, 0.7, 0.1] 
    ref_param = [clo_c, met_c, v_c]
    return ref_param


def finding_BoundaryPoints(df, i, t_ref, step, rh, comp, threshold, clo, met, v):
    """
    This method determines each boundary point of the comfort zone based on the PMV values that should be 
    between -0.5 and +0.5. It gets an initial temperature and increases/decreases this temperature based on
    the defined "step" to get to the point where the PMV would be bigger/smaller than -0.5/+0.5.

    Parameters
    ----------
    df : dataframe
        Input dataframe.
    i : integer
        Time step.
    t_ref : float
        Initial temperature, [°C] 
    step : float
        To which degree the initail temperarue should decrease or increase.
    rh : float
        relative humidity, [%]
    comp : string
        Type of operator to apply: Perform “rich comparisons” between two values, such as a and b.
            ge(a, b) is equivalent to a >= b.
            le(a, b) is equivalent to a <= b        
    threshold : float
        PMV comfort boundary.
    clo : float
        clothing insulation, [clo]
    met : float
        metabolic rate, [met]
    v : float
        air velocity, default in [m/s] in [fps] if `units` = 'IP'

    Returns
    -------
    t_new : float
        New temperature corresponding comfort condition (PMW = either -0.5 or +0.5), [°C] 

    """
    init_state = t_ref
    counter = True
    while counter:
        init_state = round(init_state + step, 2) 
        vr = v_relative(v, met)  
        if int(df['Time'][i][:2]) >= 9:            
            pmv_c = pmv(init_state, init_state, vr, rh, met, clo, standard='ASHRAE', units='SI')
        else:
            pmv_c = pmv_sleep(init_state, init_state, vr, rh, met, clo)
          
        if comp(pmv_c, threshold):        
            t_new = init_state
            counter = False
    return t_new


def shift_to_right(df, i, t_ref, clo, met, v): 
    """
    This method receives one reference boundary point of the ASHRAE AWAKE comfort zone and determines four new 
    boundary points of the shifted comfort zone to the right based on afecting parameters (clo, met, and v).
    
    Parameters
    ----------
    df : dataframe
        Input dataframe.
    i : integer
        Time step.
    t_ref : float
        Initial temperature, [°C]
    clo : float
        clothing insulation, [clo]
    met : float
        metabolic rate, [met]
    v : float
        air velocity, default in [m/s] in [fps] if `units` = 'IP'

    Returns
    -------
    list of tuples
        Coordinates of four boundary points of shifted comfort zone,
        in the form of tuples (dry bulb temperature, relative humadity).
    """
    # finding four boundary points by moving the ASHRAE comfort boundary to the right
    p_4 = finding_BoundaryPoints(df, i, t_ref, 0.1, 100, ge, -0.5, clo, met, v)
    p_1 = finding_BoundaryPoints(df, i, p_4, 0.1, 0, ge, -0.5, clo, met, v)
    # Consider cases where p1 > p3
    vr = v_relative(v, met)
    
    if int(df['Time'][i][:2]) >= 9: 
        pmv_c = pmv(p_1, p_1, vr, 100, met, clo, standard='ASHRAE', units='SI')
    else:
        pmv_c = pmv_sleep(p_1, p_1, vr, 100, met, clo)        
    
    if pmv_c >= 0.5:
        p_3 = finding_BoundaryPoints(df, i, p_1, -0.1, 100, le, 0.5, clo, met, v)
    else:
        p_3 = finding_BoundaryPoints(df, i, p_1, 0.1, 100, ge, 0.5, clo, met, v)
    p_2 = finding_BoundaryPoints(df, i, p_3, 0.1, 0, ge, 0.5, clo, met, v)
     
    state_t = [p_1, p_2, p_3, p_4]
    return [(p_1, 0), (p_2, 0), (p_3, 100), (p_4, 100)]


def shift_to_left(df, i, t_ref, clo, met, v):  
    """
    This method receives one reference boundary point of the ASHRAE AWAKE comfort zone and determines four new 
    boundary points of the shifted comfort zone to the left based on afecting parameters (clo and met).

    Parameters
    ----------
    df : dataframe
        Input dataframe.
    i : integer
        Time step.
    t_ref : float
        Initial temperature, [°C]
    clo : float
        clothing insulation, [clo]
    met : float
        metabolic rate, [met]
    v : float
        air velocity, default in [m/s] in [fps] if `units` = 'IP'

    Returns
    -------
    list of tuples
        Coordinates of four boundary points of shifted comfort zone,
        in the form of tuples (dry bulb temperature, relative humadity).
    """
    # finding four boundary points by moving the Ref comfort boundary to the left
    p_4 = finding_BoundaryPoints(df, i, t_ref, -0.1, 100, le, -0.5, clo, met, v)
    p_1 = finding_BoundaryPoints(df, i, p_4, 0.1, 0, ge, -0.5, clo, met, v)
    # Consider cases where p1 > p3
    vr = v_relative(v, met)
    
    if int(df['Time'][i][:2]) >= 9: 
        pmv_c = pmv(p_1, p_1, vr, 100, met, clo, standard='ASHRAE', units='SI')   
    else:
        pmv_c = pmv_sleep(p_1, p_1, vr, 100, met, clo)
        
    if pmv_c >= 0.5:
        p_3 = finding_BoundaryPoints(df, i, p_1, -0.1, 100, le, 0.5, clo, met, v)
    else:
        p_3 = finding_BoundaryPoints(df, i, p_1, 0.1, 100, ge, 0.5, clo, met, v)
    p_2 = finding_BoundaryPoints(df, i, p_3, 0.1, 0, ge, 0.5, clo, met, v)
     
    state_t = [p_1, p_2, p_3, p_4]
    return [(p_1, 0), (p_2, 0), (p_3, 100), (p_4, 100)]   

    
def affecting_param(j, in_point, ref): 
    """
    This method receives values of clo, met, and v corresponding to the current point of time (in_point) 
    along with values of clo, met, and v corresponding to the reference comfort zone (ref) and returns 
    values of clo, met, and v based on the j value. For example, if j == 0, it means the affecting
    parameter is clo and its value corresponding to the current point of time should be considered 
    in the calculations while met and v values corresponding to the reference comfort zone are used.
    When j == 1, it means clo is already considered and now the method moves to the second affecting
    parameter (met).

    Parameters
    ----------
    j : integer
        An integer with value 0, 1, or 2 based on the order of affecting parameter.
    in_point : list of values
        Values of clo, met, and v corresponding to the current point of time.
    ref : list of values
        Values of clo, met, and v corresponding to the reference comfort zone.

    Returns
    -------
    clo : float
        clothing insulation, [clo]
    met : float
        metabolic rate, [met]
    v : float
        air velocity, default in [m/s] in [fps] if `units` = 'IP'
    """
    if j == 0:
        [clo, met, v] = [in_point[0], ref[1], ref[2]]
    elif j == 1:
        [clo, met, v] = [in_point[0], in_point[1], ref[2]]
    else:
        [clo, met, v] = [in_point[0], in_point[1], in_point[2]]
    return clo, met, v


def comparison(df, i, comfortzone, t_ref, in_point, ref):
    """
    This method receives current comfort zone (comfortzone) and determines whether the comfort zone 
    should be shifted by comparing the current values (in_point) of the affecting parameters 
    (clo and met) with their corresponding comfort values (ref).

    Parameters
    ----------
    df : dataframe
        Input dataframe.
    i : integer
        Time step.
    comfortzone : list of tuples
        Current comfort zone.
    t_ref : float
        Initial temperature, [°C]
    in_point : list of values
        Values of clo, met, and v corresponding to the current point of time.
    ref : list of values
        Values of clo, met, and v corresponding to the reference comfort zone.

    Returns
    -------
    comfortzone : list of tuples
        New comfort zone determined by shifting the reference comfort zone, if applicable.
    """
    for j in range(len(in_point)-1):
        [clo, met, v] = affecting_param(j, in_point, ref)
        if in_point[j] > ref[j]:            
            comfortzone = shift_to_left(df, i, t_ref, clo, met, v)
            t_ref = comfortzone[3][0]
        elif in_point[j] < ref[j]:
            comfortzone = shift_to_right(df, i, t_ref, clo, met, v)
            t_ref = comfortzone[3][0]
        else:
            comfortzone = comfortzone 
    return comfortzone, t_ref


# Shifting comfort zone due to air speed bigger than 0.2
def comparison_v(df, i, comfortzone, t_ref, in_point, ref):
    """
    This method receives current comfort zone (comfortzone) and determines whether the comfort zone 
    should be shifted by comparing the current values (in_point) of the affecting parameter (v) 
    with 0.2 m/s.

    Parameters
    ----------
    df : dataframe
        Input dataframe.
    i : integer
        Time step.
    comfortzone : list of tuples
        Current comfort zone.
    t_ref : float
        Initial temperature, [°C]
    in_point : list of values
        Values of clo, met, and v corresponding to the current point of time.
    ref : list of values
        Values of clo, met, and v corresponding to the reference comfort zone.

    Returns
    -------
    comfortzone : list of tuples
        New comfort zone determined by shifting the reference comfort zone, if applicable.
    """
    j = 2
    vr = v_relative(in_point[j], in_point[j-1])
    if vr >= 0.2:
        [clo, met, v] = affecting_param(j, in_point, ref)
        comfortzone = shift_to_right(df, i, t_ref, clo, met, v)
        t_ref = comfortzone[3][0]        
    else:
        comfortzone = comfortzone
    return comfortzone


# ---------------- 3) Comfort Zone Determination Methods -------------------------------------
 
def comfort_zone(df, i):
    """
    This method receives input dataframe (df) and for each time step (i) determines corresponding
    comfort zone w/ or w/o shifting.

    Parameters
    ----------
    df : dataframe
        Input dataframe of the processed input data in the "input_data_pythonCode.csv" file.
    df_helper : dataframe
        Input dataframe, which is a copy of the original df with summer and winter clothing
        but not sleep and awake specific clothing levels.
    i : integer
        Time step, [hr]

    Returns
    -------
    comfortzone : list of tuples
        Comfort zone determined by shifting the reference comfort zone, if applicable.
    """
    # Determine ref comfort zane based on time of a day and day of a year
    comfortzone = ref_comfortzone(df, i)
    
    # Shift reference zone if applicable
    t_ref = comfortzone[3][0]
    # if j is defined within a foor loop: j + 1
    # Determine clo, met , and v of point (row of df)
    in_point = [df['clo'][i], df['met'][i], df['v'][i]]
    ref = ref_parameters(df, i)
    # Shifting comfort zone due to changes in clo and met parameters
    comfortzone_c = comparison(df, i, comfortzone, t_ref, in_point, ref)
    comfortzone = comfortzone_c[0]
    t_ref = comfortzone_c[1]
    # Shifting comfort zone due to air speed bigger than 0.2
    comfortzone = comparison_v(df, i, comfortzone, t_ref, in_point, ref)       
    
    return comfortzone


# ---------------- 4) Comfort Model Methods -------------------------------------

def pmv_model(df, i, comfortzone):
    """
    This method determines discomfort severity using the Fanger PMV-based comfort zones. 
    In this method, the discomfort severity is calculated by finding the difference between 
    the SET value of the point in question with that of closest point on the comfort zone boundary.

    Parameters
    ----------
    df : dataframe
        Input dataframe of the processed input data in the "input_data_pythonCode.csv" file.
    i : integer
        Time step, [hr]
    comfortzone : list of tuples
        Output of "latest_comfortzone" method.

    Returns
    -------
    discomfort : float
        Discomfort severity.

    """            
    poly = Polygon(comfortzone)
    # boundary of the polygon = LinearRing
    pol_ext = LinearRing(list(poly.exterior.coords))
    point_t = df['OpTem'][i]
    point_rh = df['rh'][i]
    point = Point(point_t, point_rh)
    # Need to check whether point is inside the poly or on the boundaries of the poly
    if (poly.contains(point) == True) or (pol_ext.contains(point) == True):
        position = 'inside' 
        # The distance will always be 0.0
        discomfort = poly.distance(point)  # 0.0            
    else:
        position = 'outside'
        SET_point = df['SET'][i]
        d = pol_ext.project(point)
        p = pol_ext.interpolate(d)
        closest_point_coords = list(p.coords)[0]
                
        SET_closest_point = set_tmp(list(p.coords)[0][0], list(p.coords)[0][0], df['v'][i], 
                                    list(p.coords)[0][1], df['met'][i], 
                                    df['clo'][i], wme=0, body_surface_area=BODY_SURFACE_AREA, 
                                    patm=101325, units='SI')
                
        discomfort = math.fabs(SET_point - SET_closest_point)
    
    return discomfort


def adaptive_model(df, i, t_rm):
    """
    This method determines discomfort severity using the Adaptive comfort model. 
    In this method, the discomfort severity is calculated by finding the difference 
    between the SET value of the point in question with that of closest
    point on either upper or lower 80% applicability limits.

    Parameters
    ----------
    df : dataframe
        Input dataframe of the processed input data in the "input_data_pythonCode.csv" file.
    i : integer
        Time step, [hr]
    t_rm : float
        Prevailing (running) mean outdoor temperature.

    Returns
    -------
    discomfort : float
        Discomfort severity.

    """        
    # Need to check whether point is acceptble
    if adaptive_ashrae(df['dbTem'][i], df['mrTem'][i], 
                      t_rm, df['v'][i], units="SI")['acceptability_80'] == True:
        discomfort = 0
    else:
        # Determine SET of the point in question                
        SET_point = df['SET'][i]        
        # Upper acceptable comfort temperature for 80% occupants
        tdb_80_up = adaptive_ashrae(df['dbTem'][i], df['mrTem'][i], 
                               t_rm, df['v'][i], units="SI")['tmp_cmf_80_up']
        # Lower acceptable comfort temperature for 80% occupants
        tdb_80_low = adaptive_ashrae(df['dbTem'][i], df['mrTem'][i], 
                               t_rm, df['v'][i], units="SI")['tmp_cmf_80_low']
        if df['OpTem'][i] > tdb_80_up:
            # Determine SET corresponding the upper acceptable comfort temperature for 80% occupants
            SET_80_up = set_tmp(tdb_80_up, tdb_80_up, df['v'][i], 
                            df['rh'][i], df['met'][i], 
                            df['clo'][i], wme=0, body_surface_area=BODY_SURFACE_AREA, 
                            patm=101325, units='SI')
            discomfort = math.fabs(SET_point - SET_80_up)
        else:
            # Determine SET corresponding the lower acceptable comfort temperature for 80% occupants
            SET_80_low = set_tmp(tdb_80_low, tdb_80_low, df['v'][i], 
                            df['rh'][i], df['met'][i], 
                            df['clo'][i], wme=0, body_surface_area=BODY_SURFACE_AREA, 
                            patm=101325, units='SI')
            discomfort = math.fabs(SET_point - SET_80_low)
    
    return discomfort


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        