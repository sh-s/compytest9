"""
Discomfort Severity of Different Comfort Definitions
Python version used: 3.7
Created by: Shide Salimi, Harvard University, February 2021
ssalimi@gsd.harvard.edu,
shide.salimi@gmail.com,
Version 0.0.0
"""
from pythermalcomfort.psychrometrics import *
import math



def cooling_effect(tdb, tr, vr, rh, met, clo, wme=0, units='SI'):
    """
    Returns the value of the Cooling Effect (`CE`_) calculated in compliance with the ASHRAE 55 2017 Standard [1]_. The `CE`_ of the elevated air speed is the value that,
    when subtracted equally from both the average air temperature and the mean radiant temperature, yields the same `SET`_ under still air as in the first `SET`_ calculation
    under elevated air speed.

    .. _CE: https://en.wikipedia.org/wiki/Thermal_comfort#Cooling_Effect

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] in [°F] if `units` = 'IP'
    tr : float
        mean radiant temperature, default in [°C] in [°F] if `units` = 'IP'
    vr : float
        relative air velocity, default in [m/s] in [fps] if `units` = 'IP'

        Note: vr is the relative air velocity caused by body movement and not the air speed measured by the air velocity sensor.
        It can be calculate using the function :py:meth:`pythermalcomfort.psychrometrics.v_relative`.
    rh : float
        relative humidity, [%]
    met : float
        metabolic rate, [met]
    clo : float
        clothing insulation, [clo]
    wme : float
        external work, [met] default 0
    units: str default="SI"
        select the SI (International System of Units) or the IP (Imperial Units) system.

    Returns
    -------
    ce
        Cooling Effect, default in [°C] in [°F] if `units` = 'IP'

    Examples
    --------
    .. code-block:: python

        >>> from pythermalcomfort.models import cooling_effect
        >>> ce = cooling_effect(tdb=25, tr=25, vr=0.3, rh=50, met=1.2, clo=0.5)
        >>> print(ce)
        1.64

        >>> # for users who wants to use the IP system
        >>> ce = cooling_effect(tdb=77, tr=77, vr=1.64, rh=50, met=1, clo=0.6, units="IP")
        >>> print(ce)
        3.74

    Raises
    ------
    ValueError
        If the cooling effect could not be calculated
    """

    if units.lower() == 'ip':
        tdb, tr, vr = units_converter(tdb=tdb, tr=tr, v=vr)

    still_air_threshold = 0.1

    warnings.simplefilter("ignore")
    
    ce = bisection(lambda x: set_tmp(tdb - x, tr - x, v=still_air_threshold, rh=rh, met=met, clo=clo) - set_tmp(tdb=tdb, tr=tr, v=vr, rh=rh, met=met, clo=clo), 0.0, 15.0, 150)
    if ce is None:
        raise ValueError("It could not calculate the cooling effect")
    warnings.simplefilter("always")

    if units.lower() == 'ip':
        ce = ce / 1.8 * 3.28

    return round(ce, 2)

def pmv_ppd(tdb, tr, vr, rh, met, clo, wme=0, standard='ISO', units='SI'):
    """
    Returns Predicted Mean Vote (`PMV`_) and Predicted Percentage of Dissatisfied (`PPD`_) calculated in accordance to main thermal comfort Standards. The `PMV`_ is an index that
    predicts the mean value of the thermal sensation votes (self-reported perceptions) of a large group of people on a sensation scale expressed from –3 to +3 corresponding to
    the categories \"cold,\" \"cool,\" \"slightly cool,\" \"neutral,\" \"slightly warm,\" \"warm,\" and \"hot.\"[1]_. The `PPD`_ is an index that establishes a quantitative
    prediction of the percentage of thermally dissatisfied people determined from `PMV`_ [1]_.

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] in [°F] if `units` = 'IP'
    tr : float
        mean radiant temperature, default in [°C] in [°F] if `units` = 'IP'
    vr : float
        relative air velocity, default in [m/s] in [fps] if `units` = 'IP'

        Note: vr is the relative air velocity caused by body movement and not the air speed measured by the air velocity sensor.
        It can be calculate using the function :py:meth:`pythermalcomfort.psychrometrics.v_relative`.
    rh : float
        relative humidity, [%]
    met : float
        metabolic rate, [met]
    clo : float
        clothing insulation, [clo]
    wme : float
        external work, [met] default 0
    standard: str (default="ISO")
        comfort standard used for calculation

        - If "ISO", then the ISO Equation is used
        - If "ASHRAE", then the ASHRAE Equation is used

        Note: While the PMV equation is the same for both the ISO and ASHRAE standards,
        the ASHRAE Standard Use of the PMV model is limited to air speeds below 0.20 m/s (40 fpm).
        When air speeds exceed 0.20 m/s (40 fpm), the comfort zone boundaries are adjusted based on the SET model.
    units: str default="SI"
        select the SI (International System of Units) or the IP (Imperial Units) system.

    Returns
    -------
    PMV
        Predicted Mean Vote
    PPD
        Predicted Percentage of Dissatisfied occupants, [%]

    Notes
    -----
    You can use this function to calculate the `PMV`_ and `PPD`_ in accordance with either the ASHRAE 55 2017 Standard [1]_ or the ISO 7730 Standard [2]_.

    .. _PMV: https://en.wikipedia.org/wiki/Thermal_comfort#PMV/PPD_method
    .. _PPD: https://en.wikipedia.org/wiki/Thermal_comfort#PMV/PPD_method

    Examples
    --------
    .. code-block:: python

        >>> from pythermalcomfort.models import pmv_ppd
        >>> results = pmv_ppd(tdb=25, tr=25, vr=0.1, rh=50, met=1.2, clo=0.5, wme=0, standard="ISO")
        >>> print(results)
        {'pmv': 0.08, 'ppd': 5.1}

        >>> print(results['pmv'])
        0.08

        >>> # for users who wants to use the IP system
        >>> results_ip = pmv_ppd(tdb=77, tr=77, vr=0.4, rh=50, met=1.2, clo=0.5, units="IP")
        >>> print(results_ip)
        {'pmv': 0.01, 'ppd': 5.0}

    Raises
    ------
    StopIteration
        Raised if the number of iterations exceeds the threshold
    ValueError
        The 'standard' function input parameter can only be 'ISO' or 'ASHRAE'
    """
    if units.lower() == 'ip':
        tdb, tr, vr = units_converter(tdb=tdb, tr=tr, v=vr)

    standard = standard.lower()
    if standard not in ['iso', 'ashrae']:
        raise ValueError("PMV calculations can only be performed in compliance with ISO or ASHRAE Standards")

    #check_standard_compliance(standard=standard, tdb=tdb, tr=tr, v=vr, rh=rh, met=met, clo=clo)

    # if the relative air velocity is higher than 0.2 then follow methodology ASHRAE Appendix H, H3
    if standard == 'ashrae' and vr >= 0.2:
        # calculate the cooling effect
        ce = cooling_effect(tdb=tdb, tr=tr, vr=vr, rh=rh, met=met, clo=clo, wme=wme)

        tdb = tdb - ce
        tr = tr - ce
        vr = 0.1

    pa = rh * 10 * math.exp(16.6536 - 4030.183 / (tdb + 235))

    icl = 0.155 * clo  # thermal insulation of the clothing in M2K/W
    m = met * 58.15  # metabolic rate in W/M2
    w = wme * 58.15  # external work in W/M2
    mw = m - w  # internal heat production in the human body
    if icl <= 0.078:
        fcl = 1 + (1.29 * icl)
    else:
        fcl = 1.05 + (0.645 * icl)

    # heat transf. coeff. by forced convection
    hcf = 12.1 * math.sqrt(vr)
    taa = tdb + 273
    tra = tr + 273
    tcla = taa + (35.5 - tdb) / (3.5 * icl + 0.1)

    p1 = icl * fcl
    p2 = p1 * 3.96
    p3 = p1 * 100
    p4 = p1 * taa
    p5 = (308.7 - 0.028 * mw) + (p2 * math.pow(tra / 100.0, 4))
    xn = tcla / 100
    xf = tcla / 50
    eps = 0.00015

    n = 0
    while abs(xn - xf) > eps:
        xf = (xf + xn) / 2
        hcn = 2.38 * math.pow(abs(100.0 * xf - taa), 0.25)
        if (hcf > hcn):
            hc = hcf
        else:
            hc = hcn
        xn = (p5 + p4 * hc - p2 * math.pow(xf, 4)) / (100 + p3 * hc)
        n += 1
        if n > 150:
            raise StopIteration('Max iterations exceeded')

    tcl = 100 * xn - 273

    # heat loss diff. through skin
    hl1 = 3.05 * 0.001 * (5733 - (6.99 * mw) - pa)
    # heat loss by sweating
    if mw > 58.15:
        hl2 = 0.42 * (mw - 58.15)
    else:
        hl2 = 0
    # latent respiration heat loss
    hl3 = 1.7 * 0.00001 * m * (5867 - pa)
    # dry respiration heat loss
    hl4 = 0.0014 * m * (34 - tdb)
    # heat loss by radiation
    hl5 = 3.96 * fcl * (math.pow(xn, 4) - math.pow(tra / 100.0, 4))
    # heat loss by convection
    hl6 = fcl * hc * (tcl - tdb)

    ts = 0.303 * math.exp(-0.036 * m) + 0.028
    pmv = ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6)
    ppd = 100.0 - 95.0 * math.exp(-0.03353 * pow(pmv, 4.0) - 0.2179 * pow(pmv, 2.0))

    return {'pmv': round(pmv, 2), 'ppd': round(ppd, 1)}


def pmv(tdb, tr, vr, rh, met, clo, wme=0, standard='ASHRAE', units='SI'):
    """
    Returns Predicted Mean Vote (`PMV`_) calculated in accordance to main thermal comfort Standards. The PMV is an index that predicts the mean value of the thermal sensation votes
    (self-reported perceptions) of a large group of people on a sensation scale expressed from –3 to +3 corresponding to the categories \"cold,\" \"cool,\" \"slightly cool,\"
    \"neutral,\" \"slightly warm,\" \"warm,\" and \"hot.\" [1]_

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] in [°F] if `units` = 'IP'
    tr : float
        mean radiant temperature, default in [°C] in [°F] if `units` = 'IP'
    vr : float
        relative air velocity, default in [m/s] in [fps] if `units` = 'IP'

        Note: vr is the relative air velocity caused by body movement and not the air speed measured by the air velocity sensor.
        It can be calculate using the function :py:meth:`pythermalcomfort.psychrometrics.v_relative`.
    rh : float
        relative humidity, [%]
    met : float
        metabolic rate, [met]
    clo : float
        clothing insulation, [clo]
    wme : float
        external work, [met] default 0
    standard: str (default="ISO")
        comfort standard used for calculation

        - If "ISO", then the ISO Equation is used
        - If "ASHRAE", then the ASHRAE Equation is used

        Note: While the PMV equation is the same for both the ISO and ASHRAE standards,
        the ASHRAE Standard Use of the PMV model is limited to air speeds below 0.20 m/s (40 fpm).
        When air speeds exceed 0.20 m/s (40 fpm), the comfort zone boundaries are adjusted based on the SET model.
        See ASHRAE 55 2017 Appendix H for more information [1]_.
    units: str default="SI"
        select the SI (International System of Units) or the IP (Imperial Units) system.

    Returns
    -------
    PMV : float
        Predicted Mean Vote

    Notes
    -----
    You can use this function to calculate the `PMV`_ [1]_ [2]_.

    .. _PMV: https://en.wikipedia.org/wiki/Thermal_comfort#PMV/PPD_method

    Examples
    --------
    .. code-block:: python

        >>> from pythermalcomfort.models import pmv
        >>> pmv(25, 25, 0.1, 50, 1.2, .5, wme=0)
        0.08
    """

    return pmv_ppd(tdb, tr, vr, rh, met, clo, wme, standard=standard, units=units)['pmv']


def pmv_ppd_sleep(tdb, tr, vr, rh, met, clo, wme=0):
    """
    Returns Predicted Mean Vote (`PMV`_) and Predicted Percentage of Dissatisfied (`PPD`_) for sleep conditions calculated in accordance to the research study
    done by Lin and Deng (2008):
         "A study on the thermal comfort in sleeping environments in the subtropics—Developing a thermal comfort model for sleeping environments"
         https://doi.org/10.1016/j.buildenv.2006.11.026

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] 
    tr : float
        mean radiant temperature, default in [°C] 
    vr : float
        relative air velocity, default in [m/s] 

        Note: vr is the relative air velocity caused by body movement and not the air speed measured by the air velocity sensor.
        It can be calculate using the function :py:meth:`pythermalcomfort.psychrometrics.v_relative`.
    rh : float
        relative humidity, [%]
    met : float
        metabolic rate, [met]
    clo : float
        the total thermal resistance for a bedding system, [clo]
    wme : float
        external work, [met] default 0.    

    Returns
    -------
    PMV
        Predicted Mean Vote
    PPD
        Predicted Percentage of Dissatisfied occupants, [%]

    
    Examples
    --------    
    .. code-block:: python

        >>> from compy.thermaldefinitions import pmv_ppd_sleep
        >>> results = pmv_ppd_sleep(tdb=24, tr=24, vr=0.1, rh=50, met=0.7, clo=2, wme=0)
        >>> print(results)
        {'pmv': -0.2, 'ppd': 5.8}

        >>> print(results['pmv'])
        -0.2        

    """
    # if the relative air velocity is higher than 0.2 
    if vr >= 0.2:
        # calculate the cooling effect
        ce = cooling_effect(tdb=tdb, tr=tr, vr=vr, rh=rh, met=met, clo=clo, wme=wme)

        tdb = tdb - ce
        tr = tr - ce
        vr = 0.1

    pa = rh * 10 * math.exp(16.6536 - 4030.183 / (tdb + 235))
    kpa = pa / 1000

    Rt_ins = 0.155 * clo  # thermal insulation of the clothing in M2K/W
    m = met * 58.15  # metabolic rate in W/M2
    w = wme * 58.15  # external work in W/M2
    mw = m - w  # internal heat production in the human body       
        
    if 0.15 < vr < 1.5:
        hc = 2.7 + 8.7 * math.pow(vr, 0.67)
    elif 0 < vr <= 0.15:
        hc = 5.1        
    
    # Sensible heat loss from skin
    hl1 = 34.6 - ((4.7 * tr +hc * tdb)/(4.7 + hc))
    # Evaporative heat loss
    hl2 = 0.3762 * (5.52 - kpa)
    # Latent respiration heat loss
    hl3 = 0.0173 * m * (5.867 - kpa)
    # Sensible respiration heat loss
    hl4 = 0.0014 * m * (34 - tdb)
    
    # Sensitivity coefficient for evaluating PMV
    ts = 0.303 * math.exp(-0.036 * m) + 0.028
    pmv = ts * (mw - (1/Rt_ins)*(hl1 + hl2) - hl3 - hl4)
    
    ppd = 100.0 - 95.0 * math.exp(-0.03353 * pow(pmv, 4.0) - 0.2179 * pow(pmv, 2.0))

    return {'pmv': round(pmv, 2), 'ppd': round(ppd, 1)}


def pmv_sleep(tdb, tr, vr, rh, met, clo, wme=0):
    """
    Returns Predicted Mean Vote (`PMV`_) for sleep conditions calculated in accordance to the research study
    done by Lin and Deng (2008):
         "A study on the thermal comfort in sleeping environments in the subtropics—Developing a thermal comfort model for sleeping environments"
         https://doi.org/10.1016/j.buildenv.2006.11.026

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] 
    tr : float
        mean radiant temperature, default in [°C] 
    vr : float
        relative air velocity, default in [m/s] 

        Note: vr is the relative air velocity caused by body movement and not the air speed measured by the air velocity sensor.
        It can be calculate using the function :py:meth:`pythermalcomfort.psychrometrics.v_relative`.
    rh : float
        relative humidity, [%]
    met : float
        metabolic rate, [met]
    Rt : float
        the total thermal resistance for a bedding system, [clo]
    wme : float
        external work, [met] default 0
    
    Returns
    -------
    PMV : float
        Predicted Mean Vote

    
    Examples
    --------
    .. code-block:: python

        >>> from compy.thermaldefinitions import pmv_ppd_sleep
        >>> pmv_sleep(24, 24, 0.1, 50, 0.7, 2, wme=0)
        -0.2
    """

    return pmv_ppd_sleep(tdb, tr, vr, rh, met, clo, wme=0)['pmv']


def set_tmp(tdb, tr, v, rh, met, clo, wme=0, body_surface_area=1.8258, patm=101325, units='SI'):
    """
    Calculates the Standard Effective Temperature (SET). The SET is the temperature of an imaginary environment at 50% (rh), <0.1 m/s (20 fpm) average air speed (v), and tr = tdb ,
    in which the total heat loss from the skin of an imaginary occupant with an activity level of 1.0 met and a clothing level of 0.6 clo is the same as that
    from a person in the actual environment with actual clothing and activity level.

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] in [°F] if `units` = 'IP'
    tr : float
        mean radiant temperature, default in [°C] in [°F] if `units` = 'IP'
    v : float
        air velocity, default in [m/s] in [fps] if `units` = 'IP'
    rh : float
        relative humidity, [%]
    met : float
        metabolic rate, [met]
    clo : float
        clothing insulation, [clo]
    wme : float
        external work, [met] default 0
    body_surface_area : float
        body surface area, default value 1.8258 [m2] in [ft2] if `units` = 'IP'
    patm : float
        atmospheric pressure, default value 101325 [Pa] in [atm] if `units` = 'IP'
    units: str default="SI"
        select the SI (International System of Units) or the IP (Imperial Units) system.

    Returns
    -------
    SET : float
        Standard effective temperature, [°C]

    Notes
    -----
    You can use this function to calculate the `SET`_ temperature in accordance with the ASHRAE 55 2017 Standard [1]_.

    .. _SET: https://en.wikipedia.org/wiki/Thermal_comfort#Standard_effective_temperature

    Examples
    --------
    .. code-block:: python

        >>> from pythermalcomfort.models import set_tmp
        >>> set_tmp(tdb=25, tr=25, v=0.1, rh=50, met=1.2, clo=.5)
        25.3

        >>> # for users who wants to use the IP system
        >>> set_tmp(tdb=77, tr=77, v=0.328, rh=50, met=1.2, clo=.5, units='IP')
        77.6

    """
    if units.lower() == 'ip':
        if body_surface_area == 1.8258:
            body_surface_area = 19.65
        if patm == 101325:
            patm = 1
        tdb, tr, v, body_surface_area, patm = units_converter(tdb=tdb, tr=tr, v=v, area=body_surface_area, pressure=patm)

    #check_standard_compliance(standard='ashrae', tdb=tdb, tr=tr, v=v, rh=rh, met=met, clo=clo)

    # Initial variables as defined in the ASHRAE 55-2017
    vapor_pressure = rh * p_sat_torr(tdb) / 100
    air_velocity = max(v, 0.1)
    k_clo = 0.25
    body_weight = 69.9
    met_factor = 58.2
    SBC = 0.000000056697  # Stefan-Boltzmann constant (W/m2K4)
    CSW = 170
    CDIL = 120
    CSTR = 0.5

    temp_skin_neutral = 33.7
    temp_core_neutral = 36.8
    temp_body_neutral = 36.49
    skin_blood_flow_neutral = 6.3

    temp_skin = temp_skin_neutral
    temp_core = temp_core_neutral
    skin_blood_flow = skin_blood_flow_neutral
    ALFA = 0.1
    ESK = 0.1 * met

    pressure_in_atmospheres = patm / 101325
    LTIME = 60
    RCL = 0.155 * clo

    FACL = 1.0 + 0.15 * clo  # INCREASE IN BODY SURFACE AREA DUE TO CLOTHING
    LR = 2.2 / pressure_in_atmospheres
    RM = met * met_factor
    M = met * met_factor

    if clo <= 0:
        WCRIT = 0.38 * pow(air_velocity, -0.29)
        ICL = 1.0
    else:
        WCRIT = 0.59 * pow(air_velocity, -0.08)
        ICL = 0.45

    CHC = 3.0 * pow(pressure_in_atmospheres, 0.53)
    CHCV = 8.600001 * pow((air_velocity * pressure_in_atmospheres), 0.53)
    CHC = max(CHC, CHCV)

    CHR = 4.7
    CTC = CHR + CHC
    RA = 1.0 / (FACL * CTC)
    TOP = (CHR * tr + CHC * tdb) / CTC
    TCL = TOP + (temp_skin - TOP) / (CTC * (RA + RCL))

    TCL_OLD = False
    flag = True
    i = 0
    for TIM in range(LTIME):
        while abs(TCL - TCL_OLD) > 0.01:
            if flag:
                i += 1
                TCL_OLD = TCL
                CHR = 4.0 * SBC * pow(((TCL + tr) / 2.0 + 273.15), 3.0) * 0.72
                CTC = CHR + CHC
                RA = 1.0 / (FACL * CTC)
                TOP = (CHR * tr + CHC * tdb) / CTC
            TCL = (RA * temp_skin + RCL * TOP) / (RA + RCL)
            flag = True
        flag = False
        DRY = (temp_skin - TOP) / (RA + RCL)
        HFCS = (temp_core - temp_skin) * (5.28 + 1.163 * skin_blood_flow)
        ERES = 0.0023 * M * (44.0 - vapor_pressure)
        CRES = 0.0014 * M * (34.0 - tdb)
        SCR = M - HFCS - ERES - CRES - wme
        SSK = HFCS - DRY - ESK
        TCSK = 0.97 * ALFA * body_weight
        TCCR = 0.97 * (1 - ALFA) * body_weight
        DTSK = (SSK * body_surface_area) / (TCSK * 60.0)
        DTCR = SCR * body_surface_area / (TCCR * 60.0)
        temp_skin = temp_skin + DTSK
        temp_core = temp_core + DTCR
        TB = ALFA * temp_skin + (1 - ALFA) * temp_core
        SKSIG = temp_skin - temp_skin_neutral
        WARMS = (SKSIG > 0) * SKSIG
        COLDS = ((-1.0 * SKSIG) > 0) * (-1.0 * SKSIG)
        CRSIG = (temp_core - temp_core_neutral)
        WARMC = (CRSIG > 0) * CRSIG
        COLDC = ((-1.0 * CRSIG) > 0) * (-1.0 * CRSIG)
        BDSIG = TB - temp_body_neutral
        WARMB = (BDSIG > 0) * BDSIG
        skin_blood_flow = (skin_blood_flow_neutral + CDIL * WARMC) / (1 + CSTR * COLDS)
        if skin_blood_flow > 90.0:
            skin_blood_flow = 90.0
        if skin_blood_flow < 0.5:
            skin_blood_flow = 0.5
        REGSW = CSW * WARMB * math.exp(WARMS / 10.7)
        if REGSW > 500.0:
            REGSW = 500.0
        ERSW = 0.68 * REGSW
        REA = 1.0 / (LR * FACL * CHC)
        RECL = RCL / (LR * ICL)
        EMAX = (p_sat_torr(temp_skin) - vapor_pressure) / (REA + RECL)
        PRSW = ERSW / EMAX
        PWET = 0.06 + 0.94 * PRSW
        EDIF = PWET * EMAX - ERSW
        ESK = ERSW + EDIF
        if PWET > WCRIT:
            PWET = WCRIT
            PRSW = WCRIT / 0.94
            ERSW = PRSW * EMAX
            EDIF = 0.06 * (1.0 - PRSW) * EMAX
        if EMAX < 0:
            EDIF = 0
            ERSW = 0
            PWET = WCRIT
        ESK = ERSW + EDIF
        MSHIV = 19.4 * COLDS * COLDC
        M = RM + MSHIV
        ALFA = 0.0417737 + 0.7451833 / (skin_blood_flow + .585417)

    HSK = DRY + ESK
    W = PWET
    PSSK = p_sat_torr(temp_skin)
    CHRS = CHR
    if met < 0.85:
        CHCS = 3.0
    else:
        CHCS = 5.66 * math.pow((met - 0.85), 0.39)
    if CHCS < 3.0:
        CHCS = 3.0
    CTCS = CHCS + CHRS
    RCLOS = 1.52 / ((met - wme / met_factor) + 0.6944) - 0.1835
    RCLS = 0.155 * RCLOS
    FACLS = 1.0 + k_clo * RCLOS
    FCLS = 1.0 / (1.0 + 0.155 * FACLS * CTCS * RCLOS)
    IMS = 0.45
    ICLS = IMS * CHCS / CTCS * (1 - FCLS) / (CHCS / CTCS - FCLS * IMS)
    RAS = 1.0 / (FACLS * CTCS)
    REAS = 1.0 / (LR * FACLS * CHCS)
    RECLS = RCLS / (LR * ICLS)
    HD_S = 1.0 / (RAS + RCLS)
    HE_S = 1.0 / (REAS + RECLS)

    DELTA = .0001
    dx = 100.0
    set_old = round(temp_skin - HSK / HD_S, 2)
    while abs(dx) > .01:
        ERR1 = (HSK - HD_S * (temp_skin - set_old) - W * HE_S * (PSSK - 0.5 * p_sat_torr(set_old)))
        ERR2 = (HSK - HD_S * (temp_skin - (set_old + DELTA)) - W * HE_S * (PSSK - 0.5 * p_sat_torr((set_old + DELTA))))
        set = set_old - DELTA * ERR1 / (ERR2 - ERR1)
        dx = set - set_old
        set_old = set

    if units.lower() == 'ip':
        set = units_converter(tmp=set, from_units='si')[0]

    return round(set, 4)


def adaptive_ashrae(tdb, tr, t_running_mean, v, units='SI'):
    """
    Determines the adaptive thermal comfort based on ASHRAE 55. The adaptive model relates indoor design temperatures or acceptable temperature ranges to outdoor meteorological
    or climatological parameters.

    Parameters
    ----------
    tdb : float
        dry bulb air temperature, default in [°C] in [°F] if `units` = 'IP'
    tr : float
        mean radiant temperature, default in [°C] in [°F] if `units` = 'IP'
    t_running_mean: float
        running mean temperature, default in [°C] in [°C] in [°F] if `units` = 'IP'
    v : float
        air velocity, default in [m/s] in [fps] if `units` = 'IP'
    units: str default="SI"
        select the SI (International System of Units) or the IP (Imperial Units) system.

    Returns
    -------
    tmp_cmf : float
        Comfort temperature a that specific running mean temperature, default in [°C] or in [°F]
    tmp_cmf_80_low : float
        Lower acceptable comfort temperature for 80% occupants, default in [°C] or in [°F]
    tmp_cmf_80_up : float
        Upper acceptable comfort temperature for 80% occupants, default in [°C] or in [°F]
    tmp_cmf_90_low : float
        Lower acceptable comfort temperature for 90% occupants, default in [°C] or in [°F]
    tmp_cmf_90_up : float
        Upper acceptable comfort temperature for 90% occupants, default in [°C] or in [°F]
    acceptability_80 : bol
        Acceptability for 80% occupants
    acceptability_90 : bol
        Acceptability for 90% occupants

    Notes
    -----
    You can use this function to calculate if your conditions are within the `adaptive thermal comfort region`.
    Calculations with comply with the ASHRAE 55 2017 Standard [1]_.

    Examples
    --------
    .. code-block:: python

        >>> from pythermalcomfort.models import adaptive_ashrae
        >>> results = adaptive_ashrae(tdb=25, tr=25, t_running_mean=20, v=0.1)
        >>> print(results)
        {'tmp_cmf': 24.0, 'tmp_cmf_80_low': 20.5, 'tmp_cmf_80_up': 27.5, 'tmp_cmf_90_low': 21.5, 'tmp_cmf_90_up': 26.5, 'acceptability_80': True, 'acceptability_90': False}

        >>> print(results['acceptability_80'])
        True
        # The conditions you entered are considered to be comfortable for by 80% of the occupants

        >>> # for users who wants to use the IP system
        >>> results = adaptive_ashrae(tdb=77, tr=77, t_running_mean=68, v=0.3, units='ip')
        >>> print(results)
        {'tmp_cmf': 75.2, 'tmp_cmf_80_low': 68.9, 'tmp_cmf_80_up': 81.5, 'tmp_cmf_90_low': 70.7, 'tmp_cmf_90_up': 79.7, 'acceptability_80': True, 'acceptability_90': False}

        >>> results = adaptive_ashrae(tdb=25, tr=25, t_running_mean=9, v=0.1)
        ValueError: The running mean is outside the standards applicability limits
        # The adaptive thermal comfort model can only be used
        # if the running mean temperature is higher than 10°C

    Raises
    ------
    ValueError
        Raised if the input are outside the Standard's applicability limits

    """
    if units.lower() == 'ip':
        tdb, tr, t_running_mean, vr = units_converter(tdb=tdb, tr=tr, tmp_running_mean=t_running_mean, v=v)

    #check_standard_compliance(standard='ashrae', tdb=tdb, tr=tr, v=v)

    # Define the variables that will be used throughout the calculation.
    results = dict()

    to = t_o(tdb, tr, v)

    # See if the running mean temperature is between 10 °C and 33.5 °C (the range where the adaptive model is supposed to be used)
    if 10.0 <= t_running_mean <= 33.5:

        cooling_effect = 0
        # calculate cooling effect of elevated air speed when top > 25 degC.
        if v >= 0.6 and to >= 25:
            if v < 0.9:
                cooling_effect = 1.2
            elif v < 1.2:
                cooling_effect = 1.8
            else:
                cooling_effect = 2.2

        # Figure out the relation between comfort and outdoor temperature depending on the level of conditioning.
        t_cmf = 0.31 * t_running_mean + 17.8
        tmp_cmf_80_low = t_cmf - 3.5
        tmp_cmf_90_low = t_cmf - 2.5
        tmp_cmf_80_up = t_cmf + 3.5 + cooling_effect
        tmp_cmf_90_up = t_cmf + 2.5 + cooling_effect

        def acceptability(t_cmf_lower, t_cmf_upper):
            # See if the conditions are comfortable.
            if t_cmf_lower < to < t_cmf_upper:
                return True
            else:
                return False

        acceptability_80 = acceptability(tmp_cmf_80_low, tmp_cmf_80_up)
        acceptability_90 = acceptability(tmp_cmf_90_low, tmp_cmf_90_up)

        if units.lower() == 'ip':
            t_cmf, tmp_cmf_80_low, tmp_cmf_80_up, tmp_cmf_90_low, tmp_cmf_90_up = units_converter(from_units='si', tmp_cmf=t_cmf, tmp_cmf_80_low=tmp_cmf_80_low,
                                                                                                tmp_cmf_80_up=tmp_cmf_80_up, tmp_cmf_90_low=tmp_cmf_90_low,
                                                                                                tmp_cmf_90_up=tmp_cmf_90_up)

        results = {'tmp_cmf': t_cmf, 'tmp_cmf_80_low': tmp_cmf_80_low, 'tmp_cmf_80_up': tmp_cmf_80_up,
                   'tmp_cmf_90_low': tmp_cmf_90_low, 'tmp_cmf_90_up': tmp_cmf_90_up,
                   'acceptability_80': acceptability_80, 'acceptability_90': acceptability_90, }

    else:
        raise ValueError("The running mean is outside the standards applicability limits")

    return results