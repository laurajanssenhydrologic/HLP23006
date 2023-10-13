import numpy as np


## Branches
class Branches:
    bodembreedte = None
    bodemhoogte_benedenstrooms = None
    bodemhoogte_bovenstrooms = None
    hoogte_insteek_linkerzijde = None
    hoogte_insteek_rechterzijde = None
    ruwheidhoog = 23
    ruwheidlaag = 23
    taludhelling_linkerzijde = 0
    taludhelling_rechterzijde = 0
    typeruwheid = "Bos en Bijkerk"
    water_width_index = bodembreedte


## Bridges
class Bridges:
    intreeverlies = 0.5
    lengte = np.nan
    ruwheid = 75.0
    typeruwheid = "StricklerKn"
    uittreeverlies = 0.7


## Culverts
class Culverts:
    breedteopening = np.nan
    gesloten = "yes"
    hoogtebinnenonderkantbene = np.nan#0
    hoogtebinnenonderkantbov = np.nan#0
    hoogteopening = np.nan
    intreeverlies = 0.6
    lengte = 1
    ruwheid = 75.0
    typeruwheid = "StricklerKn"
    uittreeverlies = 0.8
    vormkoker = 1


class Dambreak:
    algorithm = 2
    timetobreachtomaximumdepth = 360  # s
    f1 = 1.3
    f2 = 0.04
    ucrit = 0.2


class FixedWeirs:
    pass


class MeasuredProfiles:
    pass


class Peil:
    boven_peil = np.nan
    onder_peil = np.nan
    vast_peil = np.nan


## Pumps
class Pumps:
    doelvariabele = "waterstand"
    maximalecapaciteit = 0
    peil_marge = 0.1
    streefwaarde = np.nan


## Weirs
class Weirs:
    afvoercoefficient_stuw = 1
    afvoercoefficient_opening = 0.85
    hoogstedoorstroombreedte = np.nan
    hoogstedoorstroomhoogte = np.nan
    laagstedoorstroombreedte = np.nan#0.1
    laagstedoorstroomhoogte = np.nan#0
    overlaatonderlaat = "Overlaat"
    soortregelbaarheid = 1
    soortstuw = 11
    vormopening = 3
