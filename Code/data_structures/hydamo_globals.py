BRANCH_FRICTION_FUNCTION = {
    "Constant": "constant",
    "Discharge": "absDischarge",
    "Waterlevel": "waterLevel",
}
HYDAMO_SHAPE_NUMS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 98, 99]
HYDAMO_WEIR_TYPES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24, 25, 26, 98, 99]
MANAGEMENT_DEVICE_TYPES = {
    "niet regelbaar (vast)": 1,
    "regelbaar, niet automatisch": 2,
    "regelbaar, automatisch": 3,
    "handmatig": 4,
    "overig": 98,
    "onbekend": 99,
}
ROUGHNESS_MAPPING = {
    "Chezy": "Chezy",
    "Manning": "Manning",
    "StricklerKn": "StricklerNikuradse",
    "StricklerKs": "Strickler",
    "White Colebrook": "WhiteColebrook",
    "Bos en Bijkerk": "deBosBijkerk",
    "Onbekend": "Strickler",
    "Overig": "Strickler",
}
ROUGHNESS_MAPPING_LIST = list(ROUGHNESS_MAPPING)
WEIR_MAPPING = {
    "schotbalkstuw": 1,
    "stuw met schuif": 2,
    "stuw met klep": 3,
    "segmentstuw": 4,
    "cascadestuw": 5,
    "hevelstuw": 6,
    "meetstuw": 7,
    "meetschot": 8,
    "stuw met contra-gewicht": 9,
    "inlaat- en/of aflaatstuw": 10,
    "overlaat": 11,
    "drijverstuw": 12,
    "trommelstuw": 13,
    "gronddamstuw": 20,
    "stuwbak": 21,
    "tuimel- of kantelstuw": 22,
    "balgstuw": 23,
    "brievenbusstuw": 24,
    "knijpstuw": 25,
    "conserveringstuw": 26,
    "overig": 98,
    "onbekend": 99,
}
