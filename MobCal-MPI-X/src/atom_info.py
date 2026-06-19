import numpy as np

chemical_symbols = [
    # 0
    'X',
    # 1
    'H', 'He',
    # 2
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    # 3
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    # 4
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    # 5
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    # 6
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn',
    # 7
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc',
    'Lv', 'Ts', 'Og']

atomic_masses_legacy = np.array([
    1.00000,  # X
    1.00794,  # H
    4.00260,  # He
    6.94100,  # Li
    9.01218,  # Be
    10.81100,  # B
    12.01100,  # C
    14.00670,  # N
    15.99940,  # O
    18.99840,  # F
    20.17970,  # Ne
    22.98977,  # Na
    24.30500,  # Mg
    26.98154,  # Al
    28.08550,  # Si
    30.97376,  # P
    32.06600,  # S
    35.45270,  # Cl
    39.94800,  # Ar
    39.09830,  # K
    40.07800,  # Ca
    44.95590,  # Sc
    47.88000,  # Ti
    50.94150,  # V
    51.99600,  # Cr
    54.93800,  # Mn
    55.84700,  # Fe
    58.93320,  # Co
    58.69340,  # Ni
    63.54600,  # Cu
    65.39000,  # Zn
    69.72300,  # Ga
    72.61000,  # Ge
    74.92160,  # As
    78.96000,  # Se
    79.90400,  # Br
    83.80000,  # Kr
    85.46780,  # Rb
    87.62000,  # Sr
    88.90590,  # Y
    91.22400,  # Zr
    92.90640,  # Nb
    95.94000,  # Mo
    np.nan,  # Tc
    101.07000,  # Ru
    102.90550,  # Rh
    106.42000,  # Pd
    107.86800,  # Ag
    112.41000,  # Cd
    114.82000,  # In
    118.71000,  # Sn
    121.75700,  # Sb
    127.60000,  # Te
    126.90450,  # I
    131.29000,  # Xe
    132.90540,  # Cs
    137.33000,  # Ba
    138.90550,  # La
    140.12000,  # Ce
    140.90770,  # Pr
    144.24000,  # Nd
    np.nan,  # Pm
    150.36000,  # Sm
    151.96500,  # Eu
    157.25000,  # Gd
    158.92530,  # Tb
    162.50000,  # Dy
    164.93030,  # Ho
    167.26000,  # Er
    168.93420,  # Tm
    173.04000,  # Yb
    174.96700,  # Lu
    178.49000,  # Hf
    180.94790,  # Ta
    183.85000,  # W
    186.20700,  # Re
    190.20000,  # Os
    192.22000,  # Ir
    195.08000,  # Pt
    196.96650,  # Au
    200.59000,  # Hg
    204.38300,  # Tl
    207.20000,  # Pb
    208.98040,  # Bi
    np.nan,  # Po
    np.nan,  # At
    np.nan,  # Rn
    np.nan,  # Fr
    226.02540,  # Ra
    np.nan,  # Ac
    232.03810,  # Th
    231.03590,  # Pa
    238.02900,  # U
    237.04820,  # Np
    np.nan,  # Pu
    np.nan,  # Am
    np.nan,  # Cm
    np.nan,  # Bk
    np.nan,  # Cf
    np.nan,  # Es
    np.nan,  # Fm
    np.nan,  # Md
    np.nan,  # No
    np.nan,  # Lr
    np.nan,  # Rf
    np.nan,  # Db
    np.nan,  # Sg
    np.nan,  # Bh
    np.nan,  # Hs
    np.nan,  # Mt
    np.nan,  # Ds
    np.nan,  # Rg
    np.nan,  # Cn
    np.nan,  # Nh
    np.nan,  # Fl
    np.nan,  # Mc
    np.nan,  # Lv
    np.nan,  # Ts
    np.nan,  # Og
])

missing = 0.2
covalent_radii = np.array([
    missing,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    1.68,  # Bk
    1.68,  # Cf
    1.65,  # Es
    1.67,  # Fm
    1.73,  # Md
    1.76,  # No
    1.61,  # Lr
    1.57,  # Rf
    1.49,  # Db
    1.43,  # Sg
    1.41,  # Bh
    1.34,  # Hs
    1.29,  # Mt
    1.28,  # Ds
    1.21,  # Rg
    1.22,  # Cn
    1.36,  # Nh
    1.43,  # Fl
    1.62,  # Mc
    1.75,  # Lv
    1.65,  # Ts
    1.57,  # Og
])

vdw_radii = np.array([
    np.nan,  # X
    1.20,  # H
    1.40,  # He [1]
    1.82,  # Li [1]
    1.53,  # Be [5]
    1.92,  # B [5]
    1.70,  # C [1]
    1.55,  # N [1]
    1.52,  # O [1]
    1.47,  # F [1]
    1.54,  # Ne [1]
    2.27,  # Na [1]
    1.73,  # Mg [1]
    1.84,  # Al [5]
    2.10,  # Si [1]
    1.80,  # P [1]
    1.80,  # S [1]
    1.75,  # Cl [1]
    1.88,  # Ar [1]
    2.75,  # K [1]
    2.31,  # Ca [5]
    np.nan,  # Sc
    np.nan,  # Ti
    np.nan,  # V
    np.nan,  # Cr
    np.nan,  # Mn
    np.nan,  # Fe
    np.nan,  # Co
    1.63,  # Ni [1]
    1.40,  # Cu [1]
    1.39,  # Zn [1]
    1.87,  # Ga [1]
    2.11,  # Ge [5]
    1.85,  # As [1]
    1.90,  # Se [1]
    1.85,  # Br [1]
    2.02,  # Kr [1]
    3.03,  # Rb [5]
    2.49,  # Sr [5]
    np.nan,  # Y
    np.nan,  # Zr
    np.nan,  # Nb
    np.nan,  # Mo
    np.nan,  # Tc
    np.nan,  # Ru
    np.nan,  # Rh
    1.63,  # Pd [1]
    1.72,  # Ag [1]
    1.58,  # Cd [1]
    1.93,  # In [1]
    2.17,  # Sn [1]
    2.06,  # Sb [5]
    2.06,  # Te [1]
    1.98,  # I [1]
    2.16,  # Xe [1]
    3.43,  # Cs [5]
    2.49,  # Ba [5]
    np.nan,  # La
    np.nan,  # Ce
    np.nan,  # Pr
    np.nan,  # Nd
    np.nan,  # Pm
    np.nan,  # Sm
    np.nan,  # Eu
    np.nan,  # Gd
    np.nan,  # Tb
    np.nan,  # Dy
    np.nan,  # Ho
    np.nan,  # Er
    np.nan,  # Tm
    np.nan,  # Yb
    np.nan,  # Lu
    np.nan,  # Hf
    np.nan,  # Ta
    np.nan,  # W
    np.nan,  # Re
    np.nan,  # Os
    np.nan,  # Ir
    1.75,  # Pt [1]
    1.66,  # Au [1]
    1.55,  # Hg [1]
    1.96,  # Tl [1]
    2.02,  # Pb [1]
    2.07,  # Bi [5]
    1.97,  # Po [5]
    2.02,  # At [5]
    2.20,  # Rn [5]
    3.48,  # Fr [5]
    2.83,  # Ra [5]
    np.nan,  # Ac
    np.nan,  # Th
    np.nan,  # Pa
    1.86,  # U [1]
    np.nan,  # Np
    np.nan,  # Pu
    np.nan,  # Am
    np.nan,  # Cm
    np.nan,  # Bk
    np.nan,  # Cf
    np.nan,  # Es
    np.nan,  # Fm
    np.nan,  # Md
    np.nan,  # No
    np.nan])  # Lr