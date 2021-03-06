# 2014.10.20 12:29:04 CEST
typew = {'AROMATIC': 3.0,
 'DOUBLE': 2.0,
 'TRIPLE': 3.0,
 'SINGLE': 1.0}
heterow = {False: 2,
 True: 1}
missingfragmentpenalty = 10.0

mims = {'H': 1.0078250321,
 'He': 3.016029,
 'Li': 6.015122,
 'Be': 9.012182,
 'B': 10.012937,
 'C': 12.000000,
 'N': 14.0030740052,
 'O': 15.9949146221,
 'F': 18.9984032,
 'Ne': 19.992440,
 'Na': 22.9897692809,
 'Mg': 23.985042,
 'Al': 26.981538,
 'Si': 27.976927,
 'P': 30.97376151,
 'S': 31.97207069,
 'Cl': 34.96885271,
 'Ar': 35.967546,
 'K': 38.96370668,
 'Ca': 39.962591,
 'Sc': 44.955910,
 'Ti': 45.952629,
 'V': 49.947163,
 'Cr': 49.946050,
 'Mn': 54.938050,
 'Fe': 53.939615,
 'Co': 58.933200,
 'Ni': 57.935348,
 'Cu': 62.929601,
 'Zn': 63.929147,
 'Ga': 68.925581,
 'Ge': 69.924250,
 'As': 74.921596,
 'Se': 73.922477,
 'Br': 78.9183376,
 'Kr': 77.920386,
 'Rb': 84.911789,
 'Sr': 83.913425,
 'Y': 88.905848,
 'Zr': 89.904704,
 'Nb': 92.906378,
 'Mo': 91.906810,
 'Tc': 97.907216,
 'Ru': 95.907598,
 'Rh': 102.905504,
 'Pd': 101.905608,
 'Ag': 106.905093,
 'Cd': 105.906458,
 'In': 112.904061,
 'Sn': 111.904821,
 'Sb': 120.903818,
 'Te': 119.904020,
 'I': 126.904468,
 'Xe': 123.905896,
 'Cs': 132.905447,
 'Ba': 129.906310,
 'La': 137.907107,
 'Ce': 135.907144,
 'Pr': 140.907648,
 'Nd': 141.907719,
 'Pm': 144.912744,
 'Sm': 143.911995,
 'Eu': 150.919846,
 'Gd': 151.919788,
 'Tb': 158.925343,
 'Dy': 155.924278,
 'Ho': 164.930319,
 'Er': 161.928775,
 'Tm': 168.934211,
 'Yb': 167.933894,
 'Lu': 174.940768,
 'Hf': 173.940040,
 'Ta': 179.947466,
 'W': 179.946706,
 'Re': 184.952956,
 'Os': 183.952491,
 'Ir': 190.960591,
 'Pt': 189.959930,
 'Au': 196.966552,
 'Hg': 195.965815,
 'Tl': 202.972329,
 'Pb': 203.973029,
 'Bi': 208.980383}

Hmass = mims['H']
elmass = 0.0005486
ionmasses = {1: {'+H': mims['H'],
     '+NH4': mims['N'] + 4 * mims['H'],
     '+Na': mims['Na'],
     '-OH': -(mims['O'] + mims['H']),
     '+K': mims['K']},
 -1: {'-H': -mims['H'],
      '+Cl': mims['Cl']}}
