#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Authors: Dries Verdegem <dries.verdegem@vib-kuleuven.be>
#
# License: BSD 3 clause

'''
MAGMa+ allows for MAGMa metabolite identifications with dynamic parameter selection.

Copyright (C) 2015 Dries Verdegem and VIB - KU Leuven (BSD-3)

December 2015: Original open-source release
'''

__version__ = '1.0.0'
__author__ = "Dries verdegem"
__copyright__ = "Copyright 2015, VIB-KU Leuven (http://www.vib.be/en, http://www.kuleuven.be/english), Dries Verdegem"
__credits__ = ["Dries Verdegem", "Diether Lambrechts", "Peter Carmeliet", "Bart Ghesquière"]
__license__ = "BSD-3"
__maintainer__ = "Dries Verdegem"
__email__ = "dries.verdegem@vib-kuleuven.be"
__status__ = "Production"


import os, sys
import re
import pickle
import zlib
import sqlite3
from cStringIO import StringIO


from rdkit import Chem
import Allkeys


def version():
    return __version__


def get_features(ionization_mode):
    rfecv_pos_features = [25, 26, 42, 50, 53, 54, 57, 62, 65, 66, 69, 72, 74, 75, 76, 77, 79, 80, 82,
                          83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100,
                          101, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                          117, 118, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132,
                          133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147,
                          148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
                          163, 165, 167, 168, 169, 171, 178, 179, 180, 181, 184, 189, 190, 194, 207,
                          215, 250, 251, 252, 254, 291, 301, 303, 343, 345, 347, 348, 350, 440, 441,
                          444, 445, 446, 453, 455, 456, 462, 463, 465, 466, 467, 468, 470, 473, 474,
                          475, 477, 499, 500, 501, 502, 503, 506, 509, 510, 511, 517, 519, 520, 521,
                          522, 523, 526, 558, 592, 603, 621, 624, 632, 634, 635, 641, 647, 650, 660,
                          681, 682, 687, 705, 707, 708, 709, 724, 728, 732, 739, 746, 747, 748, 749,
                          752, 756, 757, 758, 760, 761, 764, 765, 766, 769, 772, 773, 776, 777, 779,
                          803, 824, 829, 845, 865, 886, 908, 928, 949, 950, 956, 962, 965, 968, 973,
                          1020, 1040, 1074, 1081, 1137, 1238, 1245, 1246, 1247, 1248, 1250, 1286,
                          1289, 1292, 1294, 1295, 1296, 1297, 1298, 1299, 1302, 1306, 1310, 1316,
                          1317, 1318, 1330, 1346, 1349, 1357, 1358, 1414, 1442, 1452, 1496, 1505,
                          1595, 1602, 1617, 1625, 1626, 1639, 1655, 1736, 1838, 2094, 2096, 2102,
                          2105, 2141, 2166, 2211, 2354, 2512, 2537, 2590, 2593, 2594, 2880, 3208,
                          3212, 3221, 3495, 3512, 3615, 3621, 3630, 3897, 3923, 3927, 3934, 3958,
                          3973, 4098, 4108, 4172, 4216, 4231, 4276, 4308, 4317, 4323, 4335, 4337,
                          4338, 4340, 4341, 4343, 4347, 4350, 4356, 4384, 4388, 4389, 4403, 4422,
                          4539, 4541, 4554, 4588, 4591, 4594, 4602, 4603, 4609, 4619, 4629, 4630,
                          4645, 4649, 4652, 4654, 4656, 4657, 4658, 4660, 4664, 4666, 4667, 4669,
                          4670, 4673, 4675, 4678, 4685, 4688, 4692, 4698, 4721, 4730, 4734, 4736,
                          4737, 4738, 4830, 4874, 4885, 4890, 4904, 4905, 4907, 4919, 4948, 4980,
                          5028, 5185, 5222, 5243, 5358, 5366, 5643, 5665, 5686, 5705, 5718, 5758,
                          5761, 5777, 5791]
    
    rfecv_neg_features = [19, 25, 29, 34, 37, 38, 42, 43, 48, 49, 50, 52, 53, 54, 57, 62, 63, 65, 66,
                          69, 70, 71, 72, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88,
                          89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104, 105, 106,
                          107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
                          122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
                          137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                          152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 165, 167, 168,
                          169, 170, 171, 178, 179, 180, 181, 182, 184, 189, 190, 191, 194, 207, 209,
                          215, 218, 222, 229, 231, 232, 240, 250, 251, 252, 254, 264, 265, 266, 268,
                          281, 292, 301, 302, 303, 304, 309, 313, 335, 336, 343, 345, 346, 347, 348,
                          350, 354, 403, 412, 440, 441, 442, 444, 445, 446, 449, 451, 452, 453, 455,
                          456, 457, 462, 463, 464, 465, 466, 467, 468, 470, 473, 474, 475, 477, 499,
                          500, 501, 502, 503, 506, 507, 509, 510, 511, 517, 519, 520, 521, 522, 523,
                          524, 526, 530, 544, 558, 566, 592, 603, 604, 605, 621, 624, 631, 632, 634,
                          635, 641, 647, 649, 650, 660, 662, 678, 681, 682, 687, 705, 706, 707, 708,
                          709, 710, 713, 724, 727, 728, 729, 730, 732, 734, 739, 746, 747, 748, 749,
                          750, 751, 752, 754, 755, 756, 757, 758, 760, 761, 762, 764, 765, 766, 767,
                          768, 769, 770, 772, 773, 776, 777, 778, 779, 780, 781, 782, 787, 793, 802,
                          803, 808, 823, 824, 829, 831, 844, 845, 865, 866, 886, 908, 913, 928, 929,
                          949, 950, 956, 958, 962, 964, 965, 966, 968, 973, 1020, 1022, 1039, 1040,
                          1041, 1046, 1047, 1052, 1059, 1074, 1081, 1084, 1137, 1162, 1235, 1238,
                          1245, 1246, 1247, 1248, 1250, 1252, 1255, 1259, 1263, 1273, 1286, 1289,
                          1292, 1294, 1295, 1296, 1297, 1298, 1300, 1302, 1305, 1306, 1310, 1311,
                          1314, 1315, 1316, 1318, 1323, 1330, 1331, 1346, 1349, 1357, 1392, 1414,
                          1415, 1424, 1428, 1442, 1452, 1495, 1496, 1514, 1547, 1595, 1602, 1617,
                          1620, 1625, 1626, 1627, 1628, 1631, 1637, 1638, 1639, 1642, 1655, 1745,
                          1758, 1787, 1838, 2094, 2095, 2096, 2102, 2104, 2105, 2108, 2109, 2141,
                          2166, 2176, 2177, 2188, 2193, 2211, 2212, 2213, 2221, 2243, 2353, 2354,
                          2355, 2374, 2375, 2378, 2400, 2468, 2472, 2512, 2514, 2519, 2540, 2590,
                          2593, 2594, 2601, 2672, 2715, 2717, 2868, 2879, 2880, 2881, 2985, 3061,
                          3080, 3083, 3085, 3134, 3208, 3210, 3213, 3328, 3495, 3496, 3498, 3512,
                          3513, 3543, 3546, 3614, 3615, 3621, 3643, 3659, 3660, 3669, 3719, 3803,
                          3882, 3897, 3923, 3927, 3932, 3934, 3958, 3961, 3973, 3982, 4002, 4026,
                          4087, 4091, 4097, 4098, 4099, 4104, 4108, 4128, 4130, 4158, 4172, 4216,
                          4229, 4231, 4233, 4236, 4243, 4276, 4284, 4287, 4299, 4308, 4309, 4317,
                          4319, 4321, 4323, 4329, 4330, 4335, 4336, 4337, 4338, 4340, 4341, 4342,
                          4343, 4344, 4346, 4347, 4350, 4352, 4356, 4358, 4367, 4373, 4374, 4375,
                          4376, 4378, 4382, 4384, 4388, 4389, 4403, 4422, 4438, 4468, 4477, 4482,
                          4485, 4488, 4502, 4507, 4508, 4509, 4512, 4517, 4522, 4527, 4530, 4539,
                          4541, 4544, 4546, 4554, 4556, 4557, 4583, 4588, 4591, 4593, 4594, 4602,
                          4603, 4605, 4607, 4608, 4609, 4610, 4615, 4616, 4619, 4623, 4629, 4630,
                          4633, 4636, 4640, 4645, 4647, 4649, 4650, 4651, 4652, 4654, 4656, 4657,
                          4658, 4660, 4664, 4666, 4667, 4669, 4670, 4673, 4675, 4676, 4677, 4678,
                          4681, 4683, 4685, 4688, 4690, 4692, 4694, 4698, 4705, 4721, 4723, 4729,
                          4730, 4736, 4737, 4738, 4754, 4830, 4874, 4876, 4885, 4886, 4890, 4891,
                          4894, 4896, 4904, 4905, 4907, 4915, 4919, 4945, 4946, 4948, 4953, 4963,
                          4980, 5028, 5111, 5185, 5194, 5195, 5216, 5222, 5234, 5235, 5241, 5243,
                          5279, 5284, 5332, 5356, 5358, 5364, 5366, 5401, 5438, 5443, 5444, 5467,
                          5468, 5469, 5470, 5471, 5472, 5479, 5532, 5558, 5560, 5581, 5593, 5598,
                          5611, 5612, 5613, 5615, 5616, 5626 ,5637, 5643, 5656, 5665, 5684, 5695,
                          5700, 5701, 5702, 5705, 5711, 5712, 5718, 5751, 5758, 5759, 5761, 5765,
                          5766, 5768, 5769, 5774, 5777, 5778, 5784, 5791]
    
    if ionization_mode == '1':
        return rfecv_pos_features
    elif ionization_mode == '-1':
        return rfecv_neg_features


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


def get_params(ionization_mode, opt_index):
    if ionization_mode == '1':
        if opt_index == 'A':
            DBBC = 1.7
            ABBC = 1.2
            MSP = 7.7
        elif opt_index == 'x':
            DBBC = 0.7
            ABBC = 0.9
            MSP = 6.5
    elif ionization_mode == '-1':
        if opt_index == 'A':
            DBBC = 7.0
            ABBC = 1.9
            MSP = 10.5
        elif opt_index == 'x':
            DBBC = 1.9
            ABBC = 0.9
            MSP = 7.2
    
    return DBBC, ABBC, MSP


def set_pars(magma_script, ionization_mode, opt_index):
    DBBC, ABBC, MSP = get_params(ionization_mode, opt_index)
    
    for par in magma_script.magma.pars.typew.keys():
        par_str = str(par)
        if par_str == 'DOUBLE':
            magma_script.magma.pars.typew[par] = DBBC
        elif par_str == 'AROMATIC':
            magma_script.magma.pars.typew[par] = ABBC
    
    magma_script.magma.pars.missingfragmentpenalty = MSP


def get_top_ranked(matches, topk = 3):
    top_matches = []
    for _, _, mols in matches:
        top_matches.extend(mols)
        if len(top_matches) >= topk: break
    return top_matches


def get_classifier_path(classifier_dir, ionization_mode):
    if ionization_mode == '1':
        mode = 'pos'
        label = 'D'
    elif ionization_mode == '-1':
        mode = 'neg'
        label = 'C'
    classifier_path = os.path.join(classifier_dir, 'metabolite_%s%s_rfc.pkl' % (mode, label))
    return classifier_path
    

def get_molclass(structuresA, structuresx, ionization_mode, classifier_dir):
    # if both parameter sets generated the same output, the "A" parameter is choosen
    outputA = [set(s[1]) for s in structuresA]
    outputx = [set(s[1]) for s in structuresx]
    if outputA == outputx:
        print 'Equivalent outcome for both parameter settings.'
        return 'A'
    
    # otherwise, determine the molclass of the top three ranked metabolites
    # for both parameter sets using the molecular classifier
    
    # getting the feature list to use
    feature_list = get_features(ionization_mode)
    # getting the classifier to use
    classifier_path = get_classifier_path(classifier_dir, ionization_mode)
    with open(classifier_path, 'r') as rfc_file:
        clf = pickle.load(rfc_file)
    
    molsA = get_top_ranked(structuresA)
    molsx = get_top_ranked(structuresx)
    
    featuresA = [[int(k) for k in list(Allkeys._pyGenFilteredMACCSKeys(mol, feature_list).ToBitString())] for mol in molsA]
    featuresx = [[int(k) for k in list(Allkeys._pyGenFilteredMACCSKeys(mol, feature_list).ToBitString())] for mol in molsx]
    
    membershipA = [clf.predict(features)[0] for features in featuresA]
    membershipx = [clf.predict(features)[0] for features in featuresx]
    
    print 'membershipA', membershipA
    print 'membershipx', membershipx
    
    memcA = membershipA.count(0) + membershipx.count(0)
    memcx = membershipA.count(1) + membershipx.count(1)
    
    if memcA >= memcx:
        return 'A'
    elif memcx > memcA:
        return 'x'


def parse_export_structures(structures_txt):
    # first get database file name
    with open('structure_db.txt', 'r') as structure_db_f:
        structure_db_fn = structure_db_f.readline()
    # connecting to database
    structure_db = sqlite3.connect(structure_db_fn)
    cursor = structure_db.cursor()
    # defining the metabolite id regular expression
    metabolite_id_regex = re.compile(r'\((.+)\)\s+[^\s]+$')
    # filtering the input
    try:
        structures_txt = structures_txt.split('Candidate_Score Name Smiles\n')[1]
    except:
        return []
    structures_temp = []
    for structure_info in structures_txt.split('\n'):
        metabolite_id = metabolite_id_regex.findall(structure_info)[0]
        # get item from database
        cursor.execute('''SELECT name, molblock FROM molecules WHERE id=?''', (metabolite_id,))
        try:
            name, molblock = cursor.fetchone()
        except:
            print metabolite_id, structure_db_fn
            exit()
        molblock = zlib.decompress(molblock)
        mol = Chem.MolFromMolBlock(molblock)
        score = float(structure_info.split(' ')[0])
        structures_temp.append((score, name, mol))
    # closing database
    structure_db.close()

    structures_dict = {}
    for score, name, mol in structures_temp:
        if score not in structures_dict:
            structures_dict[score] = [score, [name], [mol]]
        else:
            structures_dict[score][1].append(name)
            structures_dict[score][2].append(mol)
    
    structures = structures_dict.values()
    structures.sort(key = lambda x: x[0])
    
    return structures


def main():
    # finding out where the the classifiers are located
    try:
        magma_plus_classifier_path = os.environ['MAGMAPLUS_CLASSIFIER_PATH']
    except:
        magma_plus_classifier_path = os.getcwd()
    
    error = ''
    pos_clf_fn = get_classifier_path(magma_plus_classifier_path, '1')
    neg_clf_fn = get_classifier_path(magma_plus_classifier_path, '-1')
    fingerprint_fn = os.path.join(magma_plus_classifier_path, 'AllFingerprints.txt')
    if not os.path.isfile(pos_clf_fn):
        error += 'Positive ionization mode classifier (%s) cannot be found.\n' % pos_clf_fn
    if not os.path.isfile(neg_clf_fn):
        error += 'Negative ionization mode classifier (%s) cannot be found.\n' % neg_clf_fn
    if not os.path.isfile(fingerprint_fn):
        error += 'Fingerprint file (%s) cannot be found.\n' % fingerprint_fn
    if error:
        error += 'Either set the environment variable MAGMAPLUS_CLASSIFIER_PATH to the directory containing these files.\n'
        error += 'Or copy these files to the current directory.\n'
        error += 'See the README file.'
        raise LookupError(error)
    else:
        # initializing Allkeys
        Allkeys.InitPatts(magma_plus_classifier_path)

    # getting the tool
    import MAGMa_script
    magmacommand = MAGMa_script.MagmaCommand()
    
    # getting the arguments
    magma_arguments_l = sys.argv[1:]
    
    if len(magma_arguments_l) == 0:
        magma_arguments_l = ['-h']
    
    # check if help argument was given
    if '-h' in magma_arguments_l or '--help' in magma_arguments_l:
        magmacommand.run(magma_arguments_l)
        exit()
    
    # parsing some arguments
    action = magma_arguments_l[0]
    db_fn = magma_arguments_l[-1]
    
    # performing MAGMa identification twice
    out = {}

    for opt_index in ['A', 'x']:
                
        if action == 'read_ms_data':
            # determine the ionization mode
            # (based on the -i or --ionisation_mode option)
            ionization_arg_keys = ['-i', '--ionisation_mode']
            ionization_mode = '1'
            for iak in ionization_arg_keys:
                if iak in magma_arguments_l:
                    iaki = magma_arguments_l.index(iak)
                    try:
                        ionization_mode = magma_arguments_l[iaki + 1]
                    except:
                        ionization_mode = '1' # positive ionization mode is default
                    else:
                        break
            
            # adapting the database file name
            if ionization_mode == '1': # positive
                db_fn_sub = db_fn + '.pos.%s' % opt_index
            elif ionization_mode == '-1': # negative
                db_fn_sub = db_fn + '.neg.%s' % opt_index
        
        else:
            # determine the ionization mode
            # (based on the available db_fn)
            ionization_mode = '1'
            # get the creation date of the two possible database variants
            db_fn_pos = db_fn + '.pos.%s' % opt_index
            db_fn_neg = db_fn + '.neg.%s' % opt_index
            if os.path.isfile(db_fn_pos):
                postime = os.path.getmtime(db_fn_pos)
            else:
                postime = 0.
            if os.path.isfile(db_fn_neg):
                negtime = os.path.getmtime(db_fn_neg)
            else:
                negtime = 0.
            if postime > negtime:
                ionization_mode = '1'
                db_fn_sub = db_fn_pos
            else:
                ionization_mode = '-1'
                db_fn_sub = db_fn_neg
        
        if action == 'annotate':
            structure_db_arg_keys = ['-s', '--structure_database']
            structure_db = None
            for sdak in structure_db_arg_keys:
                if sdak in magma_arguments_l:
                    sdaki = magma_arguments_l.index(sdak)
                    try:
                        structure_db = magma_arguments_l[sdaki + 1]
                    except:
                        structure_db = None # unknown structure db
                    else:
                        break
            # open the magma_job.ini file to obtain the structure db file name
            structure_db_fn = ''
            with open('magma_job.ini', 'r') as magma_job_f:
                magma_job_content = magma_job_f.readlines()
            for line in magma_job_content:
                if structure_db == None:
                    if line.startswith('structure_database'):
                        # the first db is taken
                        structure_db_fn = line.split('=')[1].strip()
                        break
                else:
                    if line.startswith('structure_database.%s' % structure_db):
                        structure_db_fn = line.split('=')[1].strip()
                        break
            # store the structure db file name in a file
            with open('structure_db.txt', 'w') as structure_db_f:
                structure_db_f.write(structure_db_fn)


        # adapting the magma parameters
        set_pars(MAGMa_script, ionization_mode, opt_index)
    
        # running magma
        magma_arguments_l[-1] = db_fn_sub
        print 'Parameter settings %s' % opt_index
        with Capturing() as out[opt_index]:
            magmacommand.run(magma_arguments_l)

        if action != 'export_result':
            for line in out[opt_index]:
                print line

    # deciding which of the structures to output
    if action == 'export_result':
        structures = {}
        for opt_index in ['A', 'x']:
            structures[opt_index] = parse_export_structures('\n'.join(out[opt_index]))
            
        # check which parameters (A or x) to choose
        molclass = get_molclass(structures['A'], structures['x'], ionization_mode, magma_plus_classifier_path)
        print 'Final results'
        for line in out[molclass]:
            print line
                
        print 'MAGMa+ identification done using: DBBC: %s - ABBC: %s - MSP: %s' % (get_params(ionization_mode, molclass))
            

if __name__ == '__main__':
    main()