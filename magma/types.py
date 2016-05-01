# 2014.10.20 12:29:26 CEST
import os
import rdkit_engine as Chem
missingfragmentpenalty = 10 # non-optimized value

class ScanType(object):

    def __init__(self, scanid, mslevel):
        self.peaks = []
        self.scanid = scanid
        self.mslevel = mslevel


class PeakType(object):

    def __init__(self, mz, intensity, scanid, missing_fragment_score):
        self.mz = mz
        self.intensity = intensity
        self.scan = scanid
        self.childscan = None
        self.missing_fragment_score = missing_fragment_score


class HitType(object):

    def __init__(self, peak, fragment, score, bondbreaks, mass, ionmass, ion):
        self.mz = peak.mz
        self.intensity = peak.intensity
        self.intensity_weight = peak.missing_fragment_score / missingfragmentpenalty
        self.scan = peak.scan
        self.fragment = fragment
        self.score = score
        self.breaks = bondbreaks
        self.mass = mass
        self.deltaH = ionmass
        self.bonds = []
        self.allbonds = 0
        self.besthits = []
        self.atomstring = ''
        self.atomlist = []
        self.inchikey = ''
        self.formula = ''
        self.ion = ion


class MoleculeType(object):

    def __init__(self, molblock, name, prob, level, isquery = 1, mim = None, natoms = None, inchikey = None, molform = None, reference = None, logp = None):
        if inchikey == None or mim == None or molform == None or logp == None or natoms == None:
            mol = Chem.MolFromMolBlock(molblock)
            inchikey = Chem.MolToInchiKey(mol)[:14]
            (mim, molform,) = Chem.GetFormulaProps(mol)
            natoms = mol.GetNumHeavyAtoms()
            logp = Chem.LogP(mol)
        self.molblock = molblock
        self.level = level
        self.probability = prob
        self.inchikey = inchikey
        self.molformula = molform
        self.isquery = isquery
        self.name = name
        self.mim = mim
        self.natoms = natoms
        self.logp = logp
        self.reference = reference
