# $Id$
#
# Copyright (C) 2001-2011 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" SMARTS definitions for the publically available MACCS keys
and a MACCS fingerprinter

I compared the MACCS fingerprints generated here with those from two
other packages (not MDL, unfortunately). Of course there are
disagreements between the various fingerprints still, but I think
these definitions work pretty well. Some notes:

1) most of the differences have to do with aromaticity
2) there's a discrepancy sometimes because the current RDKit
definitions do not require multiple matches to be distinct. e.g. the
SMILES C(=O)CC(=O) can match the (hypothetical) key O=CC twice in my
definition. It's not clear to me what the correct behavior is.
3) Some keys are not fully defined in the MDL documentation
4) Two keys, 125 and 166, have to be done outside of SMARTS.
5) Key 1 (ISOTOPE) isn't defined

Rev history:
2006 (gl): Original open-source release
May 2011 (gl): Update some definitions based on feedback from Andrew Dalke
Oct 2015: Update for compatibility with all SMARTS-formatted fingerprints
          from external file (with filter) by Dries Verdegem

"""
import os

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs



smartsPatts = {}

def InitPatts(magma_classifier_path):
    # getting the keys
    global smartsPatts
    all_fingerprints_f = open(os.path.join(magma_classifier_path, 'AllFingerprints.txt'),'r')
    key_num = 1
    for line in all_fingerprints_f:
        line = line.strip()
        if line.startswith('#'): continue
        if not line: continue
        fingerprint, count = line.split('\t')
        count = int(count)
        smartsPatts[key_num] = (fingerprint, count)
        key_num += 1
    all_fingerprints_f.close()

maccsKeys = None

def _InitKeys(keyList,keyDict):
    """ *Internal Use Only*

       generates SMARTS patterns for the keys, run once

    """
    assert len(keyList) == len(keyDict.keys()),'length mismatch'
    for key in keyDict.keys():
        patt,count = keyDict[key]
        if patt != '?':
            try:
                sma = Chem.MolFromSmarts(patt)
            except:
                sma = None
            if not sma:
                print 'SMARTS parser error for key #%d: %s'%(key,patt)
            else:
                keyList[key-1] = sma,count
      
def _pyGenMACCSKeys(mol,**kwargs):
    """ generates the MACCS fingerprint for a molecules

     **Arguments**

     - mol: the molecule to be fingerprinted

     - any extra keyword arguments are ignored
     
     **Returns**

      a _DataStructs.SparseBitVect_ containing the fingerprint.

    >>> m = Chem.MolFromSmiles('CNO')
    >>> bv = GenMACCSKeys(m)
    >>> tuple(bv.GetOnBits())
    (24, 68, 69, 71, 93, 94, 102, 124, 131, 139, 151, 158, 160, 161, 164)
    >>> bv = GenMACCSKeys(Chem.MolFromSmiles('CCC'))
    >>> tuple(bv.GetOnBits())
    (74, 114, 149, 155, 160)

    """
    global maccsKeys
    if maccsKeys is None:
        maccsKeys = [(None,0)]*len(smartsPatts.keys())
        _InitKeys(maccsKeys,smartsPatts)
    ctor=kwargs.get('ctor',DataStructs.SparseBitVect)

    res = ctor(len(maccsKeys)+1)
    for i,(patt,count) in enumerate(maccsKeys):
        if patt is not None:
            if count==0:
                res[i+1] = mol.HasSubstructMatch(patt)
            else:
                matches = mol.GetSubstructMatches(patt)
                if len(matches) > count:
                    res[i+1] = 1
        elif (i+1)==125:
            # special case: num aromatic rings > 1
            ri = mol.GetRingInfo()
            nArom=0
            res[125]=0
            for ring in ri.BondRings():
                isArom=True
                for bondIdx in ring:
                    if not mol.GetBondWithIdx(bondIdx).GetIsAromatic():
                        isArom=False
                        break
                if isArom:
                    nArom+=1
                    if nArom>1:
                        res[125]=1
                        break
        elif (i+1)==166:
            res[166]=0
            # special case: num frags > 1
            if len(Chem.GetMolFrags(mol))>1:
                res[166]=1
          
    return res

def _pyGenFilteredMACCSKeys(mol, feature_filter=[], **kwargs):
    """ generates the MACCS fingerprint for a molecules

     **Arguments**

     - mol: the molecule to be fingerprinted

     - any extra keyword arguments are ignored
     
     **Returns**

      a _DataStructs.SparseBitVect_ containing the fingerprint.

    >>> m = Chem.MolFromSmiles('CNO')
    >>> bv = GenMACCSKeys(m)
    >>> tuple(bv.GetOnBits())
    (24, 68, 69, 71, 93, 94, 102, 124, 131, 139, 151, 158, 160, 161, 164)
    >>> bv = GenMACCSKeys(Chem.MolFromSmiles('CCC'))
    >>> tuple(bv.GetOnBits())
    (74, 114, 149, 155, 160)

    """
    global maccsKeys
    if maccsKeys is None:
        maccsKeys = [(None,0)]*len(smartsPatts.keys())
        _InitKeys(maccsKeys,smartsPatts)
    ctor=kwargs.get('ctor',DataStructs.SparseBitVect)

    res = ctor(len(feature_filter)+1)
    
    res_count = -1
    for i,(patt,count) in enumerate(maccsKeys):
        if i+1 in feature_filter:
            #print 'patt', patt, smartsPatts[i+1]
            res_count += 1
        else:
            continue
        if patt is not None:
            if count==0:
                res[res_count+1] = mol.HasSubstructMatch(patt)
            else:
                matches = mol.GetSubstructMatches(patt)
                if len(matches) > count:
                    res[res_count+1] = 1
        elif (i+1)==125:
            # special case: num aromatic rings > 1
            ri = mol.GetRingInfo()
            nArom=0
            res[res_count+1]=0
            for ring in ri.BondRings():
                isArom=True
                for bondIdx in ring:
                    if not mol.GetBondWithIdx(bondIdx).GetIsAromatic():
                        isArom=False
                        break
                if isArom:
                    nArom+=1
                    if nArom>1:
                        res[res_count+1]=1
                        break
        elif (i+1)==166:
            res[res_count+1]=0
            # special case: num frags > 1
            if len(Chem.GetMolFrags(mol))>1:
                res[res_count+1]=1
          
    return res

GenMACCSKeys = rdMolDescriptors.GetMACCSKeysFingerprint
FingerprintMol = rdMolDescriptors.GetMACCSKeysFingerprint
  
#------------------------------------
#
#  doctest boilerplate
#
def _test():
    import doctest,sys
    return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
    import sys
    failed,tried = _test()
    sys.exit(failed)