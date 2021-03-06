# -*- coding: utf-8 -*-

# 2014.10.20 12:28:15 CEST
import sys
import base64
import subprocess
import StringIO
import time
import re
import os
import sqlite3
import struct
import zlib
import gzip
import copy
import pkg_resources
import numpy
import logging
from lxml import etree
from sqlalchemy import create_engine, and_, desc, distinct
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func
from models import Base, Metabolite, Reaction, fill_molecules_reactions, Scan, Peak, Fragment, Run
import requests
import functools
#import macauthlib
from requests.auth import AuthBase
import cPickle as pickle
import types
import pars
from operator import itemgetter
import ConfigParser
config = ConfigParser.ConfigParser()
config.add_section('magma job')
config.set('magma job', 'chemical_engine', 'rdkit')
config.read(['magma_job.ini', os.path.expanduser('~/magma_job.ini')])
if config.get('magma job', 'chemical_engine') == 'rdkit':
    import rdkit_engine as Chem
elif config.get('magma job', 'chemical_engine') == 'cdk':
    import cdk_engine
    Chem = cdk_engine.engine()

class MagmaSession(object):

    def __init__(self, db_name, description = '', loglevel = 'warning'):
        engine = create_engine('sqlite:///' + db_name, echo=False)
        session = sessionmaker()
        session.configure(bind=engine)
        self.db_session = session()
        Base.metadata.create_all(engine)
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.description == None:
            rundata.description = unicode(description)
        self.db_session.add(rundata)
        self.db_session.commit()
        logging.basicConfig(format='%(levelname)s: %(message)s', level=getattr(logging, loglevel.upper()))



    def get_structure_engine(self, pubchem_names = False, call_back_url = None):
        return StructureEngine(self.db_session, pubchem_names, call_back_url)



    def get_ms_data_engine(self, ionisation_mode = 1, abs_peak_cutoff = 1000, mz_precision = 5.0, mz_precision_abs = 0.001, precursor_mz_precision = 0.005, max_ms_level = 10, call_back_url = None):
        return MsDataEngine(self.db_session, ionisation_mode, abs_peak_cutoff, mz_precision, mz_precision_abs, precursor_mz_precision, max_ms_level, call_back_url)



    def get_annotate_engine(self, skip_fragmentation = False, max_broken_bonds = 3, max_water_losses = 1, ms_intensity_cutoff = 1000000.0, msms_intensity_cutoff = 5, use_all_peaks = False, adducts = None, force_adduct = False, max_charge = 1, call_back_url = None):
        return AnnotateEngine(self.db_session, skip_fragmentation, max_broken_bonds, max_water_losses, ms_intensity_cutoff, msms_intensity_cutoff, use_all_peaks, adducts, force_adduct, max_charge, call_back_url)



    def get_select_engine(self):
        return SelectEngine(self.db_session)



    def get_data_analysis_engine(self):
        return DataAnalysisEngine(self.db_session)



    def get_call_back_engine(self, id, key):
        return CallBackEngine(id, key)



    def fill_molecules_reactions(self):
        fill_molecules_reactions(self.db_session)



    def commit(self):
        self.db_session.commit()



    def close(self):
        self.db_session.close()




class CallBackEngine(object):

    def __init__(self, url):
        self.access_token = config.get('magma job', 'macs.id')
        self.mac_key = config.get('magma job', 'macs.key')
        self.url = url
        self.update_interval = 2
        self.update_time = time.time()



    def update_callback_url(self, status, elapsed_time = None, time_limit = None, force = False):

        class HTTPMacAuth(AuthBase):
            """Attaches HTTP Basic Authentication to the given Request object."""


            def __init__(self, id, key):
                self.id = id
                self.key = key



            def __call__(self, r):
                r.headers['Authorization'] = macauthlib.sign_request(r, id=self.id, key=self.key)
                return r



        update_last = time.time() - self.update_time
        if force or update_last > self.update_interval:
            if elapsed_time != None:
                status += '<h3>Time: %02d:%02d:%02d' % (elapsed_time // 3600, elapsed_time % 3600 // 60, elapsed_time % 60)
                if time_limit != None:
                    status += ' / max. %02d:%02d:00 (%d%%)' % (time_limit // 60, time_limit % 60, elapsed_time / time_limit / 60 * 100)
                status += '</h3>'
            self.update_time = self.update_time + update_last // self.update_interval * self.update_interval
            r = requests.put(self.url, status, auth=HTTPMacAuth(self.access_token, self.mac_key))




class StructureEngine(object):

    def __init__(self, db_session, pubchem_names = False, call_back_url = None):
        self.db_session = db_session
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        self.db_session.add(rundata)
        self.db_session.commit()
        self.pubchem_names = pubchem_names
        if pubchem_names:
            self.pubchem_engine = PubChemEngine()
        if call_back_url != None:
            self.call_back_engine = CallBackEngine(call_back_url)
        else:
            self.call_back_engine = None



    def read_sdf(self, file_name, mass_filter):
        print 'READING SDF (Structure Data File)'
        metids = set([])
        for mol in Chem.SDMolSupplier(file_name):
            metids.add(self.add_structure(Chem.MolToMolBlock(mol), mol.GetProp('_Name'), None, 0, 1, mass_filter=mass_filter))

        print str(len(metids)) + ' molecules added to library\n'



    def read_smiles(self, file_name, mass_filter):
        print 'READING SMILES'
        smiles_file = open(file_name)
        metids = set([])
        nonames = 0
        for line in smiles_file:
            line = line.strip()
            if line == '':
                continue
            splitline = line.split()
            smiles = splitline[0]
            if len(splitline) > 1:
                name = '_'.join(splitline[1:])
            else:
                nonames += 1
                name = 'Noname' + str(nonames)
            try:
                mol = Chem.SmilesToMol(smiles, name)
                metids.add(self.add_structure(Chem.MolToMolBlock(mol), name, 1.0, 0, 1, mass_filter=mass_filter))
            except:
                print 'WARNING: Failed to read smiles: ' + smiles + ' (' + name + ')'

        print str(len(metids)) + ' molecules added to library\n'



    def add_structure(self, molblock, name, prob, level, isquery, mim = None, natoms = None, inchikey = None, molform = None, reference = None, logp = None, mass_filter = 9999):
        molecule = types.MoleculeType(molblock, name, prob, level, isquery, mim, natoms, inchikey, molform, reference, logp)
        return self.add_molecule(molecule, mass_filter)



    def add_molecule(self, molecule, mass_filter = 9999, check_duplicates = True, merge = False):
        if molecule.mim > mass_filter:
            return 
        mol = Chem.MolFromMolBlock(molecule.molblock)
        if mol == None:
            return 
        metab = Metabolite(mol=unicode(molecule.molblock, 'utf-8', 'xmlcharrefreplace'), level=molecule.level, probability=molecule.probability, smiles=unicode(molecule.inchikey), molformula=unicode(molecule.molformula), isquery=molecule.isquery, origin=unicode(molecule.name, 'utf-8', 'xmlcharrefreplace'), nhits=0, mim=molecule.mim, natoms=molecule.natoms, reference=unicode(molecule.reference), logp=molecule.logp)
        if check_duplicates:
            dups = self.db_session.query(Metabolite).filter_by(smiles=molecule.inchikey).all()
            if len(dups) > 0:
                if merge:
                    metab = dups[0]
                    if metab.origin == '' and molecule.name != '':
                        metab.origin = unicode(str(metab.origin) + '</br>' + molecule.name, 'utf-8', 'xmlcharrefreplace')
                    if metab.reference == None and molecule.reference != '':
                        metab.reference = unicode(molecule.reference)
                    if molecule.probability > 0:
                        metab.probability = molecule.probability
                elif dups[0].probability < molecule.probability:
                    metab.metid = dups[0].metid
                    self.db_session.delete(dups[0])
                    print 'Duplicate structure, first one removed: ' + molecule.origin
                else:
                    print 'Duplicate structure, kept first one: ' + metab.origin
                    return 
        if self.pubchem_names:
            in_pubchem = self.pubchem_engine.check_inchi(metab.mim, metab.smiles)
            if in_pubchem != False:
                (name, reference,) = in_pubchem
                if metab.origin == '':
                    metab.origin = unicode(name, 'utf-8', 'xmlcharrefreplace')
                if metab.reference == None:
                    metab.reference = unicode(reference)
        self.db_session.add(metab)
        logging.debug('Added molecule: ' + Chem.MolToSmiles(mol))
        self.db_session.flush()
        return metab.metid



    def metabolize(self, metid, metabolism, endpoints = False):
        metabolism_engine = config.get('magma job', 'metabolism_engine')
        if metabolism_engine == 'reactor':
            exec_reactor = pkg_resources.resource_filename('magma', 'script/reactor')
            exec_reactor += ' -s'
            if endpoints:
                exec_reactor += ' -m 10'
            else:
                exec_reactor += ' -a -m 1'
            metabolism_files = {'phase1': pkg_resources.resource_filename('magma', 'data/sygma_rules_GE_0.1.phase1.smirks'),
             'phase2': pkg_resources.resource_filename('magma', 'data/sygma_rules.phase2.smirks'),
             'gut': pkg_resources.resource_filename('magma', 'data/gut.smirks'),
             'digest': pkg_resources.resource_filename('magma', 'data/digest.smirks')}
        elif metabolism_engine == 'cactvs':
            cactvs_root = config.get('magma job', 'cactvs_root')
            exec_reactor = pkg_resources.resource_filename('magma', 'script/csreact')
            exec_reactor += ' ' + cactvs_root
            if endpoints:
                exec_reactor += ' endpoints'
            else:
                exec_reactor += ' parallel'
            metabolism_files = {'phase1': pkg_resources.resource_filename('magma', 'data/sygma_rules.phase1.cactvs.smirks'),
             'phase2': pkg_resources.resource_filename('magma', 'data/sygma_rules.phase2.cactvs.smirks'),
             'phase1_selected': pkg_resources.resource_filename('magma', 'data/sygma_rules.phase1_GE0.05.cactvs.smirks'),
             'phase2_selected': pkg_resources.resource_filename('magma', 'data/sygma_rules.phase2_GE0.05.cactvs.smirks'),
             'gut': pkg_resources.resource_filename('magma', 'data/gut.cactvs.smirks'),
             'glycosidase': pkg_resources.resource_filename('magma', 'data/glycosidase.cactvs.smirks'),
             'peptide': pkg_resources.resource_filename('magma', 'data/peptide.cactvs.smirks'),
             'ptm': pkg_resources.resource_filename('magma', 'data/ptm.cactvs.smirks'),
             'plant': pkg_resources.resource_filename('magma', 'data/plant.cactvs.smirks')}
        try:
            parent = self.db_session.query(Metabolite).filter_by(metid=metid).one()
        except:
            print 'Metabolite record ',
            print metid,
            print ' does not exist.'
            return 
        for m in metabolism.split(','):
            if m in metabolism_files:
                if metabolism_engine == 'reactor':
                    exec_reactor += ' -q'
                exec_reactor += ' ' + metabolism_files[m]

        logging.debug('Execute reactions: ' + exec_reactor)
        reactor = subprocess.Popen(exec_reactor, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        reactor.stdin.write(parent.mol + '$$$$\n')
        reactor.stdin.close()
        metids = set()
        line = reactor.stdout.readline()
        while line != '':
            splitline = line[:-1].split(' {>')
            if len(splitline) == 1:
                reaction = 'PARENT'
                name = line[:-1]
            elif len(splitline) == 2:
                reaction = splitline[1][:-1]
                name = splitline[0]
            else:
                reaction = metabolism
                name = splitline[0]
            mol = name + '\n'
            isquery = 0
            while line != 'M  END\n':
                line = reactor.stdout.readline()
                mol += line

            line = reactor.stdout.readline()
            while line != '$$$$\n' and line != '':
                if line == '> <ReactionSequence>\n':
                    line = reactor.stdout.readline()
                    reaction = line[:-1]
                    line = reactor.stdout.readline()
                    while line != '\n':
                        reaction += ' + ' + line[:-1]
                        line = reactor.stdout.readline()

                line = reactor.stdout.readline()

            if reaction != 'PARENT':
                molecule = types.MoleculeType(mol, '', None, None, isquery)
                new_metid = self.add_molecule(molecule, merge=True)
                metids.add(new_metid)
                reactid = self.db_session.query(Reaction.reactid).filter(Reaction.reactant == metid, Reaction.product == new_metid, Reaction.name == reaction).all()
                if len(reactid) == 0:
                    react = Reaction(reactant=metid, product=new_metid, name=reaction)
                    self.db_session.add(react)
                    self.db_session.flush()
            elif endpoints:
                metids.add(metid)
            line = reactor.stdout.readline()

        self.db_session.add(parent)
        reactor.stdout.close()
        self.db_session.commit()
        if len(metids) == 0:
            metids.add(metid)
        return metids



    def metabolize_all(self, metabolism, endpoints = False):
        print 'Metabolize all'
        parentids = self.db_session.query(Metabolite.metid).all()
        metids = set([])
        for (parentid,) in parentids:
            metids |= self.metabolize(parentid, metabolism, endpoints)

        return metids



    def run_scenario(self, scenario, time_limit = None):
        print 'RUNNING METABOLIC SCENARIO'
        if time_limit == None:
            result = self.db_session.query(Metabolite.metid).all()
            metids = {x[0] for x in result}
        else:
            metids = set([self.db_session.query(Metabolite.metid).first()[0]])
        start_time = time.time()
        for step in range(len(scenario)):
            (action, value,) = scenario[step]
            print 'Stage ' + str(step + 1) + ':'
            endpoints = False
            if action == 'mass_filter':
                result = self.db_session.query(Metabolite.metid).filter(Metabolite.mim < float(value), Metabolite.metid.in_(metids)).all()
                metids = {x[0] for x in result}
                print '   from',
                print len(metids),
                print 'compounds,',
                print len(metids),
                print 'were selected with mass <',
                print value
            else:
                prev_metids = metids
                if value == 'complete':
                    endpoints = True
                    value = 1
                print '  ',
                print len(metids),
                print 'metabolites',
                new_metids = set()
                for metid in metids:
                    new_metids |= self.metabolize(metid, action, endpoints)
                    elapsed_time = time.time() - start_time
                    if self.call_back_engine != None:
                        status = 'Transformation: %s, step 1<br>Metabolites generated: %d' % (action, len(prev_metids) + len(metids) + len(new_metids))
                        self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
                    if time_limit and elapsed_time > time_limit * 60:
                        break
                else:
                    active_metids = new_metids.difference(metids)
                    metids = new_metids
                    for i in range(1, int(value)):
                        new_metids = set()
                        for metid in active_metids:
                            new_metids |= self.metabolize(metid, action, endpoints)
                            elapsed_time = time.time() - start_time
                            if self.call_back_engine != None:
                                status = 'Transformation: %s, step %d<br>Metabolites generated: %d' % (action, i + 1, len(prev_metids) + len(metids) + len(new_metids))
                                self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
                            if time_limit and elapsed_time > time_limit * 60:
                                break

                        active_metids = new_metids.difference(metids)
                        metids |= new_metids
                        if time_limit and time.time() - start_time > time_limit * 60:
                            break


                print 'were metabolized according to',
                print action,
                print 'rules (' + str(int(value)) + ' steps)'
                print '   yielding',
                print len(metids),
                print 'metabolites'
                if not endpoints:
                    metids |= prev_metids
            if time_limit and time.time() - start_time > time_limit * 60:
                if self.call_back_engine != None:
                    self.call_back_engine.update_callback_url('Transformation stopped: time limit exceeded', force=True)
                print 'WARNING: Transformation stopped: time limit exceeded'
                break
        else:
            if self.call_back_engine != None:
                self.call_back_engine.update_callback_url('Transformations completed', force=True)

        print self.db_session.query(Metabolite).count(),
        print 'molecules in library\n'



    def retrieve_structures(self, mass):
        dbfilename = '/home/ridderl/chebi/ChEBI_complete_3star.sqlite'
        conn = sqlite3.connect(dbfilename)
        c = conn.cursor()
        result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ?', (mass - 0.01, mass + 0.01))
        for (id, mim, molblock, smiles, chebi_name,) in result:
            self.add_molecule(zlib.decompress(molblock), str(chebi_name), 1.0, 1, '', 1)




    def get_metids_with_mass_less_then(self, mass, metids = None):
        if metids == None:
            result = self.db_session.query(Metabolite.metid).filter(Metabolite.mim < mass).all()
        else:
            result = self.db_session.query(Metabolite.metid).filter(Metabolite.mim < mass, Metabolite.metid.in_(metids)).all()
        metids = {x[0] for x in result}
        return metids




class MsDataEngine(object):

    def __init__(self, db_session, ionisation_mode, abs_peak_cutoff, mz_precision, mz_precision_abs, precursor_mz_precision, max_ms_level, call_back_url = None):
        self.db_session = db_session
        mz_precision_abs = max(mz_precision_abs, 1e-06)
        precursor_mz_precision = max(precursor_mz_precision, 1e-06)
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.ionisation_mode == None:
            rundata.ionisation_mode = ionisation_mode
        if rundata.abs_peak_cutoff == None:
            rundata.abs_peak_cutoff = abs_peak_cutoff
        if rundata.max_ms_level == None:
            rundata.max_ms_level = max_ms_level
        if rundata.mz_precision == None:
            rundata.mz_precision = mz_precision
        if rundata.mz_precision_abs == None:
            rundata.mz_precision_abs = mz_precision_abs
        if rundata.precursor_mz_precision == None:
            rundata.precursor_mz_precision = precursor_mz_precision
        self.db_session.add(rundata)
        self.db_session.commit()
        if rundata.ionisation_mode == -1:
            self.polarity = '-'
        else:
            self.polarity = '+'
        self.abs_peak_cutoff = rundata.abs_peak_cutoff
        self.max_ms_level = rundata.max_ms_level
        self.mz_precision = rundata.mz_precision * 2
        self.precision = 1 + 2 * rundata.mz_precision / 1000000.0
        self.mz_precision_abs = rundata.mz_precision_abs * 2
        self.precursor_mz_precision = rundata.precursor_mz_precision
        if call_back_url != None:
            self.call_back_engine = CallBackEngine(call_back_url)
        else:
            self.call_back_engine = None



    def store_mzxml_file(self, mzxml_file, scan_filter = None, time_limit = None):
        print 'READING MZXML FILE'
        rundata = self.db_session.query(Run).one()
        if rundata.ms_filename != None:
            print 'ERROR: Attempt to read MS data twice'
            exit()
        rundata.ms_filename = unicode(mzxml_file)
        self.db_session.add(rundata)
        self.ms_filename = mzxml_file
        tree = etree.parse(mzxml_file)
        root = tree.getroot()
        namespace = '{' + root.nsmap[None] + '}'
        mzxml_query = namespace + 'msRun/' + namespace + 'scan'
        prec_scans = []
        start_time = time.time()
        for mzxmlScan in root.findall(mzxml_query):
            elapsed_time = time.time() - start_time
            if mzxmlScan.attrib['polarity'] == self.polarity and (scan_filter == None or mzxmlScan.attrib['num'] == scan_filter or int(mzxmlScan.attrib['msLevel']) > 1 and mzxmlScan.find(namespace + 'precursorMz').attrib['precursorScanNum'] in prec_scans):
                self.store_mzxml_scan(mzxmlScan, 0, namespace)
                prec_scans.append(mzxmlScan.attrib['num'])
                if self.call_back_engine != None:
                    status = 'Reading mzXML, scan: %s' % (prec_scans[-1],)
                    self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
            if time_limit and elapsed_time > time_limit * 60:
                if self.call_back_engine != None:
                    self.call_back_engine.update_callback_url('Reading mzXML stopped: time limit exceeded', force=True)
                print 'WARNING: Reading mzXML stopped: time limit exceeded'
                break
        else:
            if self.call_back_engine != None:
                self.call_back_engine.update_callback_url('Reading mzXML completed', force=True)

        self.db_session.commit()
        print self.db_session.query(Scan).count(),
        print 'spectra read from file\n'



    def store_mzxml_scan(self, mzxmlScan, precScan, namespace):
        if mzxmlScan.attrib['peaksCount'] == '0':
            return 
        try:
            lowmz = float(mzxmlScan.attrib['lowMz'])
            highmz = float(mzxmlScan.attrib['highMz'])
        except:
            lowmz = None
            highmz = None
        try:
            totioncurrent = float(mzxmlScan.attrib['totIonCurrent'])
        except:
            totioncurrent = None
        scan = Scan(scanid=int(mzxmlScan.attrib['num']), mslevel=int(mzxmlScan.attrib['msLevel']), rt=float(mzxmlScan.attrib['retentionTime'].strip('PTS')) / 60, lowmz=lowmz, highmz=highmz, basepeakmz=float(mzxmlScan.attrib['basePeakMz']), basepeakintensity=float(mzxmlScan.attrib['basePeakIntensity']), totioncurrent=totioncurrent, precursorscanid=precScan)
        logging.info('Processing scan ' + str(scan.scanid) + ' (level ' + str(scan.mslevel) + ')')
        comp = None
        for child in mzxmlScan:
            if child.tag == namespace + 'precursorMz':
                scan.precursormz = float(child.text)
                scan.precursorintensity = float(child.attrib['precursorIntensity'])
                if child.attrib.get('precursorScanNum') != None:
                    scan.precursorscanid = float(child.attrib['precursorScanNum'])
                if scan.precursorscanid == 0:
                    scan.precursorscanid = self.db_session.query(Scan.scanid).filter(Scan.mslevel == scan.mslevel - 1).filter(Scan.rt < scan.rt).order_by(desc(Scan.rt)).first()[0]
                    print 'Assigning precursor scanid ' + str(scan.precursorscanid) + ' to scan ' + str(scan.scanid)
                comp = self.db_session.query(Scan).filter(Scan.precursorscanid == scan.precursorscanid, Scan.precursormz == scan.precursormz, Scan.precursorintensity == scan.precursorintensity).all()
            if child.tag == namespace + 'peaks':
                decoded = base64.decodestring(child.text)
                try:
                    if child.attrib['compressionType'] == 'zlib':
                        decoded = zlib.decompress(decoded)
                except:
                    pass
                if comp == None or len(comp) == 0:
                    self.store_mzxml_peaks(scan, decoded)
                else:
                    self.merge_spectrum(comp[0], scan, decoded)
            if child.tag == namespace + 'scan' and int(child.attrib['msLevel']) <= self.max_ms_level:
                self.store_mzxml_scan(child, scan.scanid, namespace)

        if comp == None or len(comp) == 0:
            self.db_session.add(scan)
        self.db_session.flush()



    def merge_spectrum(self, existing_scan, newscan, decoded):
        print 'Merging scans ' + str(existing_scan.scanid) + ' and ' + str(newscan.scanid)
        if existing_scan.lowmz > newscan.lowmz:
            existing_scan.lowmz = newscan.lowmz
        if existing_scan.highmz < newscan.highmz:
            existing_scan.highmz = newscan.highmz
        if existing_scan.basepeakintensity < newscan.basepeakintensity:
            existing_scan.basepeakintensity = newscan.basepeakintensity
            existing_scan.basepeakmz = newscan.basepeakmz
        self.db_session.add(existing_scan)
        tmp_size = len(decoded) / 4
        unpack_format1 = '>%df' % tmp_size
        unpacked = struct.unpack(unpack_format1, decoded)
        for (mz, intensity,) in zip(unpacked[::2], unpacked[1::2]):
            if intensity > self.abs_peak_cutoff:
                matching_peaks = self.db_session.query(Peak).filter(Peak.scanid == existing_scan.scanid, Peak.mz.between(min(mz / self.precision, mz - self.mz_precision_abs), max(mz * self.precision, mz + self.mz_precision_abs))).all()
                if len(matching_peaks) == 0:
                    self.db_session.add(Peak(scanid=existing_scan.scanid, mz=mz, intensity=intensity))
                else:
                    replace = True
                    for p in matching_peaks:
                        if intensity > p.intensity:
                            self.db_session.delete(p)
                        else:
                            replace = False
                            intensity = p.intensity

                    if replace:
                        self.db_session.add(Peak(scanid=existing_scan.scanid, mz=mz, intensity=intensity))




    def store_mzxml_peaks(self, scan, decoded):
        tmp_size = len(decoded) / 4
        unpack_format1 = '>%df' % tmp_size
        unpacked = struct.unpack(unpack_format1, decoded)
        for (mz, intensity,) in zip(unpacked[::2], unpacked[1::2]):
            if intensity > self.abs_peak_cutoff:
                self.db_session.add(Peak(scanid=scan.scanid, mz=mz, intensity=intensity))




    def store_peak_list(self, scanid, rt, precursormz, basepeak, peaklist):
        self.db_session.add(Scan(scanid=scanid, mslevel=1, rt=rt, lowmz=precursormz, highmz=precursormz, basepeakmz=precursormz, basepeakintensity=basepeak[1], precursorscanid=0))
        self.db_session.add(Peak(scanid=scanid, mz=precursormz, intensity=basepeak[1]))
        self.db_session.add(Scan(scanid=scanid + 1, mslevel=2, rt=rt, lowmz=peaklist[0][0], highmz=peaklist[-1][0], basepeakmz=basepeak[0], basepeakintensity=basepeak[1], precursormz=precursormz, precursorintensity=basepeak[1], precursorscanid=scanid))
        for peak in peaklist:
            self.db_session.add(Peak(scanid=scanid + 1, mz=peak[0], intensity=peak[1]))

        self.db_session.commit()



    def store_manual_tree(self, manual_tree, tree_type):
        tree_string = ''.join(open(manual_tree).read().split())
        tree_list = re.split('([\\,\\(\\)])', tree_string)
        self.global_scanid = 1
        self.store_manual_subtree(tree_list, 0, 0, 0, 1, tree_type)



    def store_manual_subtree(self, tree_list, precursor_scanid, precursor_mz, precursor_intensity, mslevel, tree_type):
        lowmz = None
        highmz = None
        basepeakmz = None
        basepeakintensity = None
        scanid = self.global_scanid
        npeaks = 0
        while len(tree_list) > 0 and tree_list[0] != ')':
            tree_item = tree_list.pop(0)
            if tree_item.find(':') >= 0:
                (mz, intensity,) = tree_item.split(':')
                if tree_type != 0:
                    mz = self.mass_from_formula(mz) - tree_type * pars.elmass
                self.db_session.add(Peak(scanid=scanid, mz=mz, intensity=intensity))
                npeaks += 1
                if lowmz == None or mz < lowmz:
                    lowmz = mz
                if highmz == None or mz > highmz:
                    highmz = mz
                if basepeakintensity == None or intensity > basepeakintensity:
                    basepeakmz = mz
                    basepeakintensity = intensity
            elif tree_item == '(':
                self.global_scanid += 1
                self.store_manual_subtree(tree_list, scanid, mz, intensity, mslevel + 1, tree_type)
            elif tree_item != ',' and tree_item != '':
                print 'Corrupt Tree format ...'
                exit()

        if npeaks > 0:
            self.db_session.add(Scan(scanid=scanid, mslevel=mslevel, lowmz=lowmz, highmz=highmz, basepeakmz=basepeakmz, basepeakintensity=basepeakintensity, precursorscanid=precursor_scanid, precursormz=precursor_mz, precursorintensity=precursor_intensity))
        self.db_session.commit()
        if len(tree_list) > 0:
            tree_list.pop(0)



    def mass_from_formula(self, form):
        mass = 0.0
        while len(form) > 0:
            if form[:2] in pars.mims:
                m = pars.mims[form[:2]]
                form = form[2:]
            elif form[:1] in pars.mims:
                m = pars.mims[form[:1]]
                form = form[1:]
            else:
                print 'ERROR: Element not allowed in formula tree: ' + form
                exit()
            x = 0
            while len(form) > x and form[x] in '0123456789':
                x += 1

            if x > 0:
                n = int(form[:x])
            else:
                n = 1
            mass += m * n
            form = form[x:]

        return mass




class AnnotateEngine(object):

    def __init__(self, db_session, skip_fragmentation, max_broken_bonds, max_water_losses, ms_intensity_cutoff, msms_intensity_cutoff, use_all_peaks, adducts = None, force_adduct = False, max_charge = 1, call_back_url = None):
        self.db_session = db_session
        try:
            rundata = self.db_session.query(Run).one()
        except:
            rundata = Run()
        if rundata.skip_fragmentation == None:
            rundata.skip_fragmentation = skip_fragmentation
        if rundata.max_broken_bonds == None:
            rundata.max_broken_bonds = max_broken_bonds
        if rundata.max_water_losses == None:
            rundata.max_water_losses = max_water_losses
        if rundata.ms_intensity_cutoff == None:
            rundata.ms_intensity_cutoff = ms_intensity_cutoff
        if rundata.msms_intensity_cutoff == None:
            rundata.msms_intensity_cutoff = msms_intensity_cutoff
        if rundata.use_all_peaks == None:
            rundata.use_all_peaks = use_all_peaks
        self.db_session.add(rundata)
        self.db_session.commit()
        self.ionisation_mode = rundata.ionisation_mode
        self.skip_fragmentation = rundata.skip_fragmentation
        self.max_broken_bonds = rundata.max_broken_bonds
        self.max_water_losses = rundata.max_water_losses
        self.ms_intensity_cutoff = rundata.ms_intensity_cutoff
        self.msms_intensity_cutoff = rundata.msms_intensity_cutoff
        if rundata.mz_precision == None:
            print 'ERROR: No MS data parameters read.'
            exit()
        self.mz_precision = rundata.mz_precision
        self.precision = 1 + rundata.mz_precision / 1000000.0
        self.mz_precision_abs = rundata.mz_precision_abs
        self.precursor_mz_precision = rundata.precursor_mz_precision
        self.use_all_peaks = rundata.use_all_peaks
        self.scans = []
        if call_back_url != None:
            self.call_back_engine = CallBackEngine(call_back_url)
        else:
            self.call_back_engine = None
        print 'SETTINGS FOR MATCHING PRECURSOR IONS AND CANDIDATES:'
        print 'Ionisation mode:',
        print self.ionisation_mode
        print 'Precursor intensity threshold:',
        print self.ms_intensity_cutoff
        print 'Maximum relative m/z error (ppm):',
        print self.mz_precision
        print 'Maximum absolute m/z error (Da):',
        print self.mz_precision_abs
        if self.ionisation_mode == 1:
            iontypes = ['+H']
        if self.ionisation_mode == -1:
            iontypes = ['-H']
        if adducts != None:
            if force_adduct: iontypes = []
            for i in adducts.split(','):
                if i == 'OH':
                    iontypes.append('-OH')
                else:
                    iontypes.append('+' + i)

        self.ions = self.generate_ions(iontypes, max_charge)



    def generate_ions(self, iontypes, maxcharge):
        ions = [{0: ''}]
        for c in range(0, maxcharge):
            ions.append({})
            for ionmass in ions[c]:
                for i in iontypes:
                    if i not in pars.ionmasses[self.ionisation_mode]:
                        print 'ERROR: Invalid adduct: ' + i + ' for ionisation mode: ' + str(self.ionisation_mode)
                        exit()
                    ions[(c + 1)][ionmass + pars.ionmasses[self.ionisation_mode][i]] = ions[c][ionmass] + i



        for c in range(maxcharge + 1):
            for ionmass in ions[c]:
                ions[c][ionmass] = '[M' + ions[c][ionmass] + ']' + str(c) * (c > 1) + '-' * (self.ionisation_mode < 0) + '+' * (self.ionisation_mode > 0)


        for c in range(1, maxcharge + 1):
            print 'Adducts with charge ' + str(self.ionisation_mode * c) + ': ' + str(ions[c].values())

        print ''
        return ions



    def build_spectrum(self, dbscan):
        scan = types.ScanType(dbscan.scanid, dbscan.mslevel)
        logging.info('Building scan ' + str(dbscan.scanid))
        if scan.mslevel == 1:
            cutoff = self.ms_intensity_cutoff
        else:
            cutoff = dbscan.basepeakintensity * self.msms_intensity_cutoff / 100
        dbpeaks = self.db_session.query(Peak).filter(Peak.scanid == scan.scanid).filter(Peak.intensity >= cutoff).all()
        for dbpeak in dbpeaks:
            scan.peaks.append(types.PeakType(dbpeak.mz, dbpeak.intensity, scan.scanid, pars.missingfragmentpenalty * dbpeak.intensity ** 0.5))

        dbchildscans = self.db_session.query(Scan).filter(Scan.precursorscanid == scan.scanid).all()
        max_depth = 0
        for dbchildscan in dbchildscans:
            prec_intensity = 0.0
            for peak in scan.peaks:
                if peak.intensity > prec_intensity and -self.precursor_mz_precision < dbchildscan.precursormz - peak.mz < self.precursor_mz_precision:
                    prec_peak = peak
                    prec_intensity = peak.intensity

            if prec_intensity == 0.0:
                if dbchildscan.precursorintensity >= cutoff:
                    scan.peaks.append(types.PeakType(dbchildscan.precursormz, dbchildscan.precursorintensity, scan.scanid, pars.missingfragmentpenalty * dbchildscan.precursorintensity ** 0.5))
                    prec_peak = scan.peaks[-1]
                else:
                    continue
            (prec_peak.childscan, depth,) = self.build_spectrum(dbchildscan)
            for childpeak in prec_peak.childscan.peaks:
                prec_peak.missing_fragment_score += childpeak.missing_fragment_score

            max_depth = max(depth, max_depth)

        max_depth += 1
        return (scan, max_depth)



    def build_spectra(self, scans = 'all'):
        print 'BUILDING SPECTRAL TREES'
        ndepths = {}
        if scans == 'all':
            queryscans = self.db_session.query(Scan).filter(Scan.mslevel == 1).all()
        else:
            queryscans = self.db_session.query(Scan).filter(Scan.mslevel == 1).filter(Scan.scanid.in_(scans)).all()
        for dbscan in queryscans:
            (spectrum, depth,) = self.build_spectrum(dbscan)
            if depth > 1:
                if depth in ndepths:
                    ndepths[depth] += 1
                else:
                    ndepths[depth] = 1
            self.scans.append(spectrum)

        print str(len(self.scans)) + ' MS1 spectra'
        for depth in ndepths:
            print str(ndepths[depth]),
            print 'spectral trees of depth',
            print depth

        print ''
        self.indexed_peaks = {}
        for scan in self.scans:
            for peak in scan.peaks:
                if not (not self.use_all_peaks and peak.childscan == None):
                    int_mass = int(round(peak.mz))
                    if int_mass not in self.indexed_peaks:
                        self.indexed_peaks[int_mass] = set([])
                    self.indexed_peaks[int_mass].add(peak)





    def write_tree(self, scanid):
        for scan in self.scans:
            if scan.scanid == scanid:
                for peak in scan.peaks:
                    if not (not self.use_all_peaks and peak.childscan == None):
                        self.write_peak(peak)





    def write_peak(self, peak):
        peak_string = '%.6f: %i' % (peak.mz, peak.intensity)
        if peak.childscan != None:
            peak_string += ' ('
            n = 0
            for childpeak in peak.childscan.peaks:
                if n > 0:
                    peak_string += ', '
                peak_string += self.write_peak(childpeak)
                n += 1

            peak_string += ')'
        return peak_string



    def get_chebi_candidates(self):
        dbfilename = '/home/ridderl/chebi/ChEBI_complete_3star.sqlite'
        conn = sqlite3.connect(dbfilename)
        c = conn.cursor()
        db_candidates = {}
        for scan in self.scans:
            for peak in scan.peaks:
                mass = peak.mz - self.ionisation_mode * pars.Hmass
                if not (not self.use_all_peaks and peak.childscan == None):
                    result = c.execute('SELECT * FROM molecules WHERE mim BETWEEN ? AND ?', (mass / self.precision, mass * self.precision))
                    for (id, mim, molblock, chebi_name,) in result:
                        db_candidates[id] = [molblock, chebi_name]

                    print str(mass) + ' --> ' + str(len(db_candidates)) + ' candidates'


        return db_candidates



    def get_db_candidates(self, query_engine, max_mim = ''):
        print 'RETRIEVING CANDIDATE MOLECULES FROM: ' + str(query_engine.name)
        if max_mim == '':
            max_mim = '1200'
        mmim = int(float(max_mim) * 1000000.0)
        print 'Mass limit:',
        print max_mim
        struct_engine = StructureEngine(self.db_session)
        mzs = []
        for scan in self.scans:
            for peak in scan.peaks:
                if not (not self.use_all_peaks and peak.childscan == None):
                    mzs.append(peak.mz)


        mzs.sort()
        candidate = {}
        for mz in mzs:
            for charge in range(1, len(self.ions)):
                for ionmass in self.ions[charge]:
                    ql = int(1000000.0 * ((min(mz / self.precision, mz - self.mz_precision_abs) - ionmass) * charge + self.ionisation_mode * pars.elmass))
                    qh = int(1000000.0 * ((max(mz * self.precision, mz + self.mz_precision_abs) - ionmass) * charge + self.ionisation_mode * pars.elmass))
                    result = query_engine.query_on_mim(ql, qh, 0)
                    for molecule in result:
                        candidate[molecule.inchikey] = molecule

                    logging.info(str(ql) + ',' + str(qh) + ' --> ' + str(len(candidate)) + ' candidates')

                for ionmass in self.ions[(charge - 1)]:
                    ql = int(1000000.0 * ((min(mz / self.precision, mz - self.mz_precision_abs) - ionmass) * charge + self.ionisation_mode * pars.elmass))
                    qh = int(1000000.0 * ((max(mz * self.precision, mz + self.mz_precision_abs) - ionmass) * charge + self.ionisation_mode * pars.elmass))
                    result = query_engine.query_on_mim(ql, qh, self.ionisation_mode)
                    for molecule in result:
                        candidate[molecule.inchikey] = molecule

                    logging.info(str(ql) + ',' + str(qh) + ' --> ' + str(len(candidate)) + ' candidates')



        check_duplicates = self.db_session.query(Metabolite.metid).count() > 0
        logging.debug('check duplicates: ' + str(check_duplicates))
        metids = set([])
        for molecule in candidate.itervalues():
            metid = struct_engine.add_molecule(molecule, check_duplicates=check_duplicates, merge=True)
            metids.add(metid)

        self.db_session.commit()
        print str(len(metids)) + ' candidates retrieved\n'
        return metids



    def search_structures(self, metids = None, fast = False, time_limit = None):
        global fragid
        print 'MATCHING CANDIDATE MOLECULES'
        fragid = self.db_session.query(func.max(Fragment.fragid)).scalar()
        if fragid == None:
            fragid = 0
        if metids == None:
            if time_limit == None:
                metabdata = self.db_session.query(Metabolite.metid).order_by(desc(Metabolite.metid)).all()
            else:
                metabdata = self.db_session.query(Metabolite.metid).order_by(Metabolite.probability).all()
            metids = [ x[0] for x in metabdata ]
        total_frags = 0
        total_metids = len(metids)
        start_time = time.time()
        count = 0
        while len(metids) > 0:
            ids = set([])
            while len(ids) < 500 and len(metids) > 0:
                ids.add(metids.pop())

            structures = self.db_session.query(Metabolite).filter(Metabolite.metid.in_(ids)).all()
            jobs = []
            for structure in structures:
                if self.db_session.query(Fragment.fragid).filter(Fragment.metid == structure.metid).count() > 0:
                    logging.warn('Metabolite ' + str(structure.metid) + ': Already annotated, skipped')
                    continue
                molcharge = 0
                molcharge += 1 * (structure.molformula[-1] == '-' and self.ionisation_mode == -1 or structure.molformula[-1] == '+' and self.ionisation_mode == 1)
                peaks = set([])
                for charge in range(1, len(self.ions)):
                    for ionmass in self.ions[(charge - molcharge)]:
                        int_mass = int(round((structure.mim + ionmass) / charge))
                        try:
                            peaks = peaks.union(self.indexed_peaks[int_mass])
                        except:
                            pass
                        try:
                            peaks = peaks.union(self.indexed_peaks[(int_mass - 1)])
                        except:
                            pass
                        try:
                            peaks = peaks.union(self.indexed_peaks[(int_mass + 1)])
                        except:
                            pass


                if len(peaks) == 0:
                    logging.info('Metabolite ' + str(structure.metid) + ': No match')
                    continue
                (hits, frags,) = search_structure(structure.mol, structure.mim, molcharge, peaks, self.max_broken_bonds, self.max_water_losses, self.precision, self.mz_precision_abs, self.use_all_peaks, self.ionisation_mode, self.skip_fragmentation, fast and structure.natoms <= 64, config.get('magma job', 'chemical_engine'), self.ions)
                total_frags += frags
                logging.debug(' -> ' + str(frags) + ' fragments')
                structure.nhits = len(hits)
                self.db_session.add(structure)
                if len(hits) == 0:
                    logging.info('Metabolite ' + str(structure.metid) + ': No match')
                else:
                    print 'Metabolite',
                    print str(structure.metid) + ':',
                    print structure.origin.encode('utf-8')
                    for hit in hits:
                        print 'Scan: ' + str(hit.scan) + ' - Mz: ' + str(hit.mz) + ' - ' + 'Score:',
                        print self.store_hit(hit, structure.metid, 0)

                self.db_session.flush()
                count += 1
                elapsed_time = time.time() - start_time
                if self.call_back_engine != None:
                    status = 'Annotation: %d / %d candidate molecules processed  (%d%%)' % (total_metids - len(metids) - len(ids) + count, total_metids, 100.0 * (total_metids - len(metids) - len(ids) + count) / total_metids)
                    self.call_back_engine.update_callback_url(status, elapsed_time, time_limit)
                if time_limit and elapsed_time > time_limit * 60:
                    metids = []
                    if self.call_back_engine != None:
                        self.call_back_engine.update_callback_url('Annotation stopped: time limit exceeded', force=True)
                    print 'WARNING: Annotation stopped: time limit exceeded'
                    break

            self.db_session.commit()

        if self.call_back_engine != None:
            self.call_back_engine.update_callback_url('Annotation completed', force=True)
        logging.info(str(total_frags) + ' fragments generated in total.')
        print self.db_session.query(Fragment.metid).filter(Fragment.parentfragid == 0).distinct().count(),
        print 'Molecules matched with',
        print self.db_session.query(Fragment.scanid).filter(Fragment.parentfragid == 0).distinct().count(),
        print 'precursor ions, in total\n'



    def store_hit(self, hit, metid, parentfragid):
        global fragid
        fragid += 1
        currentFragid = fragid
        score = hit.score
        deltappm = None
        if score != None:
            score = score / hit.intensity_weight
            charge = 1
            if hit.ion[-2] in '123456789':
                charge = int(hit.ion[-2])
            deltappm = (hit.mz - (hit.mass + hit.deltaH) / charge + self.ionisation_mode * pars.elmass) / hit.mz * 1000000.0
        self.db_session.add(Fragment(metid=metid, scanid=hit.scan, mz=hit.mz, mass=hit.mass, score=score, parentfragid=parentfragid, atoms=unicode(hit.atomstring), inchikey=unicode(hit.inchikey), deltah=hit.deltaH, deltappm=deltappm, formula=unicode(hit.formula + '<br>' + hit.ion)))
        if len(hit.besthits) > 0:
            for childhit in hit.besthits:
                if childhit != None:
                    self.store_hit(childhit, metid, currentFragid)

        return score




class PubChemEngine(object):

    def __init__(self, dbfilename = '', max_64atoms = False, incl_halo = '', min_refscore = ''):
        self.name = 'PubChem'
        if dbfilename == '':
            dbfilename = config.get('magma job', 'structure_database.pubchem')
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory = str
        self.c = self.conn.cursor()
        self.incl_halo = False
        if incl_halo != '' and incl_halo != 'False':
            self.incl_halo = True
            if incl_halo == 'True':
                halo_filename = config.get('magma job', 'structure_database.pubchem_halo')
            else:
                halo_filename = incl_halo
            self.connh = sqlite3.connect(halo_filename)
            self.connh.text_factory = str
            self.ch = self.connh.cursor()
        self.where = ''
        if min_refscore != '':
            self.where += ' AND refscore >= ' + min_refscore
        if max_64atoms == True:
            self.where += ' AND natoms <= 64'



    def query_on_mim(self, low, high, charge):
        molecules = []
        result = self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
        self.add_result2molecules(result, molecules)
        if self.incl_halo:
            result = self.ch.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
            self.add_result2molecules(result, molecules)
        return molecules



    def add_result2molecules(self, result, molecules):
        for (cid, mim, charge, natoms, molblock, inchikey, molform, name, refscore, logp,) in result:
            molecules.append(types.MoleculeType(molblock=zlib.decompress(molblock), name=name + ' (' + str(cid) + ')', mim=float(mim / 1000000.0), natoms=natoms, molform=molform, inchikey=inchikey, prob=refscore, level=1, isquery=1, reference='<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&cmd=Link&LinkName=pccompound_pccompound_sameisotopic_pulldown&from_uid=' + str(cid) + '" target="_blank">' + str(cid) + ' (PubChem)</a>', logp=float(logp) / 10.0))




    def check_inchi(self, mim, inchikey):
        self.c.execute('SELECT cid,name FROM molecules WHERE charge IN (-1,0,1) AND mim between ? and ? and inchikey = ?', (int(mim * 1000000.0) - 1, int(mim * 1000000.0) + 1, inchikey))
        result = self.c.fetchall()
        if len(result) > 0:
            (cid, name,) = result[0]
            reference = '<a href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=pccompound&cmd=Link&LinkName=pccompound_pccompound_sameisotopic_pulldown&from_uid=' + str(cid) + '" target="_blank">' + str(cid) + ' (PubChem)</a>'
            return [name, reference]
        else:
            return False




class KeggEngine(object):

    def __init__(self, dbfilename = '', max_64atoms = False, incl_halo = ''):
        self.name = 'Kegg'
        if dbfilename == '':
            dbfilename = config.get('magma job', 'structure_database.kegg')
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory = str
        self.c = self.conn.cursor()
        self.incl_halo = False
        if incl_halo != '' and incl_halo != 'False':
            self.incl_halo = True
            if incl_halo == 'True':
                halo_filename = config.get('magma job', 'structure_database.kegg_halo')
            else:
                halo_filename = incl_halo
            self.connh = sqlite3.connect(halo_filename)
            self.connh.text_factory = str
            self.ch = self.connh.cursor()
        self.where = ''
        if max_64atoms == True:
            self.where += ' AND natoms <= 64'



    def query_on_mim(self, low, high, charge):
        molecules = []
        result = self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
        self.add_result2molecules(result, molecules)
        if self.incl_halo:
            result = self.ch.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
            self.add_result2molecules(result, molecules)
        return molecules



    def add_result2molecules(self, result, molecules):
        for (cid, mim, charge, natoms, molblock, inchikey, molform, name, reference, logp,) in result:
            keggids = reference.split(',')
            keggrefs = '<a href="http://www.genome.jp/dbget-bin/www_bget?cpd:' + keggids[0] + '" target="_blank">' + keggids[0] + ' (Kegg)</a>'
            for keggid in keggids[1:]:
                keggrefs += '<br><a href="http://www.genome.jp/dbget-bin/www_bget?cpd:' + keggid + '" target="_blank">' + keggid + ' (Kegg)</a>'

            molecules.append(types.MoleculeType(molblock=zlib.decompress(molblock), name=name + ' (' + str(cid) + ')', mim=float(mim / 1000000.0), natoms=natoms, molform=molform, inchikey=inchikey, prob=None, level=1, isquery=1, reference=keggrefs, logp=float(logp) / 10.0))





class HmdbEngine(object):

    def __init__(self, dbfilename = '', max_64atoms = False):
        self.name = 'Human Metabolite Database'
        if dbfilename == '':
            dbfilename = config.get('magma job', 'structure_database.hmdb')
        self.where = ''
        if max_64atoms == True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory = str
        self.c = self.conn.cursor()



    def query_on_mim(self, low, high, charge):
        result = self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
        molecules = []
        for (cid, mim, charge, natoms, molblock, inchikey, molform, name, reference, logp,) in result:
            hmdb_ids = reference.split(',')
            hmdb_refs = '<a href="http://www.hmdb.ca/metabolites/' + hmdb_ids[0] + '" target="_blank">' + hmdb_ids[0] + ' (HMDB)</a>'
            for hmdb_id in hmdb_ids[1:]:
                hmdb_refs += '<br><a href="http://www.hmdb.ca/metabolites/' + hmdb_id + '" target="_blank">' + hmdb_id + ' (HMDB)</a>'

            molecules.append(types.MoleculeType(molblock=zlib.decompress(molblock), name=name + ' (' + str(cid) + ')', mim=float(mim / 1000000.0), natoms=natoms, molform=molform, inchikey=inchikey, prob=None, level=1, isquery=1, reference=hmdb_refs, logp=float(logp) / 10.0))

        return molecules



class MetlinEngine(object):

    def __init__(self, dbfilename = '', max_64atoms = False):
        self.name = 'Metlin Database'
        if dbfilename == '':
            dbfilename = config.get('magma job', 'structure_database.metlin')
        self.where = ''
        if max_64atoms == True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory = str
        self.c = self.conn.cursor()



    def query_on_mim(self, low, high, charge):
        result = self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
        molecules = []
        for (mid, mim, charge, natoms, molblock, inchikey, molform, name, reference, logp,) in result:
            metlin_id = reference
            metlin_ref = '<a href="http://metlin.scripps.edu/metabo_info.php?molid=' + metlin_id + '</a>'
            molecules.append(types.MoleculeType(molblock=zlib.decompress(molblock), name=name + ' (' + str(mid) + ')', mim=float(mim / 1000000.0), natoms=natoms, molform=molform, inchikey=inchikey, prob=None, level=1, isquery=1, reference=metlin_ref, logp=float(logp) / 10.0))

        return molecules



class MetaCycEngine(object):

    def __init__(self, dbfilename = '', max_64atoms = False):
        self.name = 'MetaCyc'
        if dbfilename == '':
            dbfilename = config.get('magma job', 'structure_database.metacyc')
        self.where = ''
        if max_64atoms == True:
            self.where += ' AND natoms <= 64'
        self.conn = sqlite3.connect(dbfilename)
        self.conn.text_factory = str
        self.c = self.conn.cursor()



    def query_on_mim(self, low, high, charge):
        result = self.c.execute('SELECT * FROM molecules WHERE charge = ? AND mim BETWEEN ? AND ? %s' % self.where, (charge, low, high))
        molecules = []
        for (cid, mim, charge, natoms, molblock, inchikey, molform, name, reference, logp,) in result:
            metacyc_ids = reference.split(',')
            metacyc_refs = '<a href="http://www.biocyc.org/META/NEW-IMAGE?type=COMPOUND&object=' + metacyc_ids[0] + '" target="_blank">' + metacyc_ids[0] + ' (MetaCyc)</a>'
            for metacyc_id in metacyc_ids[1:]:
                metacyc_refs += '<br><a href="http://www.biocyc.org/META/NEW-IMAGE?type=COMPOUND&object=' + metacyc_id + '" target="_blank">' + metacyc_id + ' (MetaCyc)</a>'

            molecules.append(types.MoleculeType(molblock=zlib.decompress(molblock), name=name, mim=float(mim / 1000000.0), natoms=natoms, molform=molform, inchikey=inchikey, prob=None, level=1, isquery=1, reference=metacyc_refs, logp=float(logp) / 10.0))

        return molecules




class SelectEngine(object):

    def __init__(self, db_session):
        self.db_session = db_session
        (mz_precision, mz_precision_abs,) = self.db_session.query(Run.mz_precision, Run.mz_precision_abs).all()[0]
        self.precision = 1 + mz_precision / 1000000.0
        self.mz_precision_abs = mz_precision_abs



    def select_fragment(self, fragid):
        child_frags = self.db_session.query(Fragment.fragid).filter(Fragment.parentfragid == fragid).all()
        if len(child_frags) == 0:
            exit('Fragment not fragmented')
        fragmented_fragids = self.db_session.query(distinct(Fragment.parentfragid))
        smiles = self.db_session.query(Fragment.inchikey).filter(Fragment.fragid == fragid).all()
        sys.stderr.write(str(smiles) + '\n')
        metids = self.db_session.query(Fragment.metid).filter(Fragment.inchikey == smiles[0][0], Fragment.fragid.in_(fragmented_fragids))
        self.db_session.query(Fragment).filter(~Fragment.metid.in_(metids)).delete(synchronize_session='fetch')
        self.db_session.query(Metabolite).filter(~Metabolite.metid.in_(metids)).delete(synchronize_session='fetch')
        self.db_session.commit()
        ref_scan = self.db_session.query(distinct(Fragment.scanid)).filter(Fragment.parentfragid == fragid).all()[0][0]
        fragids = self.db_session.query(Fragment.fragid).filter(Fragment.inchikey == smiles[0][0], Fragment.fragid.in_(fragmented_fragids))
        query_scans = self.db_session.query(distinct(Fragment.scanid)).filter(Fragment.parentfragid.in_(fragids)).all()
        for (query_scan,) in query_scans:
            dot_product = self.dot_product_scans(ref_scan, query_scan)
            compounds = self.db_session.query(Metabolite).filter(Metabolite.metid.in_(self.db_session.query(Fragment.metid).filter(Fragment.scanid == query_scan))).all()
            self.db_session.add(compound)

        self.db_session.commit()



    def dot_product_scans(self, ref_scan, query_scan):
        ref_peaks = self.db_session.query(Peak.mz, Peak.intensity).filter(Peak.scanid == ref_scan).all()
        query_peaks = self.db_session.query(Peak.mz, Peak.intensity).filter(Peak.scanid == query_scan).all()
        ref_peaks_sorted = sorted(ref_peaks, key=itemgetter(1), reverse=True)
        query_peaks_sorted = sorted(query_peaks, key=itemgetter(1), reverse=True)
        dot_product = None
        if len(query_peaks) > 0:
            matched_intensities = []
            for pr in ref_peaks_sorted:
                ll = min(pr[0] / self.precision, pr[0] - self.mz_precision_abs)
                hl = max(pr[0] * self.precision, pr[0] + self.mz_precision_abs)
                min_delta_intensity = None
                match = None
                for pq in query_peaks_sorted:
                    if ll < pq[0] < hl:
                        delta_intensity = abs(pq[1] - pr[1])
                        if min_delta_intensity == None or min_delta_intensity > delta_intensity:
                            min_delta_intensity = delta_intensity
                            match = pq

                if match == None:
                    matched_intensities.append((pr[1], 0))
                else:
                    matched_intensities.append((pr[1], match[1]))
                    query_peaks_sorted.remove(match)

            for pq in query_peaks_sorted:
                matched_intensities.append((0, pq[1]))

            srq = 0
            sr = 0
            sq = 0
            for match in matched_intensities:
                srq += (match[0] * match[1]) ** 0.5
                sr += match[0]
                sq += match[1]

            dot_product = srq ** 2 / (sr * sq)
        return dot_product




class DataAnalysisEngine(object):

    def __init__(self, db_session):
        self.db_session = db_session



    def get_scores(self, scanid):
        return self.db_session.query(Fragment.score, Metabolite.smiles, Metabolite.molformula).join((Metabolite, and_(Fragment.metid == Metabolite.metid))).filter(Fragment.parentfragid == 0).filter(Fragment.scanid == scanid).all()



    def get_num_peaks(self, scanid):
        return self.db_session.query(Peak).filter(Peak.scanid == scanid).count()



    def export_assigned_molecules(self, name):
        for (metabolite, peak,) in self.db_session.query(Metabolite, Peak).filter(Metabolite.metid == Peak.assigned_metid):
            print metabolite.origin.splitlines()[0]
            print metabolite.mol[(metabolite.mol.find('\n') + 1):-1]
            print '> <ScanID>\n' + str(peak.scanid) + '\n'
            print '> <mz>\n' + str(peak.mz) + '\n'
            print '> <intensity>\n' + str(peak.intensity) + '\n'
            print '> <rt>\n' + str(self.db_session.query(Scan.rt).filter(Scan.scanid == peak.scanid).all()[0][0]) + '\n'
            print '> <molecular formula>\n' + metabolite.molformula + '\n'
            print '$$$$'




    def write_SDF(self, file = sys.stdout, molecules = None, columns = None, sortcolumn = None, descend = False):
        if molecules == None:
            if descend:
                molecules = self.db_session.query(Metabolite).order_by(desc(sortcolumn)).all()
            else:
                molecules = self.db_session.query(Metabolite).order_by(sortcolumn).all()
        for molecule in molecules:
            file.write(molecule.mol)
            if columns == None:
                columns = dir(molecule)
            for column in columns:
                if column[:1] != '_' and column != 'mol' and column != 'metadata':
                    file.write('> <' + column + '>\n' + str(molecule.__getattribute__(column)) + '\n\n')

            file.write('$$$$\n')




    def write_smiles(self, file = sys.stdout, molecules = None, columns = None, sortcolumn = None, descend = False):
        if molecules == None:
            if descend:
                molecules = self.db_session.query(Metabolite).order_by(desc(sortcolumn)).all()
            else:
                molecules = self.db_session.query(Metabolite).order_by(sortcolumn).all()
        for molecule in molecules:
            file.write(str(molecule.origin).split()[-1][1:-1] + '_' + str(molecule.smiles) + ' ')
            file.write(Chem.MolToSmiles(Chem.MolFromMolBlock(str(molecule.mol))) + '\n')




    def write_ranked_list(self, scan = 1):
        print 'Candidate_Score Name Smiles'
        for (name, smiles, score,) in self.db_session.query(Metabolite.origin, Fragment.inchikey, Fragment.score).join(Fragment).filter(Fragment.scanid == scan).order_by(Fragment.score).all():
            print score,
            print name.encode('utf-8'),
            print smiles




    def write_network1(self, filename):
        f = open(filename + '.sif', 'w')
        assigned_metids = self.db_session.query(distinct(Peak.assigned_metid))
        written_metids = set([])
        for (reactant, product,) in self.db_session.query(Reaction.reactant, Reaction.product).filter(Reaction.reactant.in_(assigned_metids) | Reaction.product.in_(assigned_metids)).all():
            f.write(str(reactant) + ' pp ' + str(product) + '\n')
            written_metids.add(reactant)
            written_metids.add(product)

        f.close()
        f = open(filename + '.txt', 'w')
        assigned = assigned_metids.all()
        for metid in written_metids:
            if (metid,) in assigned:
                f.write(str(metid) + ' assigned\n')
            else:
                f.write(str(metid) + ' unassigned\n')




    def write_network2(self, filename):
        nodes = {}
        for (reactant, product,) in self.db_session.query(Reaction.reactant, Reaction.product).all():
            r = int(reactant)
            p = int(product)
            if r not in nodes:
                nodes[r] = set([])
            if p not in nodes:
                nodes[p] = set([])
            nodes[r].add(p)
            nodes[p].add(r)

        result = self.db_session.query(distinct(Peak.assigned_metid)).all()
        assigned_metids = [ x[0] for x in result ]
        start_compound = {}
        result = self.db_session.query(Metabolite.metid, Metabolite.isquery).all()
        for (metid, isquery,) in result:
            start_compound[metid] = isquery

        print result
        print assigned_metids
        nnodes = len(nodes) + 1
        while len(nodes) < nnodes:
            nnodes = len(nodes)
            print nnodes
            nodekeys = nodes.keys()
            for n in nodekeys:
                if n not in assigned_metids and start_compound[n] == False:
                    if len(nodes[n]) == 1:
                        nodes[list(nodes[n])[0]].remove(n)
                        del nodes[n]
                    elif len(nodes[n]) == 2:
                        tmpnode = list(nodes[n])
                        if len(nodes[tmpnode[0]] & nodes[tmpnode[1]]) > 1 or len(nodes[tmpnode[0]] & nodes[n]) > 0 or len(nodes[tmpnode[1]] & nodes[n]) > 0:
                            for c in nodes[n]:
                                nodes[c].remove(n)

                            del nodes[n]


        f = open(filename + '.sif', 'w')
        written_metids = set([])
        connections = []
        for n in nodes:
            for c in nodes[n]:
                l = set([n, c])
                if l not in connections:
                    f.write(str(n) + ' pp ' + str(c) + '\n')
                    connections.append(l)

            written_metids.add(n)

        f.close()
        f = open(filename + '.txt', 'w')
        for metid in written_metids:
            if metid in assigned_metids:
                f.write(str(metid) + ' ' + str(start_compound[metid] + 2) + '\n')
            else:
                f.write(str(metid) + ' ' + str(start_compound[metid] + 0) + '\n')





def search_structure(mol, mim, molcharge, peaks, max_broken_bonds, max_water_losses, precision, mz_precision_abs, use_all_peaks, ionisation_mode, skip_fragmentation, fast, chem_engine, ions):
    if fast:
        import fragmentation_cy as Fragmentation
    else:
        import fragmentation_py as Fragmentation

    def massmatch(peak, mim, molcharge):
        lowmz = min(peak.mz / precision, peak.mz - mz_precision_abs)
        highmz = max(peak.mz * precision, peak.mz + mz_precision_abs)
        for charge in range(1, len(ions)):
            for ionmass in ions[(charge - molcharge)]:
                if lowmz <= (mim + ionmass) / charge - ionisation_mode * pars.elmass <= highmz:
                    return [ionmass, ions[(charge - molcharge)][ionmass]]


        return False



    def gethit(peak, fragment, score, bondbreaks, mass, ionmass, ion):
        try:
            hit = types.HitType(peak, fragment, score, bondbreaks, mass, ionmass, ion)
        except:
            hit = magma.types.HitType(peak, fragment, score, bondbreaks, mass, ionmass, ion)
        if fragment > 0 and peak.childscan != None and len(peak.childscan.peaks) > 0:
            n_child_peaks = len(peak.childscan.peaks)
            total_score = 0.0
            total_count = 0.0
            for childpeak in peak.childscan.peaks:
                besthit = gethit(childpeak, 0, None, 0, 0, 0, '')
                mz_neutral = childpeak.mz + ionisation_mode * pars.elmass
                for (childfrag, childscore, childbbreaks, childmass, childH,) in fragment_engine.find_fragments(mz_neutral, fragment, precision, mz_precision_abs):
                    if childfrag & fragment == childfrag:
                        childhit = gethit(childpeak, childfrag, childscore * childpeak.intensity ** 0.5, childbbreaks, childmass, childH * pars.Hmass, '[X' + '+' * (childH > 0) + '-' * (childH < 0) + str(abs(childH)) * (not -2 < childH < 2) + 'H' * (childH != 0) + ']' + '+' * (ionisation_mode > 0) + '-' * (ionisation_mode < 0))
                        if besthit.score == None or besthit.score > childhit.score or besthit.score == childhit.score and abs(besthit.deltaH) > abs(childhit.deltaH) or fragment_engine.score_fragment_rel2parent(besthit.fragment, fragment) > fragment_engine.score_fragment_rel2parent(childhit.fragment, fragment):
                            besthit = childhit

                if besthit.score == None:
                    total_score += childpeak.missing_fragment_score
                else:
                    hit.besthits.append(besthit)
                    total_score += min(besthit.score, childpeak.missing_fragment_score)

            hit.score = hit.score + total_score
        return hit



    def add_fragment_data_to_hit(hit):
        if hit.fragment != 0:
            (hit.atomstring, hit.atomlist, hit.formula, hit.inchikey,) = fragment_engine.get_fragment_info(hit.fragment, hit.deltaH)
            if len(hit.besthits) > 0:
                for childhit in hit.besthits:
                    if childhit != None:
                        add_fragment_data_to_hit(childhit)



    Fragmented = False
    hits = []
    frags = 0
    for peak in peaks:
        if not (not use_all_peaks and peak.childscan == None):
            i = massmatch(peak, mim, molcharge)
            if i != False:
                if not Fragmented:
                    fragment_engine = Fragmentation.FragmentEngine(mol, max_broken_bonds, max_water_losses, ionisation_mode, skip_fragmentation, molcharge)
                    if fragment_engine.accepted():
                        frags = fragment_engine.generate_fragments()
                    Fragmented = True
                if fragment_engine.accepted():
                    hit = gethit(peak, (1 << fragment_engine.get_natoms()) - 1, 0, 0, mim, i[0], i[1])
                    add_fragment_data_to_hit(hit)
                    hits.append(hit)

    return (hits, frags)
