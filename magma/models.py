# 2014.10.20 12:28:51 CEST
"""
Sqlalchemy models for magma result database
"""
import json
from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import Unicode
from sqlalchemy import Float
from sqlalchemy import Boolean
from sqlalchemy import ForeignKey
from sqlalchemy import TypeDecorator
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.schema import ForeignKeyConstraint
Base = declarative_base()

class ReactionSequence(TypeDecorator):
    """List of reactions.
    
        Reactions are grouped by relation to the current row or molecule:
    
        * products, reactions which have current row as reactant
        * reactants, reactions which have current row as product
    
        Reactions is a dict with key as the reaction name and the value a dict with keys:
    
        * nr, number of molecules which are product/reactant
            of current molecule with this reaction
        * nrp, number of molecules which are product/reactant
            of current molecule with this reaction and which have been matched to at least one scan
    
        Example:
    
        .. code-block:: javascript
    
              {
                'reactants': {
                   'esterase': {'nr': 123, 'nrp': 45}
                },
                'products': {
                   'sulfation': {'nr': 678, 'nrp': 90}
                }
              }
    
        Stored in database as json serialized string.
        """

    impl = Unicode

    def process_bind_param(self, value, dialect):
        if value is not None:
            value = json.dumps(value)
        return value



    def process_result_value(self, value, dialect):
        if value is u'':
            value = u'{}'
        if value is not None:
            try:
                value = json.loads(value)
            except ValueError:
                value = {}
        return value




class Metabolite(Base):
    """Metabolite model for metabolites table"""

    __tablename__ = 'metabolites'
    metid = Column(Integer, primary_key=True, autoincrement=True)
    mol = Column(Unicode)
    level = Column(Integer)
    probability = Column(Float)
    reactionsequence = Column(ReactionSequence, default={})
    smiles = Column(Unicode, unique=True)
    molformula = Column(Unicode)
    isquery = Column(Boolean)
    origin = Column(Unicode)
    nhits = Column(Integer)
    mim = Column(Float)
    natoms = Column(Integer)
    logp = Column(Float)
    reference = Column(Unicode)
    fragments = relationship('Fragment', backref='metabolite')


class Reaction(Base):
    """Reaction model for reactions table"""

    __tablename__ = 'reactions'
    reactid = Column(Integer, primary_key=True, autoincrement=True)
    reactant = Column(Integer, ForeignKey('metabolites.metid'))
    product = Column(Integer, ForeignKey('metabolites.metid'))
    name = Column(Unicode)


def fill_molecules_reactions(session):
    """Fills the reactionsequence column in the molecules table with info from reactions table.
    
        The molecules query will become to complex when reactionsequence is queried from reactions table.
        So we fill reaction sequence with a json serialized struct which can be used during rendering.
    
        from magmaweb.job import JobFactory
        factory = JobFactory('data/jobs')
        session = factory._makeJobSession('58f05077-aad8-4fc9-a497-310495ab8b62')
        from magmaweb.models import fill_molecules_reactionsequence
        fill_molecules_reactionsequence(session)
    
        """
    from sqlalchemy.sql import func
    reactions = {}
    for (metid, rname, nr,) in session.query(Reaction.product, Reaction.name, func.count('*')).group_by(Reaction.product, Reaction.name):
        if metid not in reactions:
            reactions[metid] = {}
        if 'reactants' not in reactions[metid]:
            reactions[metid]['reactants'] = {}
        reactions[metid]['reactants'][rname] = {'nr': nr}

    for (metid, rname, nrp,) in session.query(Reaction.product, Reaction.name, func.count('*')).join(Metabolite, Metabolite.metid == Reaction.reactant).filter(Metabolite.nhits > 0).group_by(Reaction.product, Reaction.name):
        reactions[metid]['reactants'][rname]['nrp'] = nrp

    for (metid, rname, nr,) in session.query(Reaction.reactant, Reaction.name, func.count('*')).group_by(Reaction.reactant, Reaction.name):
        if metid not in reactions:
            reactions[metid] = {}
        if 'products' not in reactions[metid]:
            reactions[metid]['products'] = {}
        reactions[metid]['products'][rname] = {'nr': nr}

    for (metid, rname, nrp,) in session.query(Reaction.reactant, Reaction.name, func.count('*')).join(Metabolite, Metabolite.metid == Reaction.product).filter(Metabolite.nhits > 0).group_by(Reaction.reactant, Reaction.name):
        reactions[metid]['products'][rname]['nrp'] = nrp

    for mol in session.query(Metabolite):
        if mol.metid in reactions:
            reaction = reactions[mol.metid]
        else:
            reaction = {}
        mol.reactionsequence = reaction

    session.commit()



class Scan(Base):
    """Scan model for scans table"""

    __tablename__ = 'scans'
    scanid = Column(Integer, primary_key=True)
    mslevel = Column(Integer)
    rt = Column(Float)
    lowmz = Column(Float)
    highmz = Column(Float)
    basepeakmz = Column(Float)
    basepeakintensity = Column(Float)
    totioncurrent = Column(Float)
    precursormz = Column(Float)
    precursorintensity = Column(Float)
    precursorscanid = Column(Integer, ForeignKey('scans.scanid'))
    peaks = relationship('Peak', backref='scan')
    products = relationship('Scan', backref=backref('precursor', remote_side=[scanid]))
    fragments = relationship('Fragment', backref='scan')


class Peak(Base):
    """Peak model for peaks table"""

    __tablename__ = 'peaks'
    scanid = Column(Integer, ForeignKey('scans.scanid'), primary_key=True)
    mz = Column(Float, primary_key=True)
    intensity = Column(Float)
    assigned_metid = Column(Integer, ForeignKey('metabolites.metid'), index=True)


class Fragment(Base):
    """Fragment model for fragments table"""

    __tablename__ = 'fragments'
    fragid = Column(Integer, primary_key=True, autoincrement=True)
    metid = Column(Integer, ForeignKey('metabolites.metid'), index=True)
    scanid = Column(Integer, ForeignKey('scans.scanid'), index=True)
    mz = Column(Float)
    mass = Column(Float)
    score = Column(Float)
    parentfragid = Column(Integer, ForeignKey('fragments.fragid'))
    atoms = Column(Unicode)
    deltah = Column(Float)
    deltappm = Column(Float)
    inchikey = Column(Unicode)
    formula = Column(Unicode)
    children_backref = backref('parent', remote_side=[fragid])
    children = relationship('Fragment', backref=children_backref, lazy='joined', join_depth=1)
    __table_args__ = (ForeignKeyConstraint(['scanid', 'mz'], ['peaks.scanid', 'peaks.mz']), {})


class Run(Base):
    """Run model for run table"""

    __tablename__ = 'run'
    runid = Column(Integer, primary_key=True, autoincrement=True)
    description = Column(Unicode)
    ms_filename = Column(Unicode)
    abs_peak_cutoff = Column(Float)
    max_ms_level = Column(Integer)
    ionisation_mode = Column(Integer)
    skip_fragmentation = Column(Boolean)
    max_broken_bonds = Column(Integer)
    max_water_losses = Column(Integer)
    ms_intensity_cutoff = Column(Float)
    msms_intensity_cutoff = Column(Float)
    mz_precision = Column(Float)
    mz_precision_abs = Column(Float)
    precursor_mz_precision = Column(Float)
    use_all_peaks = Column(Boolean)
