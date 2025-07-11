import numpy as np
import pytest
import logging
from phasebo.phase_field import PhaseField
from phasebo.phase_field_bo import PhaseFieldBO

@pytest.fixture
def dummy_data():
    # Minimal example data
    compositions = np.array([
    ['Li3 B1 S3',0.5],
    ['Li1 B1 O2',0.5],
    ['B1',-6.70565352255],
    ['Li1',-1.890104425],
    ['O1',-4.94850939875],
    ['S1',-4.12858174875]
    ], dtype=object)
    
    references = np.array([
    ['B1' ],
    ['Li1'],
    ['O1'],
    ['S1']
    ], dtype=object)

    ions = {'Li':1,'B':3, 'O':-2,'S':-2} 
    exceptions = []

    return compositions, references, ions, exceptions

def test_phase_field_init(dummy_data):
    logger = logging.getLogger('test_logger')
    compositions, references, ions, exceptions = dummy_data
    pf = PhaseField(compositions, references, ions, exceptions, logger=logger)
    assert isinstance(pf, PhaseField)
    assert len(pf.compositions) > 0
    assert hasattr(pf, "energies")

def test_phase_fieldbo_init(dummy_data):
    logger = logging.getLogger('test_logger')
    compositions, references, ions, exceptions = dummy_data
    pfbo = PhaseFieldBO(compositions, references, ions, exceptions=exceptions, logger=logger)
    assert isinstance(pfbo, PhaseFieldBO)

def test_get_candidates(dummy_data):
    logger = logging.getLogger('test_logger')
    compositions, references, ions, exceptions = dummy_data
    pf = PhaseField(compositions, references, ions, exceptions, logger=logger)
    candidates = pf.candidates
    # Make sure candidates do not include references
    for c in pf.references:
        assert c not in candidates

def test_get_random_seeds(dummy_data):
    logger = logging.getLogger('test_logger')
    compositions, references, ions, exceptions = dummy_data
    pf = PhaseField(compositions, references, ions, exceptions, logger=logger)
    fc, energies = pf.get_random_seeds(1)
    assert len(fc) == 1
    assert len(energies) == 1

def test_get_seeds_from_segments(dummy_data):
    logger = logging.getLogger('test_logger')
    compositions, references, ions, exceptions = dummy_data
    pf = PhaseField(compositions, references, ions, exceptions, logger=logger)
    fc, energies = pf.get_seeds_from_segments(disect=2)
    assert len(fc) == len(pf.seeds)
    assert len(energies) == len(pf.seeds)
