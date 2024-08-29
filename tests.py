from main import (add_molecule, substructure_search, delete_molecule_by_id, retrieve_molecule, retrieve_molecule_by_id,
                  update_molecule)
from fastapi import HTTPException
import pytest


@pytest.mark.parametrize("molecule_id, mol, expected", [(5, "CC(=O)", {"molecule_id": 5, "molecule": "CC(=O)"}),
                                                        (6, "c1ccccc1", {"molecule_id": 6, "molecule": "c1ccccc1"})])
def test_add_molecule(molecule_id, mol, expected):
    result = add_molecule(molecule_id, mol)
    assert result == expected

    with pytest.raises(HTTPException):
        add_molecule(1, "CCO")


@pytest.mark.parametrize("molecule_id, mol, limit, expected", [(5, "CC(=O)", 2, [{"molecule_id": 3, "molecule": "CC(=O)O"},
                                                                       {"molecule_id": 4,
                                                                        "molecule": "CC(=O)Oc1ccccc1C(=O)O"}]),
                                                        (6, "c1ccccc1", 1, [{"molecule_id": 2, "molecule": "Cc1ccccc1"}])])
def test_substructure_search(molecule_id, mol, limit, expected):
    try:
        add_molecule(molecule_id, mol)
    except HTTPException:
        print("The molecule is already in the molecule database")
    finally:
        result = substructure_search(molecule_id, mol, limit)
        assert result == expected

    with pytest.raises(HTTPException):
        substructure_search(10, "O", 2)


def test_substructure_search2():
    try:
        substructure_search(10, 'CC(=O)', 1)
    except HTTPException:
        print("Molecule is not added")


@pytest.mark.xfail
def test_substructure_search3():
    assert substructure_search(5, 'CC(=O)', 1) == [{"molecule_id": 3, "molecule": "CC(=O)O"}]


@pytest.mark.parametrize("molecule_id, expected", [(5, {"molecule_id": 5, "molecule": "CC(=O)"}),
                                                   (6, {"molecule_id": 6, "molecule": "c1ccccc1"})])
def test_delete_molecule_by_id(molecule_id, expected):
    result = delete_molecule_by_id(molecule_id)
    assert result == expected

    with pytest.raises(HTTPException):
        delete_molecule_by_id(10)


def test_retrieve():
    assert retrieve_molecule() == [
        {"molecule_id": 1, "molecule": "CCO"},
        {"molecule_id": 2, "molecule": "Cc1ccccc1"},
        {"molecule_id": 3, "molecule": "CC(=O)O"},
        {"molecule_id": 4, "molecule": "CC(=O)Oc1ccccc1C(=O)O"}
    ]


@pytest.mark.parametrize("molecule_id, expected", [(1, {"molecule_id": 1, "molecule": "CCO"}),
                                                   (2, {"molecule_id": 2, "molecule": "Cc1ccccc1"}),
                                                   (3, {"molecule_id": 3, "molecule": "CC(=O)O"}),
                                                   (4, {"molecule_id": 4, "molecule": "CC(=O)Oc1ccccc1C(=O)O"})])
def test_retrieve_molecule_by_id(molecule_id, expected):
    result = retrieve_molecule_by_id(molecule_id)
    assert result == expected


def test_update_molecule():
    assert update_molecule(1, {"molecule_id": 1, "molecule": "O"}) == [
        {"molecule_id": 1, "molecule": "O"},
        {"molecule_id": 2, "molecule": "Cc1ccccc1"},
        {"molecule_id": 3, "molecule": "CC(=O)O"},
        {"molecule_id": 4, "molecule": "CC(=O)Oc1ccccc1C(=O)O"},
    ]
