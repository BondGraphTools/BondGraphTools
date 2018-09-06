import pytest
from BondGraphTools import new
from BondGraphTools.exceptions import InvalidComponentException


def test_uri_tree():

    model = new(name="Model")

    internal_model = new(Name="Model")
    r_i = new('R')
    internal_model.add(r_i)

    # should raise
    with pytest.raises(InvalidComponentException) as ex:
        model.add(model)

    model.add(internal_model)

    assert internal_model.uri == "/internal_model"

    c_o = new("C", name="C_1")
    model.add(c_o)

    assert c_o.uri == "/C_1"


    assert internal_model in model.components
    assert r_i in internal_model
    assert r_i not in model

    assert r_i.uri == "/internal_model/R1"
    assert r_i.name == "R1"

    c_i = new("C", name="C_1")
    internal_model.add(
        c_i
    )

    assert c_i.uri == "/internal_model/C_1"

    model.remove(internal_model)
    assert internal_model not in model.components
    assert c_i.uri == "/internal_model/C_1"

