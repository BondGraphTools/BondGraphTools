import pytest


class TestPort:
    def test_in(self):
        from BondGraphTools.base import Port
        from BondGraphTools.actions import new

        c = new("C")
        one = new("1")

        port_c = Port(c, 0)

        assert c in port_c
        assert c is not port_c
        assert one is not port_c
        assert one not in port_c

    def test_cmp(self):
        from BondGraphTools.base import Port
        from BondGraphTools.actions import new

        c = new("C")
        one = new("1")

        port_c = Port(c, 0)
        port_one = Port(one, 0)

        assert port_c is not port_one
        assert port_c != port_one


class TestIntegrations:

    def test_get_exposed_port(self):
        from BondGraphTools.actions import new, expose
        model = new()
        port = new('SS')
        model.add(port)
        expose(port, label='port_1')

        port = model.get_port()
        assert port.index == 0
        assert port.name == "port_1"


class TestPortIntegration(object):
    def test_failstate(self):
        # see issue 110
        from BondGraphTools import new, expose, connect
        from BondGraphTools.exceptions import InvalidPortException
        model = new()
        port = new('SS')
        model.add(port)
        expose(port, label='port_1')

        test_model = new()
        component = new("0")  # could be anything
        test_model.add(component, model)

        target_port = (model, 'port_')     # typo in port name
        with pytest.raises(InvalidPortException):
            connect(component, target_port)

    def test_connection_feedback(self):

        from BondGraphTools import new, connect, add
        from BondGraphTools.exceptions import InvalidPortException
        model = new()

        c1 = new('C')
        c2 = new('C')
        c3 = new('C')

        add(model, c1, c2, c3)
        connect(c1, (c2, 0))
        with pytest.raises(InvalidPortException) as ex:
            connect((c2, 0), c3)

        assert "Port is already connected: " in str(ex)
