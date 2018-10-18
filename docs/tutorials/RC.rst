Tutorial: Driven Filter Circuit
===============================

Goal: build a simple passive filter using basic bond graph modelling
techniques.

Difficulty: Beginner

Requirement: BondGraphTools, jupyter.

How to follow: Enter each block of code in consecutive cells in a jupyter
notebook.

Part 1: Basic Use
-----------------

First, import `BondGraphTools` and create a new model with the name "RC"::

    import BondGraphTools as bgt
    model = bgt.new(name="RC")

Now create a new generalised linear resistor ('R' component), a generalised
linear capacitor ("C" component) with resistance and capacitance both set to 1,
and an equal effort ("0" junction) conservation law through which these
components share energy.::

    C = bgt.new("C", value=1)
    R = bgt.new("R", value=1)
    zero_law = bgt.new("0")

Add the newly created components to the model::

    bgt.add(model, R, C, zero_law)

Once the components are added to the model, connect the components and the law
together. (Note the first argument is the tail of the energy bond, the second
is the head)::

    bgt.connect(R, zero_law)
    bgt.connect(zero_law, C)

Draw the model to make sure everything is wired up::

    bgt.draw(model)

which produces a sketch of the network topology.

.. figure:: images/RC_1.svg
    :scale: 50 %
    :align: center


In order to know::

    timespan = [0, 5]
    x0 = [1]
    t, x = simulate(model, timespan=timespan, x0=x0)
    from matplotlib.pyplot import plot
    fig = plot(t,x)

.. figure:: images/RC_2.svg
    :scale: 50 %
    :align: center

Part 2: Control
---------------

Add Current source::

    Sf = new('Sf')
    add(model, Sf)
    connect(Sf, KCL)
    draw(model)

.. figure:: images/RC_3.svg
    :scale: 50 %
    :align: center

Equations:
    model.control_vars
    model.constitutive_relations

Simulating with constant effort::

    t, x = simulate(model, timespan=timespan, x0=x0, control_vars={'u_0':2})
    plot(t,x)

.. figure:: images/RC_4.svg
    :scale: 50 %
    :align: center


Simulating with time varying forcing::

    t, x = simulate(model, timespan=timespan, x0=x0, control_vars={'u_0':'sin(2*t)'})
    plot(t,x)


.. figure:: images/RC_5.svg
    :scale: 50 %
    :align: center

With a step function::

    step_fn = 't < 1 ? 1 : 0' # if t < 0 then 1 else 0
    t, x = simulate(model, timespan=timespan, x0=x0, control_vars={'u_0':step_fn})
    plot(t,x)

.. figure:: images/RC_6.svg
    :scale: 50 %
    :align: center

Inside a loop for different harmonics::

    fig = plt.figure()
    for i in range(4):
        func_text = "cos({i}t)".format(i=i)
        t_i, x_i = simulate(model, timespan=timespan, x0=x0, control_vars={'u_0':func_text})
        plot(t_i,x_i)

.. figure:: images/RC_7.svg
    :scale: 50 %
    :align: center
