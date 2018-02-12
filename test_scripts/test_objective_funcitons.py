import numpy as np
import pytest
import BondGraphTools.core.layout_manager as lm

nodes_0 = ((0, 0), (0, 0))
edges_0 = ((0, 1),)

nodes_1 = ((0, 0), (0, 0), (0, 0))
edges_1 = ((0, 1), (1, 2), (2, 0))

nodes_2 = ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))
edges_2 = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5))


class TestGraph0():

    def test_objective_zero(self):

        adj, gr_dist, (x, y) = lm._initial_conditions(nodes_0, edges_0)
        X = np.zeros_like(x)
        Y = np.zeros_like(y)

        PHI, (dx, dy) = lm._objective(
            X, Y, adj, gr_dist, w_0=1, w_1=0, eps_m=0.01
        )

        assert PHI == 100
        for xp in dx:
            assert xp == 0
        for yp in dy:
            assert yp == 0

    # def test_objective_one(self):
    #     adj, gr_dist, (x, y) = lm._initial_conditions(nodes_0, edges_0)
    #     x[0] = 0
    #     x[1] = 1
    #     Y = np.zeros_like(y)
    #
    #     PHI, (dx, dy) = lm._objective(
    #         X, Y, adj, gr_dist, w_0=1, w_1=0, eps_m=0.01
    #     )


    @pytest.mark.xfail()
    def test_simulated_annealing(self):

        pos = lm.simulated_annealing(nodes_0, edges_0)

        dist = lm._distance_matrix(pos)

        assert dist.max() < 2
        n, m = dist.shape
        for i in range(n):
            for j in range(m):
                if i != j:
                    assert dist[i, j] != 0, str((i,j))



