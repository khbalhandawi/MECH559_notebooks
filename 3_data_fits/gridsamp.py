import numpy as np
from typing import Union, List

def gridsamp(bounds: np.ndarray, q: Union[int, np.ndarray, List[int]]) -> np.ndarray:
    """
    GRIDSAMP  n-dimensional grid over given range

    Parameters
    ----------
    bounds : np.ndarray
        2*n matrix with lower and upper limits
    q : np.ndarray
        n-vector, q(j) is the number of points
        in the j'th direction.
        If q is a scalar, then all q(j) = q

    Returns
    -------
    S : np.ndarray
        m*n array with points, m = prod(q)
    """

    mr, n = np.shape(bounds)
    dr = np.diff(bounds, axis=0)[0]  # difference across rows
    if mr != 2 or any([item < 0 for item in dr]):
        raise Exception('bounds must be an array with two rows and bounds(1,:) <= bounds(2,:)')

    if q.ndim > 1 or any([item <= 0 for item in q]):
        raise Exception('q must be a vector with non-negative elements')

    p = len(q)
    if p == 1:
        q = np.tile(q, (1, n))[0]
    elif p != n:
        raise Exception('length of q must be either 1 or %d' % n)

    # Check for degenerate intervals
    i = np.where(dr == 0)[0]
    if i.size > 0:
        q[i] = 0 * q[i]

    # Recursive computation
    if n > 1:
        a = gridsamp(bounds[:, 1::], q[1::])  # Recursive call
        [m, _] = np.shape(a)
        q = int(q[0])
        s = np.concatenate((np.zeros((m * q, 1)), np.tile(a, (q, 1))), axis=1)
        y = np.linspace(bounds[0, 0], bounds[1, 0], q)

        k = range(m)
        for i in range(q):
            aug = np.tile(y[i], (m, 1))
            aug = np.reshape(aug, s[k, 0].shape)

            s[k, 0] = aug
            k = [item + m for item in k]
    else:
        s = np.linspace(bounds[0, 0], bounds[1, 0], int(q[-1]))
        s = np.transpose([s])

    return s