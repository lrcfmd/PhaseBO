import random

import numpy as np
import torch

from phasebo.utils.other import set_seeds


def test_set_seeds_repeats_every_generator():
    def draw():
        return random.random(), np.random.rand(), torch.rand(1).item()

    set_seeds(0)
    first = draw()
    set_seeds(0)
    assert draw() == first
