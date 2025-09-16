#!/usr/bin/env python

### Tests

import unittest

from simulation import simulate_outbreak
from poisson_process_generator import PoissonProcessGenerator

class TestSimulateOutbreak(unittest.TestCase):

    def test_extinct_1(self):
        poisson_presets = [0]
        poisson_process_generator = PoissonProcessGenerator(shape=None, scale=None, presets=poisson_presets)
        
        target = 20

        actual = simulate_outbreak(poisson_process_generator, target)
        expected_extinct = True
        expected_history = [[1],
                            [0]]
        expected_population = 1
        
        self.assertEqual(actual["history"], expected_history)
        self.assertEqual(actual["extinct"], expected_extinct)
        self.assertEqual(actual["population"], expected_population)

    def test_extinct_gen3(self):
        poisson_presets = [9,2,0,4,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0]
        poisson_process_generator = PoissonProcessGenerator(shape=None, scale=None, presets=poisson_presets)
        
        target = 20

        actual = simulate_outbreak(poisson_process_generator, target)
        expected_extinct = True
        expected_history = [[1],
                            [9],
                            [2,0,4,0,0,0,0,2,0],
                            [0,0,0,0,0,0,0,0]]
        expected_population = 18
        
        self.assertEqual(actual["history"], expected_history)
        self.assertEqual(actual["extinct"], expected_extinct)
        self.assertEqual(actual["population"], expected_population)

    def test_passes(self):
        poisson_presets = [9,2,0,4,0,0,0,0,2,0,1,0,1,1,3,5,0,6,3,0,5,0,0,1,0,2,2,0,5,2,1,0,0,1,0,0,0,2,2,1,0,0,3,1,0,0,3,1,0,0,1,2,0,2,0,3,0,0,0,1,9,0,2,1,0,3,0,0,2,1,1,1,5,1,0,0,8,1,0,4,1,1,1,0,0,5,4,0,0,0,0,1,1,6,0,3,1,0,3,4,0,0,1,1,6,6,0,9,0,1,10,2,1,4,1,0,0,0,0,1,2,7,1,1,1,0,0,0,8,0,0,2,0,0,0,2,0,0,1,1,0,3,0,1,1,0,2,0,0,2,3,5,2,1,2,0,3,4,1,0,1,1,0,8,2,1,0,3,0,5,0,1,0,4,0,0,3,0,0,0,0,2,0,1,0,0,1,0,2,1,0,1,1,1,2,0,1,1,0,0,0,2,3,4,0,0,2,0,2,0,1,1,0,13,3,1,4,1,1,0,0,3,0,0,2,0,0,2,1,1,0,2,2,0,1,0,0,0,7,0,5,1,2,0,0,3,1,0,1,5,2,3,2,0,1,5,1,1,4,3,1,1,0,2,0,3,1,0,1,0,3,1,0,4,2,0,0,0,1,1,0,0,0,0,0,3,6,2,2,4,1,0,2,0,1,2,0,0,0,1,1,0,8,1,1,3,2,1,1,0,0,0,1,2,0,4,1,0,0,1,2,3,1,0,0,0,5,2,1,0,0,6,1,0,1,2,0,2,0,0,2,0,0,1,4,2,1,2,4,1,3,2,1,0,0,2,1,2,0,0,4,3,1,3,1,0,2,4,1,4,5,0,1,2,0,3,2,0,6,0,5,1,1,0,1,0,0,0,0,1,3,2,1,1,0,0,0,3,2,3,0,2,0,1,1,1,1,2,1,2,0,0,1,7,3,1,3,0,1,1,2,0,5,0,0,0,1,0,2,3,2,0,4,0,0,1,1,0,0,0,1,1,1,1,1,2,1,8,0,0,0,3,4,2,3,5,0,2,1,0,1,1,0,0,0,0,5,1,0,1,4,2,1,1,1,0,1,0,0,1,0,5,0,0,0,2,0,3,0,1,3,2,0,6,3,0,1,5,2,4,2,2,0,2,2,1,1,0,1,2,2,7,9,0,0,3,4,2,4,1,2,0,1,1,2,1,1,2,4,0,3,1,0,0,1,0,2,0,0,0,3,0,0,0,1,2,2,2,0,2,0,0,1,0,0,1,4,2,4,1,5,0,1,3,0,0,1,0,0,0,0,0,1,4,1,1,0,2,0,0,0,0,2,0,4,0,0,6,3,0,0,2,1,1,3,0,4,0,0,2,0,2,1,1,1,0,5,4,1,1,9,2,3,2,0,0,0,0,1,1,5,3,1,3,3,1,0,1,0,2,2,0,0,0,4,1,0,0,1,5,1,0,2,3,0,1,2,0,3,1,0,0,0,0,3,0,8,1,0,1,1,1,2,5,1,0,1,1,0,3,1,0,2,3,3,0,8,1,1,2,1,0,5,0,6,1,1,0,1,2,1,0,1,0,2,1,0,3,0,1,4,1,0,5,2,12,1,2,0,1,2,1,2,0,1,1,1,8,3,0,0,1,0,0,2,5,0,4,3,0,0,2,1,4,1,0,4,1,1,1,0,3,0,0,0,0,0,1,0,2,1,2,2,1,2,2,2,1,0,3,1,0,0,0,0,3,3,2,0,1,0,3,4,3,1,2,0,2,2,0,0,1,0,1,0,1,2,0,0,11,0,2,1,0,0,4,10,1,2,0,1,1,6,1,1,0,0,0,1,0,3,1,1,0,0,1,0,0,2,1,0,2,0,0,1,5,2,1,7,5,1,2,1,1,2,0,1,0,0,1,1,1,3,1,3,3,1,0,4,1,2,0,0,0,0,1,1,2,0,5,0,2,2,2,0,8,8,0,3,2,0,0,2,8,7,2,1,1,0,0,1,0,0,0,2,0,0,0,5,1,3,1,1,0,0,1,0,2,0,2,4,0,2,0,4,0,2,1,9,0,0,8,0,0,0,0,2,2,0,4,0,0,2,1,0,2,5,0,0,0,14,4,1,0,1,1,0,1,0,0,2,1,0,2,0,1,0,2,1,0,1,4,2,1,2,1,1,1,2,0,0,2,1,4,1,1,0,0,1,2,5,1,1,1,2,0,0,0,0,2,6,2,0,5,0,1,0,3,1,3,1,0,0,0,0,1,4,0,3,0]        
        poisson_process_generator = PoissonProcessGenerator(shape=None, scale=None, presets=poisson_presets)
        
        target = 20
        
        actual = simulate_outbreak(poisson_process_generator, target)
        expected_extinct = False
        expected_history = [[1],
                            [9],
                            [2,0,4,0,0,0,0,2,0],
                            [1,0,1,1,3,5,0,6],
                            [3,0,5,0,0,1,0,2,2,0,5,2,1,0,0,1,0]]        
        expected_population = 57
                
        self.assertEqual(actual["history"], expected_history)
        self.assertEqual(actual["extinct"], expected_extinct)
        self.assertEqual(actual["population"], expected_population)

if __name__ == '__main__':
    unittest.main()