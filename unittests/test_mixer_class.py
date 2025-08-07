import unittest
from Mixer import LXMixer

class TestLXMixer(unittest.TestCase):
    def test_compute_family_of_graphs(self):
        correct_family_of_valid_graphs = {
            0b0010: [(0, 1)],
            0b0101: [(1, 2)],
            0b0111: [(0, 2), (3, 4)],
            0b1000: [(1, 3)],
            0b1010: [(0, 3), (2, 4)],
            0b1101: [(0, 4), (2, 3)],
            0b1111: [(1, 4)]
        }
        
        B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
        lxmixer = LXMixer(B, 4)
        lxmixer.compute_family_of_valid_graphs()
        
        # Sort edges for comparison
        for key in correct_family_of_valid_graphs:
            correct_family_of_valid_graphs[key] = sorted(correct_family_of_valid_graphs[key])
        
        for key in lxmixer.family_of_valid_graphs:
            lxmixer.family_of_valid_graphs[key] = sorted(lxmixer.family_of_valid_graphs[key])
        
        self.assertEqual(lxmixer.family_of_valid_graphs, correct_family_of_valid_graphs)
    
    def test_compute_all_orbits(self):
        correct_orbits = {
            (0, 2, 3, 4) : [0b1101, 0b1010],
            (0, 1) : [0b0010],
            (1, 2) : [0b0101],
            (1, 3) : [0b1000],
            (1, 4) : [0b1111]
        }
        
        B = [0b1110, 0b1100, 0b1001, 0b0100, 0b0011]
        lxmixer = LXMixer(B, 4)
        lxmixer.family_of_valid_graphs = {
            0b0010: [(0, 1)],
            0b0101: [(1, 2)],
            0b0111: [(0, 2), (3, 4)],
            0b1000: [(1, 3)],
            0b1010: [(0, 3), (2, 4)],
            0b1101: [(0, 4), (2, 3)],
            0b1111: [(1, 4)]
        }
        lxmixer.compute_all_orbits()
        
        orbits = {}
        for nodes, orbit in lxmixer.orbits.items():
            orbits[nodes] = orbit.Xs
        
        self.assertEqual(orbits, correct_orbits)

class TestStabilizer(unittest.TestCase):
    def setUp(self):
        self.test_cases = [
            {
                "B": [0b1011, 0b1100, 0b0111, 0b0000, 0b1110, 0b1001, 0b0010, 0b0101],
                "n": 4,
                "orbit_dictionary": {(0,1,2,3,4,5,6,7): Orbit(Xs={11, 14, 7})},
                "expected_orbits": {(0,1,2,3,4,5,6,7) : [11, 14, 7]},
                "expected_mgs": {(0,1,2,3,4,5,6,7):[(1, 13)]},
                "expected_projectors": {(0,1,2,3,4,5,6,7): [(1, 0), (1, 13)]}
            },
            {
                "B": [0b11000, 0b00100, 0b01101, 0b10001],
                "n": 5,
                "orbit_dictionary": {(0,1,2,3): Orbit(Xs={9, 28})},
                "expected_orbits": {(0,1,2,3): [9, 28]},
                "expected_mgs": {(0,1,2,3) : [(1, 2), (-1, 20), (1, 25)]},
                "expected_projectors": {(0,1,2,3) : [(1, 0), (1, 25), (-1, 20), (-1, 13), (1, 2), (1, 27), (-1, 22), (-1, 15)]}
            }
            # }, #TODO the one below is not working....
            # {
            #     "B": [0b1110, 0b1100, 0b1001, 0b0100, 0b0011],
            #     "n": 4,
            #     "orbit_dictionary": {(0,2,3,4): Orbit(Xs={10, 13, 7}), (0,1): Orbit(Xs={2}), (1,2): Orbit(Xs={5}), (1,3): Orbit(Xs={8}), (1,4): Orbit(Xs={15})},
            #     "expected_orbits": {(0,1,2,3,4,5,6,7) : [11, 14, 7]},
            #     "expected_mgs": {(0,1,2,3,4,5,6,7):[(1, 13)]},
            #     "expected_projectors": {(0,2,3,4): [(1,0), (-1, 14), (-1, 5), (1, 11)], (0,1): [(1, 0), (1, 1), (-1, 4), (-1, 5), (-1, 8), (-1, 9), (1, 12), (1, 13)], (1,2): [(1, 0), (-1, 5), (1, 2), (-1, 7), (-1, 8), (1, 13), (-1, 10), (1, 15)], (1,3): [(1, 0), (1, 1), (1, 2), (1, 3), (-1, 4), (-1, 5), (-1, 6), (-1, 7)], (1,4): [(1, 0), (-1, 9), (-1, 10), (1, 3), (1, 12), (-1, 5), (-1, 6), (1, 15)]}
            # }
        ]

    def test_compute_minimal_generating_sets(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(B=case["B"], n=case["n"], orbit_dictionary=case["orbit_dictionary"])
                stab.compute_minimal_generating_sets()
                for orbit_key, orbit in stab.orbit_dictionary.items():
                    with self.subTest(orbit_key=orbit_key):
                        self.assertEqual(orbit.Zs, case["expected_mgs"][orbit_key])
    
    def test_compute_projector_stabilizers(self):
        for case in self.test_cases:
            with self.subTest(case=case):
                case = copy.deepcopy(case)
                stab = Stabilizer(B=case["B"], n=case["n"], orbit_dictionary=case["orbit_dictionary"])
                for orbit_key, orbit in stab.orbit_dictionary.items():
                    with self.subTest(orbit_key=orbit_key):
                        orbit.Zs = case["expected_mgs"][orbit_key]
                
                stab.compute_projector_stabilizers()
                for orbit_key, orbit in stab.orbit_dictionary.items():
                    with self.subTest(orbit_key=orbit_key):
                        self.assertEqual(orbit.Zs, case["expected_projectors"][orbit_key])

    
if __name__ == '__main__':
    unittest.main()