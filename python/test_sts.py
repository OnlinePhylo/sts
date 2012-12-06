from __future__ import division
import unittest
import sts
import os

class StateLogExpectationTestMixin(object):
    expected_average_survival = None
    expected_average_mrca_depth = None

    def setUp(self):
        super(StateLogExpectationTestMixin, self).setUp()
        with open(os.path.join('test-data', self.datafile)) as infile:
            self.state_log = sts.StateLog.of_json_file(infile)

    def test_expected_average_survival(self):
        self.assertAlmostEqual(self.state_log.average_survival(),
                               self.expected_average_survival)

    def test_expected_average_mrca_depth(self):
        self.assertAlmostEqual(self.state_log.average_mrca_depth(),
                               self.expected_average_mrca_depth)

class StateLogExpectationTest1(StateLogExpectationTestMixin, unittest.TestCase):
    datafile = '1.json'
    expected_average_survival = 1
    expected_average_mrca_depth = 11 / 6

class StateLogExpectationTest2(StateLogExpectationTestMixin, unittest.TestCase):
    datafile = '2.json'
    expected_average_survival = 5 / 9
    expected_average_mrca_depth = 5 / 3
