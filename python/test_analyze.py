import unittest
import analyze
import os

class StateLogExpectationTestMixin(object):
    expected_average_survival = None
    expected_average_mrca_depth = None

    def setUp(self):
        super(StateLogExpectationTestMixin, self).setUp()
        with open(os.path.join('test-data', self.datafile)) as infile:
            self.state_log = analyze.StateLog.of_json_file(infile)

    def test_expected_average_survival(self):
        self.assertAlmostEqual(self.state_log.average_survival(),
                               self.expected_average_survival)

    def test_expected_average_mrca_depth(self):
        self.assertAlmostEqual(self.state_log.average_mrca_depth(),
                               self.expected_average_mrca_depth)

class StateLogExpectationTest1(StateLogExpectationTestMixin, unittest.TestCase):
    datafile = '1.json'

class StateLogExpectationTest2(StateLogExpectationTestMixin, unittest.TestCase):
    datafile = '2.json'
