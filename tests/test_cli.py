import os
import tempfile
import shutil

from nose.tools import nottest
import click.testing

import data
from codon import cli

class _Base(object):
    def setup(self):
        self.tempdir = os.path.join(tempfile.gettempdir(), 'codon')
        try:
            os.mkdir(self.tempdir)
        except OSError:
            pass

    def teardown(self):
        try:
            shutil.rmtree(self.tempdir)
        except OSError:
            pass

def compare_files(file1, file2):
    file1 = open(file1, 'r')
    file2 = open(file2, 'r')
    for line1, line2 in zip(file1, file2):
        print line1
        print line2
        assert line1.strip() == line2.strip(), line1+line2

class TestFit(_Base):
    def testFit(self):
        ''' fit should fit a model '''
        datadir = data.get_data_dir()
        aln = os.path.join(datadir, 'aln.fasta')
        tree = os.path.join(datadir, 'tree.nwk')
        correct_result = os.path.join(datadir, 'MG94GTR.txt')
        test_result = os.path.join(self.tempdir, 'MG94GTR.txt')
        args = ['fit', '--model', 'MG94GTR', aln, tree, test_result]
        runner = click.testing.CliRunner()
        result = runner.invoke(cli.main, args)
        compare_files(test_result, correct_result)

class TestBootstrap(_Base):
    @nottest
    def testFit(self):
        ''' fit should fit a model '''
        datadir = data.get_data_dir()
        aln = os.path.join(datadir, 'aln.fasta')
        tree = os.path.join(datadir, 'tree.nwk')
        fit_result = os.path.join(datadir, 'MG94GTR.json')
        correct_result = os.path.join(datadir, 'MG94GTR.bootstrap')
        test_result = os.path.join(self.tempdir, 'MG94GTR.bootstrap')
        args = ['bootstrap', '--num_bootstraps', '1', fit_result, test_result]
        runner = click.testing.CliRunner()
        result = runner.invoke(cli.main, args)
        compare_files(test_result, correct_result)

class TestOmega(_Base):
    def testFit(self):
        ''' omega should fit a model with omega constraints'''
        datadir = data.get_data_dir()
        aln = os.path.join(datadir, 'aln.fasta')
        tree = os.path.join(datadir, 'tree.nwk')
        correct_result = os.path.join(datadir, 'Y98.txt')
        test_result = os.path.join(self.tempdir, 'Y98.txt')
        args = ['omega', '--model', 'Y98', '--outgroup', 'Mouse', aln, tree, 
                test_result]
        runner = click.testing.CliRunner()
        result = runner.invoke(cli.main, args)
        compare_files(test_result, correct_result)

class TestClock(_Base):
    def testClock(self):
        ''' clock should fit a clock type model '''
        datadir = data.get_data_dir()
        aln = os.path.join(datadir, 'aln.fasta')
        tree = os.path.join(datadir, 'tree.nwk')
        correct_result = os.path.join(datadir, 'MG94GTRClock.txt')
        test_result = os.path.join(self.tempdir, 'MG94GTRClock.txt')
        args = ['clock', '--model', 'MG94GTR', aln, tree, 'Mouse', test_result]
        runner = click.testing.CliRunner()
        result = runner.invoke(cli.main, args)
        compare_files(test_result, correct_result)

#class TestRooted(_Base):
#    pass
