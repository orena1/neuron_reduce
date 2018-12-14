from .subtree_reductor_func import subtree_reductor
import os


def run_tests():
    
    pp = os.path.realpath(__file__)
    pp = os.path.dirname(pp)
    print(pp)
    execfile(pp + '/../tests/run_all_tests.py')
