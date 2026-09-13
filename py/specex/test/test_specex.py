# Test of whether specex runs without crashing

from specex.specex import run_specex
import specex.io
import sys
import os
import tempfile
import shutil
import unittest

inputroot_pre = None
inputroot_exp = None

inputs_exist = False
if 'DESI_ROOT' in os.environ:
    desi_root = os.environ['DESI_ROOT']
    for specprod in ['fuji', 'iron', 'loa']:
        inputroot_pre = f'{desi_root}/spectro/redux/{specprod}/preproc/20201216/00068217'
        inputroot_exp = f'{desi_root}/spectro/redux/{specprod}/exposures/20201216/00068217'

        if os.path.exists(inputroot_pre) and os.path.exists(inputroot_exp):
            inputs_exist = True
            break

class TestSpecex(unittest.TestCase):

    def setUp(self):
        self.testdir = tempfile.mkdtemp()
        self.outpsf = self.testdir+'/psf-test.fits'

    def tearnDown(self):
        shutil.rmtree(self.testdir)

    @unittest.skipUnless(inputs_exist, 'Input image data not found')
    def test_specex(self):
        first_fiber='100'
        last_fiber='104'

        lamp_lines_file = specex.io.get_desi_linelist_file()
        self.assertTrue(os.path.exists(lamp_lines_file))

        #- Confirm that the output file doesn't already exist
        self.assertFalse(os.path.exists(self.outpsf))

        com=['desi_psf_fit', '-a', inputroot_pre+'/preproc-b1-00068217.fits',
             '--in-psf', inputroot_exp+'/shifted-input-psf-b1-00068217.fits',
             '--out-psf', self.outpsf,
             '--lamp-lines', lamp_lines_file,
             '--first-bundle', '4', '--last-bundle', '4', '--first-fiber', first_fiber, '--last-fiber', last_fiber,
             '--legendre-deg-wave', '3','--fit-continuum','--debug']

        run_specex(com)

        self.assertTrue(os.path.exists(self.outpsf))
