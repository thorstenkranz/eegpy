import os
from nose.tools import assert_true, assert_false, assert_equal
from eegpy.temp import UnlinkingTempfile

class TestTemp:
    def test_tempfile_is_unlinked(self):
        with UnlinkingTempfile() as fn:
            with open(fn,"w") as fh:
                fh.write("Test")
            assert_true(os.path.exists(fn))
        assert_false(os.path.exists(fn))
    
    def test_file_is_never_created(self):
        with UnlinkingTempfile() as fn:
            assert_false(os.path.exists(fn))
        assert_false(os.path.exists(fn))
        
    def test_suffix_is_used(self):
        with UnlinkingTempfile(suffix=".xyz") as fn:
            assert_true(fn.endswith(".xyz"))
        
    def test_prefix_is_used(self):
        with UnlinkingTempfile(prefix="xyz") as fn:
            file_part = os.path.split(fn)[1]
            assert_true(file_part.startswith("xyz"))
        
    def test_dir_is_used(self):
        dir_ = "/path/to/somewhere"
        with UnlinkingTempfile(dir=dir_) as fn:
            assert_equal(os.path.normpath(dir_),
                         os.path.normpath(os.path.split(fn)[0]))
            