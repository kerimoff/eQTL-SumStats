import os
import pytest
import sumstats.snp.loader as loader
from tests.snp.test_constants import *


class TestUnitLoader(object):
    h5file = ".testfile.h5"
    f = None

    def setup_method(self, method):
        # open h5 file in read/write mode
        self.f = h5py.File(self.h5file, mode='a')

    def teardown_method(self, method):
        os.remove(self.h5file)

    def test_open_with_empty_array(self):
        other_array = []

        dict = {"snp": snpsarray, "pval": pvalsarray, "chr": chrarray, "or": orarray, "bp": bparray,
                "effect": effectarray, "other": other_array, 'freq': frequencyarray}

        with pytest.raises(ValueError):
            loader.Loader(None, self.h5file, "PM001", dict)

    def test_open_with_None_array(self):

        other_array = None

        dict = {"snp": snpsarray, "pval": pvalsarray, "chr": chrarray, "or": orarray, "bp": bparray,
                "effect": effectarray, "other": other_array}

        with pytest.raises(ValueError):
            loader.Loader(None, self.h5file, "PM001", dict)

    def test_create_dataset(self):
        random_group = self.f.create_group("random_group")
        data = 'string1'
        dset_name = STUDY_DSET
        loader.create_dataset(random_group, dset_name, data)
        dset = random_group.get(dset_name)
        assert dset is not None
        dataset = dset[:]
        assert len(dataset) == 1
        assert dataset[0] == data

        data = 1
        dset_name = BP_DSET
        loader.create_dataset(random_group, dset_name, data)
        dset = random_group.get(dset_name)
        assert dset is not None
        dataset = dset[:]
        assert len(dataset) == 1
        assert dataset[0] == data

        data = 0.2
        dset_name = OR_DSET
        loader.create_dataset(random_group, dset_name, data)
        dset = random_group.get(dset_name)
        assert dset is not None
        dataset = dset[:]
        assert len(dataset) == 1
        assert dataset[0] == data

        dset_name = "random name"
        with pytest.raises(KeyError):
            loader.create_dataset(random_group, dset_name, data)

    def test_expand_dataset(self):
        random_group = self.f.create_group("random group")

        data = 'string1'
        dset_name = STUDY_DSET
        loader.create_dataset(random_group, dset_name, data)
        data2 = 'random string4'
        loader.expand_dataset(random_group, dset_name, data2)

        dset = random_group.get(dset_name)
        assert dset is not None
        assert len(dset) == 2
        dataset = dset[:]
        assert dataset[0] == 'string1'
        assert dataset[1] == 'random string4'

        data = 1
        dset_name = CHR_DSET
        loader.create_dataset(random_group, dset_name, data)
        data2 = 2
        loader.expand_dataset(random_group, dset_name, data2)

        dset = random_group.get(dset_name)
        assert dset is not None
        assert len(dset) == 2
        dataset = dset[:]
        assert dataset[0] == data
        assert dataset[1] == data2

        data = 0.1
        dset_name = MANTISSA_DSET
        loader.create_dataset(random_group, dset_name, data)
        data2 = 0.2
        loader.expand_dataset(random_group, dset_name, data2)

        dset = random_group.get(dset_name)
        assert dset is not None
        assert len(dset) == 2
        dataset = dset[:]
        assert dataset[0] == data
        assert dataset[1] == data2

    def test_expand_not_existing_dataset(self):
        random_group = self.f.create_group("random group")

        data = 'string1'
        dset_name = STUDY_DSET
        loader.expand_dataset(random_group, dset_name, data)
        dset = random_group.get(dset_name)

        assert dset is not None
        dataset = dset[:]
        assert len(dataset) == 1
        assert dataset[0] == 'string1'

        data2 = 'str2'
        loader.expand_dataset(random_group, dset_name, data2)
        dset = random_group.get(dset_name)
        dataset = dset[:]
        assert len(dataset) == 2

    def test_already_loaded_snp_not_in_file(self):
        load = prepare_load_object_with_study(self.h5file, "PM001")
        load.load()

        snpsarray_new = ['rs1', 'rs1', 'rs1', 'rs1']
        dict = {"snp": snpsarray_new, "pval": pvalsarray, "chr": chrarray, "or": orarray, "bp": bparray,
                "effect": effectarray, "other": otherarray, 'freq': frequencyarray}

        load = loader.Loader(None, self.h5file, 'PM001', dict)
        assert not load.already_loaded()

    def test_already_loaded_snp_group_exists_but_with_no_data(self):
        snp = 'rs1'
        load = prepare_load_object_with_study(self.h5file, "PM001")
        load.file.create_group(snp)

        assert not load.already_loaded()

    def test_already_loaded_study_not_loaded(self):
        load = prepare_load_object_with_study(self.h5file, "PM003")
        load.load()

        load = prepare_load_object_with_study(self.h5file, "PM001")
        assert not load.already_loaded()

    def test_already_loaded_study_already_loaded(self):
        load = prepare_load_object_with_study(self.h5file, "PM003")
        load.load()

        load = prepare_load_object_with_study(self.h5file, "PM003")
        assert load.already_loaded()

    def test_study_not_loaded_correctly_snp_not_saved(self):
        # load whole dataset except last one
        dict = {"snp": snpsarray[:-1], "pval": pvalsarray[:-1], "chr": chrarray[:-1], "or": orarray[:-1], "bp": bparray[:-1],
                "effect": effectarray[:-1], "other": otherarray[:-1], 'freq': frequencyarray[:-1]}

        load = loader.Loader(None, self.h5file, 'PM003', dict)
        load.load()

        # give the object the correct dictionary
        load = prepare_load_object_with_study(self.h5file, "PM003")

        assert not load.load_completed()

    def test_study_not_loaded_correctly_info_missing(self):
        load = prepare_load_object_with_study(self.h5file, "PM003")
        load.load()

        load = prepare_load_object_with_study(self.h5file, "PM001")

        assert not load.load_completed()


def prepare_load_object_with_study(h5file, study):
    dict = {"snp": snpsarray, "pval": pvalsarray, "chr": chrarray, "or": orarray, "bp": bparray,
            "effect": effectarray, "other": otherarray, 'freq': frequencyarray}

    return loader.Loader(None, h5file, study, dict)