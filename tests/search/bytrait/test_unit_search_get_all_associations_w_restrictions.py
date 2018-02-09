import os
import shutil

import sumstats.search as search
import tests.search.search_test_constants as search_arrays
import sumstats.trait.loader as loader
from tests.prep_tests import *
from sumstats.trait.constants import *
from tests.search.test_utils import *
from sumstats.utils.interval import FloatInterval


class TestLoader(object):

    output_location = './output/bytrait/'
    file = None
    start = 0
    size = 20

    def setup_method(self, method):
        # output is always stored under a directory called: 'output'
        os.makedirs('./output/bytrait')

        # loaded s1/t1 -> 50 associations
        # loaded s2/t1 -> 50 associations
        # loaded s3/t2 -> 50 associations
        # loaded s4/t2 -> 50 associations
        # total associations loaded : 200

        search_arrays.chrarray = [1 for _ in range(50)]
        search_arrays.pvalsarray = ["0.00001" for _ in range(50)]
        h5file = self.output_location + 'file_t1.h5'
        load = prepare_load_object_with_study_and_trait(h5file=h5file, study='s1', trait='t1', loader=loader, test_arrays=search_arrays)
        load.load()

        search_arrays.chrarray = [2 for _ in range(50)]
        search_arrays.pvalsarray = ["0.0001" for _ in range(50)]
        load = prepare_load_object_with_study_and_trait(h5file=h5file, study='s2', trait='t1', loader=loader, test_arrays=search_arrays)
        load.load()

        h5file = self.output_location + 'file_t2.h5'
        search_arrays.chrarray = [1 for _ in range(50)]
        search_arrays.pvalsarray = ["0.05" for _ in range(50)]
        load = prepare_load_object_with_study_and_trait(h5file=h5file, study='s3', trait='t2', loader=loader, test_arrays=search_arrays)
        load.load()

        search_arrays.chrarray = [2 for _ in range(50)]
        search_arrays.pvalsarray = ["0.1" for _ in range(50)]
        load = prepare_load_object_with_study_and_trait(h5file=h5file, study='s4', trait='t2', loader=loader, test_arrays=search_arrays)
        load.load()

        # initialize searcher with local path
        self.searcher = search.Search(path="./output")

    def teardown_method(self, method):
        shutil.rmtree('./output')

    def test_get_all_loop_through_size_20_restrict_pval_to_s1(self):
        start = 0
        size = 200
        pval_interval = FloatInterval().set_tuple(0.00001, 0.00001)
        datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)

        assert_studies_from_list(datasets, ['s1'])
        assert_datasets_have_size(datasets, TO_QUERY_DSETS, 50)

    def test_get_40_60_restrict_pval_to_s2(self):
        start = 40
        size = 20
        # p-value range of s2
        pval_interval = FloatInterval().set_tuple(0.00005, 0.0001)
        datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)

        assert_studies_from_list(datasets, ['s2'])
        assert_datasets_have_size(datasets, TO_QUERY_DSETS, 20)

    def test_get_all_restrict_pval_to_s4(self):
        start = 0
        size = 200
        pval_interval = FloatInterval().set_tuple(0.06, 0.1)
        datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)

        assert_studies_from_list(datasets, ['s4'])
        assert_datasets_have_size(datasets, TO_QUERY_DSETS, 50)

    def test_get_10_between_studies_s3_and_s4(self):
        start = 145
        size = 10
        pval_interval = FloatInterval().set_tuple(0.04, 0.3)

        datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)
        assert_datasets_have_size(datasets, TO_QUERY_DSETS, 10)
        assert_studies_from_list(datasets, ['s3', 's4'])
        assert_number_of_times_study_is_in_datasets(datasets, 's3', 5)
        assert_number_of_times_study_is_in_datasets(datasets, 's4', 5)

    def test_get_10_between_studies_s3_and_s4_restrict_to_s3(self):
        start = 145
        size = 10
        pval_interval = FloatInterval().set_tuple(0.04, 0.06)

        datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)
        print(datasets)
        assert_datasets_have_size(datasets, TO_QUERY_DSETS, 5)
        assert_studies_from_list(datasets, ['s3'])
        assert_number_of_times_study_is_in_datasets(datasets, 's3', 5)

    def test_loop_through_w_restrinction_and_always_get_size_20_results(self):
        start = 0
        size = 20

        looped_through = 1

        # s2 and s3 p-value limits
        pval_interval = FloatInterval().set_tuple(0.00002, 0.06)

        while True:
            datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)
            if len(datasets[REFERENCE_DSET]) <= 0:
                break
            if looped_through <= 2:
                assert_studies_from_list(datasets, ['s2'])
                assert_datasets_have_size(datasets, TO_QUERY_DSETS, 20)
            elif looped_through == 3:
                assert_studies_from_list(datasets, ['s2', 's3'])
                assert_number_of_times_study_is_in_datasets(datasets, 's2', 10)
                assert_number_of_times_study_is_in_datasets(datasets, 's3', 10)
            else:
                assert_studies_from_list(datasets, ['s3'])
                assert_datasets_have_size(datasets, TO_QUERY_DSETS, 20)
            looped_through += 1
            start = start + next_index

    def test_loop_through_w_restrinction_and_always_get_size_45_results(self):
        start = 0
        size = 45

        looped_through = 1

        # s2 and s3 p-value limits
        pval_interval = FloatInterval().set_tuple(0.00002, 0.06)

        while True:
            datasets, next_index = self.searcher.search_all_assocs(start=start, size=size, pval_interval=pval_interval)
            if len(datasets[REFERENCE_DSET]) <= 0:
                break
            if looped_through <= 1:
                assert_studies_from_list(datasets, ['s2'])
                assert_datasets_have_size(datasets, TO_QUERY_DSETS, 45)
            elif looped_through == 2:
                assert_studies_from_list(datasets, ['s2', 's3'])
                assert_number_of_times_study_is_in_datasets(datasets, 's2', 5)
                assert_number_of_times_study_is_in_datasets(datasets, 's3', 40)
            else:
                assert_studies_from_list(datasets, ['s3'])
                assert_datasets_have_size(datasets, TO_QUERY_DSETS, 10)
            looped_through += 1
            start = start + next_index

    def test_90_110_contains_s2_and_s3(self):
        start = 90
        size = 20
        datasets, next_index = self.searcher.search_all_assocs(start=start, size=size)
        print(datasets[STUDY_DSET])
        assert_studies_from_list(datasets, ['s2', 's3'])
        assert_number_of_times_study_is_in_datasets(datasets, 's2', 10)
        assert_number_of_times_study_is_in_datasets(datasets, 's3', 10)