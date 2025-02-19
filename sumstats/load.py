import sys
import subprocess
import glob
import os
import argparse
import pandas as pd
import tables as tb
from sumstats.common_constants import *
from sumstats.utils.properties_handler import properties
from sumstats.utils import properties_handler
from sumstats.utils import filesystem_utils as fsutils
import sumstats.utils.sqlite_client as sq
import sumstats.reindex as ri


class Loader():
    def __init__(self, tsv=None, var_file=None, qtl_group=None, quant_method=None,
                 study_dir=None, study=None, trait=None, hdf_path=None,
                 chromosome=None, expr_file=None, sqldb=None, loader=None, trait_dir=None, tissue=None,
                 tissue_ont=None, treatment=None, treatment_ont=None, condition_label=None):
        self.tsv = tsv
        self.var_file = var_file
        self.qtl_group = qtl_group
        self.quant_method = quant_method
        self.tissue = tissue
        self.study = study
        self.trait_file = trait
        self.chromosome = chromosome
        self.hdf_path = hdf_path
        self.study_dir = study_dir
        self.trait_dir = trait_dir
        self.expr_file = expr_file
        self.max_string = 255
        self.tissue_ont = tissue_ont
        self.treatment = treatment
        self.treatment_ont = treatment_ont
        self.condition_label = condition_label

        self.filename = None
        if self.tsv:
            self.filename = os.path.splitext(os.path.basename(self.tsv))[0]
        self.traits = trait

        self.sqldb = sqldb

    def load_bystudy(self):
        print(self.tsv)
        identifier = self.study + "+" + self.qtl_group + "+" + self.quant_method
        #group = "/{study}".format(study=self.study.replace('-','_'))
        #hdf_store = fsutils.create_h5file_path(path=self.hdf_path, file_name=identifier, dir_name=self.study_dir + "/" + self.chromosome)
        hdf_store = os.path.join(self.chromosome, "file_" + identifier + ".h5")
        self.write_csv_to_hdf(hdf_store, identifier)


    def load_traits(self):
        hdf_store = fsutils.create_h5file_path(path=self.hdf_path, file_name="phen_meta", dir_name=self.trait_dir)
        dftrait = pd.read_csv(self.trait_file, sep="\t")
        self.write_traitfile_to_hdf(hdf_store, os.path.basename(self.trait_file))


    def write_traitfile_to_hdf(self, hdf, group):
        """Read in the trait file"""
        dftrait = pd.read_csv(self.trait_file, sep="\t")
        #mapping = {'X': 23, 'x': 23, 'Y': 24, 'y': 24}
        dftrait[CHR_DSET].replace({'X': 23, 'x': 23, 'Y': 24, 'y': 24, 'MT': 25}, inplace=True)
        dftrait = dftrait[pd.to_numeric(dftrait[CHR_DSET], errors='coerce').notnull()]
        dftrait[CHR_DSET] = dftrait[CHR_DSET].astype(int)

        with pd.HDFStore(hdf) as store:
            """store in hdf5 as below"""
            dftrait.to_hdf(store, group,
                        complib='blosc',
                        complevel=9,
                        format='table',
                        mode='a',
                        data_columns=['phenotype_id', 'gene_id']
                        )


    def write_csv_to_hdf(self, hdf, group):
        """Read in the sumstats files in chunks"""

        dfss = pd.read_csv(self.tsv, sep="\t",
                           dtype=DSET_TYPES,
                           usecols=list(TO_LOAD_DSET_HEADERS_DEFAULT),
                           float_precision='high',
                           chunksize=1000000)


        with pd.HDFStore(hdf) as store:
            """store in hdf5 as below"""
            for chunk in dfss:
                chunk = chunk[(chunk[CHR_DSET] == self.chromosome)]
                if not chunk.empty:
                    print("store chunk in hdf5")
                    chunk = chunk[sorted(TO_LOAD_DSET_HEADERS_DEFAULT)]
                    for field in [EFFECT_DSET, OTHER_DSET]:
                        self.placeholder_if_allele_string_too_long(df=chunk, field=field)
                    self.placeholder_if_variant_id_too_long(df=chunk, field=SNP_DSET)
                    chunk.to_hdf(store, group,
                                complib='blosc',
                                complevel=9,
                                format='table',
                                mode='a',
                                append=True,
                                data_columns=list(TO_INDEX),
                                #expectedrows=num_rows,
                                min_itemsize={OTHER_DSET: self.max_string,
                                              EFFECT_DSET: self.max_string,
                                              PHEN_DSET: self.max_string,
                                              GENE_DSET: self.max_string,
                                              MTO_DSET: self.max_string,
                                              RSID_DSET: 24,
                                              CHR_DSET: 2,
                                              MUTATION_DSET: self.max_string,
                                              SNP_DSET: self.max_string},
                                index = False
                                )

                    """Store study specific metadata"""
                    store.get_storer(group).attrs.study_metadata = {'study': self.study,
                                                                    'qtl_group': self.qtl_group,
                                                                    'quant_method': self.quant_method}
            indexer = ri.H5Indexer(hdf)
            indexer.reindex_file()





    def placeholder_if_allele_string_too_long(self, df, field):
        mask = df[field].str.len() <= self.max_string
        df[field].where(mask, "LONG_STRING", inplace=True)


    def placeholder_if_variant_id_too_long(self, df, field):
        mask = df[field].str.len() <= self.max_string
        df[field].where(mask, df[field].str.replace(r'(.*_)(.*_)(.*_.*)$', r'\1\2LONG_STRING'), inplace=True)


    def load_study_info(self):
        """
            CREATE TABLE study_info (
            study blob not null,
            identifier blob not null,
            qtl_group blob not null,
            tissue blob not null,
            trait_file blob not null,
            tissue_ontology blob not null,
            treatment blob,
            treatment_ontology blob not null,
            quant_method blob,
            condition_label blob not null
            UNIQUE (identifier)
            );
        """
        sql = sq.sqlClient(self.sqldb)
        identifier = self.study + "+" + self.qtl_group + "+" + self.quant_method
        trait_file_id = os.path.basename(self.trait_file)
        data = [self.study, identifier, self.qtl_group, self.tissue, trait_file_id, self.tissue_ont, self.treatment, self.treatment_ont, self.quant_method , self.condition_label]
        sql.cur.execute("insert or ignore into study_info values (?,?,?,?,?,?,?,?,?, ?)", data)
        sql.cur.execute('COMMIT')
        #Schmiedel_2018|Schmiedel_2018+Treg_memory+exon|Treg_memory|Treg|Ensembl_96_phenotype_metadata.tsv.gz|CL_0002678|Treg_memory|None|exon|naive


def update_metadate_from_tsv(metadata, ss_file):
    meta_df = pd.read_csv(metadata, sep='\t')
    meta_for_file = meta_df[meta_df['filename'] == os.path.basename(ss_file)]
    if meta_for_file.empty:
        print("filename {} not found in metadata".format(ss_file))
    elif len(meta_for_file) > 1:
        print("filename {} found more than once in metadata".format(ss_file))
    else:
        study = meta_for_file['study'].values[0]
        qtl_group = meta_for_file['qtl_group'].values[0]
        quant_method = meta_for_file['quant_method'].values[0]
        tissue_ont = meta_for_file['tissue_ontology_id'].values[0]
        return study, qtl_group, quant_method, tissue_ont




def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', help='The path to the summary statistics file to be processed', required=False)
    argparser.add_argument('-var', help='The path to the variant/genotype metadata file', required=False)
    argparser.add_argument('-phen', help='The path to the trait/phenotype metadata file', required=False)
    argparser.add_argument('-expr', help='The path to the gene expression file', default=None, required=False)
    argparser.add_argument('-study', help='The study identifier', required=False)
    argparser.add_argument('-qtl_group', help='The qtl group e.g. "LCL"', required=False)
    argparser.add_argument('-quant', help='The quantification method e.g. "gene counts"', required=False)
    argparser.add_argument('-tissue', help='The tissue', required=False)
    argparser.add_argument('-chr', help='The chromosome the data belongs to', required=False)
    argparser.add_argument('-loader', help='The loader type', choices=['study', 'trait', 'study_info'], default=None, required=True)
    argparser.add_argument('-tissue_ont', help='The tissue ontology term', required=False)
    argparser.add_argument('-treatment', help='The treatment', required=False)
    argparser.add_argument('-treatment_ont', help='The treatment ontology term', required=False)
    argparser.add_argument('-condition_label', help='The condition label e.g. "naive"', required=False)
    argparser.add_argument('-metadata', help='path to tsv containing metadata for the study. Columns "filename", "study", "qtl_group", "tissue_ont", "quant" must be present and will override/replace and of the cmd line arguments given. "filename" refers the basename of the summary statistics file that was passed to the "-f" parameter', required=False)


    args = argparser.parse_args()

    properties_handler.set_properties()  # pragma: no cover
    h5files_path = properties.h5files_path # pragma: no cover
    tsvfiles_path = properties.tsvfiles_path  # pragma: no cover
    database = properties.sqlite_path
    chr_dir = properties.chr_dir
    trait_dir = properties.trait_dir
    study_dir = properties.study_dir
    

    ss_file = args.f
    var_file = args.var
    phen_file = args.phen
    study = args.study
    qtl_group = args.qtl_group
    quant_method = args.quant
    tissue = args.tissue
    expr_file = args.expr
    loader_type = args.loader
    chromosome = args.chr
    tissue_ont = args.tissue_ont
    treatment = args.treatment
    treatment_ont = args.treatment_ont
    condition_label = args.condition_label
    metadata = args.metadata

    if metadata is not None:
        study, qtl_group, quant_method, tissue_ont = update_metadate_from_tsv(metadata, ss_file)

    print( study, qtl_group, quant_method, tissue_ont)

    if loader_type == 'study':
        if chromosome is None:
            print("You must specify the '-chr'...exiting")
        else:
            loader = Loader(tsv=ss_file, expr_file=expr_file,  var_file=var_file, qtl_group=qtl_group, quant_method=quant_method, study_dir=study_dir, study=study, trait=phen_file, hdf_path=h5files_path, chromosome=chromosome, sqldb=database, loader=loader_type)
            #loader = Loader(filename, tsvfiles_path, chr_dir, study_dir, study, traits, h5files_path, chromosome, database, loader_type)
            loader.load_bystudy()
    elif loader_type == 'trait':
            loader = Loader(trait=phen_file, hdf_path=h5files_path, loader=loader_type, trait_dir=trait_dir)
            loader.load_traits()
    elif loader_type == "study_info":
            loader = Loader(qtl_group=qtl_group, quant_method=quant_method, study=study, sqldb=database, trait=phen_file, loader=loader_type, tissue=tissue, tissue_ont=tissue_ont, treatment=treatment, treatment_ont=treatment_ont, condition_label=condition_label)
        #loader = Loader(filename, tsvfiles_path, chr_dir, study_dir, study, traits, h5files_path, chromosome, database, loader_type)
            loader.load_study_info()
    else:
        print("You must specify the '-loader'...exiting")



if __name__ == "__main__":
    main()
