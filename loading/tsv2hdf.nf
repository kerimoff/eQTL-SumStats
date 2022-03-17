// nextflow run tsv2hdf.nf

/*
================================================================================
   Convert eQTL summary statistics to HDF5 for using in the eQTL SumStats API
================================================================================
*/


// Input files _must_ be of the *.tsv.gz residing in the tsv_in directory
def tsv_glob = []
params.quant_methods.each {
                tsv_glob.add(new File(params.tsv_in, "*_$it*.tsv.gz"))
        }
tsv_to_process = Channel.fromPath(tsv_glob).view()
meta_table_ch_1 = Channel.fromPath(params.meta_table)
// meta_table_ch_2 = Channel.fromPath(params.meta_table)

/* Any previously generated HDF5 files in the hdf5_study_dir will be included
   in the chromosome + quant_method files.
*/


/*
================================================================================
              Split tsvs by chromosome, convert to hdf5 and index
================================================================================
*/

process study_tsv_to_hdf5 {
  scratch true
  stageInMode "copy"
  // containerOptions "--bind $params.tsv_in --bind $params.meta_table --bind $params.hdf5_study_dir --bind $params.hdf5_chrom_dir"
  // storeDir "/gpfs/hpc/projects/eQTLCatalogue/HDF5"

  input:
  each chr from params.chromosomes
  file tsv from tsv_to_process
  file meta_table from meta_table_ch_1.collect()

  output:
  tuple val(chr), file("${chr}/*.h5") into study_ch

  """
  mkdir $chr;
  eqtl-load -f $tsv -metadata $meta_table -chr $chr -loader study;
  """
}


/*
================================================================================
Consolidate all chromosome + quant method combinations into their own HDF5 files
================================================================================
*/

// process consolidate_hdfs_by_chrom {
//   // publishDir "/gpfs/space/home/kerimov/HDF5_data_conv/loading/data/hdf/consolidate_hdfs_by_chrom", mode: 'move'

//   input:
//   each method from params.quant_methods
//   tuple val(chr), file(per_chr) from study
//   file meta_table from meta_table_ch_2.collect()

//   output:
//   tuple val("${per_chr.baseName}"), file("${chr}.${method}.h5") into hdf5_chrom

//   script:
//   """
//   eqtl-consolidate -in_file $per_chr -out_file ${chr}.${method}.h5 -meta $meta_table -quant $method -chrom $chr
//   """

//   // script:
//   // """
//   // hdf="${chr}.${method}.h5"
//   // out="file_${chr}.${method}.h5"
//   // for f in $params.hdf5_study_dir/$chr/*+"$method".h5; do
//   //       if [ -f \$f ]; then
//   //               echo in_if
//   //               eqtl-consolidate -in_file \$f -out_file \$hdf -meta $meta_table -quant $method -chrom $chr
//   //       fi
//   // done
//   // """
// }

/*
================================================================================
                        Index the consolodated hdf5
================================================================================
*/

process index_consolidated_hdfs {
  publishDir "$params.hdf5_study_dir", mode: 'move'

  input:
  tuple val(chr), file(hdf5) from study_ch

  output:
  file "chr${chr}_${hdf5.baseName - 'file_'}.h5" 

  script:
  """
  eqtl-reindex -f $hdf5
  ptrepack --chunkshape=auto --propindexes --complevel=9 --complib=blosc $hdf5 chr${chr}_${hdf5.baseName - 'file_'}.h5
  """

}


