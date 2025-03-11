//
// Check input samplesheet and get data channels
//

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv(header:true, sep:'\t')
        .map { row -> input_channel(row) }
        .branch { m, d -> 
            assembly: m.assembly
            reads: m.reads
        }
        .set { data }

    emit:
    reads       = data.reads // channel: [ val(meta), [ reads ] ]
    assemblies  = data.assembly
}

def input_channel(LinkedHashMap row) {
    meta = [:]

    if (!row.sample) {
        exit 1, "ERROR: Please check input samplesheet -> no sample column found!\n"
    }

    meta.sample_id    = row.sample
    meta.assembly     = false
    meta.reads        = false

    // data is an assembly
    if (row.assembly) {

        meta.assembly = true

        if (!file(row.assembly).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> assembly does not exist!\n${row.assembly}"
        }

        return [ meta, file(row.assembly)]

    // data is raw reads
    } else {

        meta.reads      = true
        meta.platform   = row.platform
        meta.single_end = true


        array = []

        valid_platforms = [ 'ILLUMINA', 'NANOPORE', 'PACBIO']

        if (!valid_platforms.contains(row.platform)) {
            exit 1, "ERROR: Please check input samplesheet -> incorrect platform provided!\n${row.platform}"
        }

        if (!file(row.fq1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fq1}"
        }

        /*
        Library ID has no real function beyond folder naming, so we make it optional and
        fill lthe field with the file name otherwise
        */

        meta.library_id = row.library_id ? row.library_id : file(row.fq1).getSimpleName()

        if (row.fq2) {
            meta.single_end = false
            if (!file(row.fq2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fq2}"
            }
            array = [ meta, [ file(row.fq1), file(row.fq2) ] ]
        } else {
            array = [ meta, [ file(row.fq1)]]
        }

        return array
    }
}
