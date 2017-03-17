# PROFQ - Processing FASTQ for format conversion or demultiplexing

# USAGE

SYNOPSIS
       FASTQ_Format_Conversion_v4.pl [options]

OPTIONS

    --task
        (Required.) User specified task to perform. Input options are
        "changeidformat" OR "makeread2umi" OR "demultiplex". Reform to
        script manual (--man option) for details.

    --read1fastq
        (Required.) Full path to the read 1 FASTQ file.

    --outputfolder
        (Required.) Full path to the output folder where the outfile should
        be written.

    --uniquestr
        (Optional.) Experiment name, which will be used to generate output
        file names. DEFAULT: "Processed"

    --read2fastq
        (Required Only for paired end data) Full path to the read 2 FASTQ
        file. Required if working with paired end data. Ignore if single end
        data.

    --barcodelength
        (Optional.) Length of the barcode sequenced. ONLY Required for the
        "makeread2umi" task. DEFAULT: 6

    --flowcellid
        (Optional.) Flowcell ID on which the sample was sequenced. This is
        usually present in the run folder name generated by the sequencer OR
        some FASTQ file read ID formats. Optional. DEFAULT: "Unknown"

    --listofbarcodes
        (Required Only for demultiplex task) List of barcodes separate by a
        comma. EXAMPLE: TAGTGC,ACGTGA,TGACGT

    --options
        Prints a brief help message and exits.

    --help
        Print a brief help message and exits.

    --man
        Prints a brief description and help message.

DESCRIPTION
       This program can perform the following 2 tasks:

       1. "changeidformat" - Converts FASTQ read header/indentifier format

               Automatically identifies the format of the read id and changes it between the 2 following 2 formats
               read id format 1: @HWUSI-EAS100R:6:73:941:1973#ATCACG/1
               read id format 2: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

               FASTQ files from WI Genome core end with ";1" or ";2" in format 1 specified above - This script can handle this variation

       2. "makeread2umi" - Extract the UMI sequence from the index read and put it in a separate FASTQ file. The output can be used to identify duplicates using the UMI FASTQ file and alignment (SAM) file for read 1.

       3. "demultiplex" - Demultiplex single FASTQ file into multiple FASTQ files based the expected barcodes to be sequenced.

