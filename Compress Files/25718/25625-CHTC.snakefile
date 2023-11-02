import os
import pysam
import pandas as pd
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# this workflow inverts the paradigm I've used previously
# it takes a FASTA file of IPD MHC exon 2 reference sequences (created in this experiment)
# and immunoWES data (WGS should work too, but slower)
# high quality FASTA immunoWES reads are extracted from the total immunoWES dataset
# these are assembled with SPAdes to form contigs that are slightly longer than exon 2
# contigs that contain IPD exon 2 sequences provide genotypes that are written to a labkey table
# 
# after genotyping, orphaned contigs from the exact match search
# and contigs formed from reads that do not exactly match MHC sequences in IPD
# are combined into a file of potential exon 2 sequences that could be named in the futre

## Config ##

sample_name = config['sample_name']
r1_fastq = config['r1_fastq']
r2_fastq = config['r2_fastq']
ref_fasta = config['ref_fasta']
experiment_number = config['experiment_number']

## Functions ##

def find_matching_sequences(experiment_number, immunowes_fasta, ref_fasta, sample_name, output_csv, output_unused_fasta, output_used_fasta):
    '''
    given an MES contig FASTA file and a reference FASTA file
    create a list containing MES contig sequences that fully
    contain reference FASTA sequences.

    also return FASTA of sequences that do not match any in IPD
    '''

    # create dictionary of sequences in reference FASTA
    ref_dict = {}
    
    for seq_record in SeqIO.parse(ref_fasta, "fasta"):
        ref_dict.update({seq_record.id : str(seq_record.seq)})

    # create dictionary of sequences in MES FASTA
    mes_dict = {}
    
    for seq_record in SeqIO.parse(immunowes_fasta, "fasta"):
        mes_dict.update({seq_record.id : str(seq_record.seq)})

    # create genotyping list
    genotypes = []

    # create FASTA file for immunoWES reads that do not include IPD exact matches
    with open(output_unused_fasta, "w") as output_unused_handle, open(output_used_fasta, "w") as output_used_handle:

        # iterate over MES FASTA sequences
        for mes_key, mes_val in mes_dict.items():
            # set counter to zero initially
            # if still zero after iterating through all reference sequences write sequence to output FASTA
            mes_used = 0

            # iterate over reference FASTA sequences
            for ref_key, ref_val in ref_dict.items():

                # test if MES FASTA contains reference FASTA
                if ref_val in mes_val:

                    # save values to add to list
                    genotypes.append([str(experiment_number), sample_name, ref_key, mes_key, mes_val, ref_fasta])

                    # add to output_used_fasta
                    # create biopython object
                    record = SeqRecord(
                    Seq(mes_val),
                    id=mes_key,
                    description=''
                    )
                    SeqIO.write(record, output_used_handle, "fasta")

                    # change flag so sequence isn't added again to unused
                    mes_used += 1
            
            if mes_used == 0: 
                # create biopython object
                record = SeqRecord(
                Seq(mes_val),
                id=mes_key,
                description=''
                )
                SeqIO.write(record, output_unused_handle, "fasta")

    # create CSV from genotypes
    with open(output_csv, 'w') as f:
      
        # using csv.writer method from CSV package
        write = csv.writer(f)
        
        write.writerow(['experiment_number', 'animal_id', 'genotype', 'contig_name', 'contig_sequence', 'ref_db'])
        write.writerows(genotypes)

    return(genotypes)   

rule all:
    input:
        # genotype CSV that is written to labkey
        sample_name + '.' + ref_fasta + '.genotypes.csv',

        # FASTA of SPAdes contigs used in genotyping CSV
        sample_name + '.' + ref_fasta + '.used_exact.fasta',

        # FASTA of SPAdes contigs from novel sequences plus contigs unused in genotyping
        sample_name + '.' + ref_fasta + '.novel.fasta',
        sample_name + '.labkey.tkn'
    run:
        # remove temporary labkey file
        shell('rm {input[3]}')

rule map_to_reference_hla:
    '''map reads to reference HLA'''
    input:
        ref_fasta,
        r1_fastq,
        r2_fastq   
    output:
        temp(sample_name + '.' + ref_fasta + '.aln.R1.fastq.gz'),
        temp(sample_name + '.' + ref_fasta + '.aln.R2.fastq.gz'),
    run:
        # map with bbmap in semiperfect mode
        # this will capture all hits to IPD FASTA sequences
        # since we are not performing allele discovery, all matches should exist in IPD already
        shell('bbmap.sh \
            ref={input[0]} \
            in={input[1]} \
            in2={input[2]} \
            outm={output[0]} \
            outm2={output[1]} \
            semiperfectmode=t ambig=random -Xmx3g')

rule filter_hq_exact_reads:
    '''
    vsearch to save only highest quality reads
    '''
    input:
        sample_name + '.' + ref_fasta + '.aln.R1.fastq.gz',
        sample_name + '.' + ref_fasta + '.aln.R2.fastq.gz',
    output:
        temp(sample_name + '.' + ref_fasta + '.hq.R1.fastq'),
        temp(sample_name + '.' + ref_fasta + '.hq.R2.fastq')  
    run:
        shell('vsearch \
        --fastx_filter {input[0]} \
        --reverse {input[1]} \
        --fastqout {output[0]} \
        --fastqout_rev {output[1]} \
        -fastq_maxee 1.0')

rule assemble_exact_fastq:
    '''
    run SPAdes in multi-cell mode
    tested other modes in Geneious, this mode works best
    '''
    input:
        sample_name + '.' + ref_fasta + '.hq.R1.fastq',
        sample_name + '.' + ref_fasta + '.hq.R2.fastq' 
    output:
        sample_name + '.' + ref_fasta + '.spades.fasta'
    params:
        SN = sample_name,
        LOCUS = ref_fasta
    run:
        # run SPAdes
        # use k=127 because this worked better in testing
        shell('spades.py \
            --pe1-1 {input[0]} \
            --pe1-2 {input[1]} \
            --isolate \
            --phred-offset 33 \
            -k 127 \
            -o {params.SN}.{params.LOCUS}')

        # copy to output location
        shell('cp {params.SN}.{params.LOCUS}/scaffolds.fasta \
        {output[0]}')

        # delete spades directory
        shell('rm -rf {params.SN}.{params.LOCUS}')

rule vsearch__exact_orient:
    '''
    orient sequences relative to trim_fasta
    otherwise genotyping doesn't work correctly
    '''
    input:
        ref_fasta,
        sample_name + '.' + ref_fasta + '.spades.fasta'
    output:
        temp(sample_name + '.' + ref_fasta + '.exact.oriented.fasta')
    run:
        shell('vsearch --orient \
        {input[1]} \
        --db {input[0]} \
        --fastaout {output[0]}')

rule rename_exact:
    '''
    rename FASTA sequences to include sample name
    '''
    input:
        sample_name + '.' + ref_fasta + '.exact.oriented.fasta'
    output:
        temp(sample_name + '.' + ref_fasta + '.exact.fasta')
    run:
        shell('rename.sh \
        in={input[0]} \
        out={output[0]} \
        addprefix=t -Xmx3g \
        prefix=' + sample_name)

rule genometyping:
    '''
    find all reference sequences fully contained by IPD sequences
    '''
    input:
        sample_name + '.' + ref_fasta + '.exact.fasta',
        ref_fasta
    output:
        sample_name + '.' + ref_fasta + '.genotypes.csv',
        temp(sample_name + '.' + ref_fasta + '.unused_exact.fasta'),
        sample_name + '.' + ref_fasta + '.used_exact.fasta'
    run:
        find_matching_sequences(experiment_number, input[0], input[1], sample_name, output[0], output[1], output[2])

rule remove_genometyped_sequences:
    '''
    this begins a second pass analysis that will find additional exon 2 sequences
    that are not in the IPD database

    first remove sequences that exactly match known genometypes
    '''
    input:
        sample_name + '.' + ref_fasta + '.exact.fasta',
        r1_fastq,
        r2_fastq   
    output:
        temp(sample_name + '.' + ref_fasta + '.unmapped.R1.fastq.gz'),
        temp(sample_name + '.' + ref_fasta + '.unmapped.R2.fastq.gz'),
    run:
        shell('bbmap.sh \
            ref={input[0]} \
            in={input[1]} \
            in2={input[2]} \
            outu={output[0]} \
            outu2={output[1]} \
            semiperfectmode=t ambig=random -Xmx3g')

rule collect_novel_exon_2_sequences:
    '''
    find MHC exon 2 sequneces in reads that did not match genometyping results

    this allows mismatches relative to reference sequences, since these are potentially novel alleles
    '''
    input: 
        ref_fasta,
        sample_name + '.' + ref_fasta + '.unmapped.R1.fastq.gz',
        sample_name + '.' + ref_fasta + '.unmapped.R2.fastq.gz'
    output:
        temp(sample_name + '.' + ref_fasta + '.novel.aln.R1.fastq.gz'),
        temp(sample_name + '.' + ref_fasta + '.novel.aln.R2.fastq.gz'),  
    run:
        shell('bbmap.sh \
            ref={input[0]} \
            in={input[1]} \
            in2={input[2]} \
            outm={output[0]} \
            outm2={output[1]} \
            ambig=random -Xmx3g')

rule filter_hq_novel_reads:
    '''
    vsearch to save only highest quality reads
    '''
    input:
        sample_name + '.' + ref_fasta + '.novel.aln.R1.fastq.gz',
        sample_name + '.' + ref_fasta + '.novel.aln.R2.fastq.gz', 
    output:
        temp(sample_name + '.' + ref_fasta + '.novel.hq.R1.fastq'),
        temp(sample_name + '.' + ref_fasta + '.novel.hq.R2.fastq')  
    run:
        shell('vsearch \
        --fastx_filter {input[0]} \
        --reverse {input[1]} \
        --fastqout {output[0]} \
        --fastqout_rev {output[1]} \
        -fastq_maxee 1.0')

rule assemble_novel_fastq:
    '''
    run SPAdes in multi-cell mode
    '''
    input:
        sample_name + '.' + ref_fasta + '.novel.hq.R1.fastq',
        sample_name + '.' + ref_fasta + '.novel.hq.R2.fastq'  
    output:
        temp(sample_name + '.' + ref_fasta + '.novel.spades.fasta')
    params:
        SN = sample_name,
        LOCUS = ref_fasta
    run:
        # run SPAdes
        # use k=127 because this worked better in testing
        shell('spades.py \
            --pe1-1 {input[0]} \
            --pe1-2 {input[1]} \
            --isolate \
            --phred-offset 33 \
            -k 127 \
            -o {params.SN}.{params.LOCUS}')

        # copy to output location
        shell('cp {params.SN}.{params.LOCUS}/scaffolds.fasta \
        {output[0]}')

        # delete spades directory
        shell('rm -rf {params.SN}.{params.LOCUS}')

rule vsearch__novel_orient:
    '''
    orient sequences relative to trim_fasta
    otherwise genotyping doesn't work correctly
    '''
    input:
        ref_fasta,
        sample_name + '.' + ref_fasta + '.novel.spades.fasta'
    output:
        temp(sample_name + '.' + ref_fasta + '.novel.oriented.fasta')
    run:
        shell('vsearch --orient \
        {input[1]} \
        --db {input[0]} \
        --fastaout {output[0]}')

rule rename_novel:
    '''
    rename FASTA sequences to include sample name
    '''
    input:
        sample_name + '.' + ref_fasta + '.novel.oriented.fasta'
    output:
        temp(sample_name + '.' + ref_fasta + '.novel.renamed.fasta')
    run:
        shell('rename.sh \
        in={input[0]} \
        out={output[0]} \
        addprefix=t -Xmx3g \
        prefix=' + sample_name)

rule combine_unused:
    '''
    there are two types of novel sequences:
    1. Exact search contigs that do not match IPD sequences
    2. Contigs from novel assembly

    Merge together into a final novel FASTA file that can be used for novel
    sequence discovery
    '''
    input:
        sample_name + '.' + ref_fasta + '.novel.renamed.fasta',
        sample_name + '.' + ref_fasta + '.unused_exact.fasta'
    output:
        temp(sample_name + '.' + ref_fasta + '.novel.merged.fasta')
    run:
        shell('cat {input[0]} {input[1]} > {output[0]}')

rule remove_near_identical:
    '''
    only retain potential "novel" sequences that are more than 1% divergent
    from sequences that are exact matches to known IPD sequences
    '''
    input:
        sample_name + '.' + ref_fasta + '.used_exact.fasta',
        sample_name + '.' + ref_fasta + '.novel.merged.fasta'
    output:
        sample_name + '.' + ref_fasta + '.novel.fasta'
    run:
        shell('mapPacBio.sh \
        ref={input[0]} \
        in={input[1]} \
        minid=0.99 \
        outu={output[0]}')

rule labkey_import:
    '''import high quality (hq) genotypes into labkey immunogenome table
    https://dholk.primate.wisc.edu/list/dho/projects/genomics/immunogenomics/grid.view?listId=3'''
    input:
        sample_name + '.' + ref_fasta + '.genotypes.csv'
    output:
        sample_name + '.labkey.tkn'
    run:       
        import pandas as pd
        import datetime as dt
        from pathlib import Path
        import os
        from labkey.api_wrapper import APIWrapper

        print("Create an APIWrapper")
        labkey_server = 'dholk.primate.wisc.edu'
        project_name = 'dho/projects/genomics/immunogenomics'  # Project folder name
        schema = 'lists'
        table = 'immunoWES'
        api = APIWrapper(labkey_server, project_name, api_key='apikey|4e80e6d4f9099b6d70f6fbbffc071e46', use_ssl=True)

        # import genotypes as dataframe
        # only import if input is not empty
        
        if os.stat(input[0]).st_size > 0:
            # create dataframe
            df = pd.read_csv(input[0], usecols=[0, 1, 2, 3, 4, 5], error_bad_lines=False)

            # rename columns to use labkey immunogenome naming convention
            df.columns = ['experiment_number', 
                        'animal_id',
                        'genotype',
                        'contig_name',
                        'contig_sequence',
                        'ref_db']

            # add column with current time
            df['genotype_date'] = dt.datetime.today().strftime("%m/%d/%Y")
            labkey_data = df.to_dict('records')

            # import dataframe to labkey
            result = api.query.insert_rows(schema, table, labkey_data)
        
        # create output file for snakemake
        Path(output[0]).touch()