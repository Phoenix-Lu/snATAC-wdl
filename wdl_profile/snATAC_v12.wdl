version 1.0

workflow snATAC {
    input {
        # FILES
        String data_dir
        String SRC_DIR
        # INFOs
        String experiment_name # 大的项目名称
        String genome_size


        # SOFTWARE PATHS
        File? picard_path

        # OPTIONS
        Boolean do_cut_for_5
        Array[Int]? adapter_5_len
        Boolean do_QCcheck
        String cut_5_mode
        Boolean preserve_sam
        Boolean do_ATAC_shift


        # options
        Boolean remove_chrM
        Boolean remove_unmapped
        Boolean remove_dup
        Boolean remove_unproper_pair
        String mpQ_threhold

        # REFs
        String bowtie2_index

        File? peak_control

        File? black_list
        String? chrom_sizes



        # GLOBAL
        Boolean ATAC_DEFAULT
        Boolean ZEN_MODE
        Boolean do_QCcheck

        
        # EFFICIENCY
        Int? threads


        }


    call get_sample_json{
        input: 
        dir = data_dir,
        SRC_DIR = SRC_DIR
        }
    
    scatter (idx in range(length(get_sample_json.sample))) {
        File read1 = get_sample_json.read1[idx]
        File read2 = get_sample_json.read2[idx]
        String sample_name = get_sample_json.sample[idx]
        if (!ZEN_MODE) {
            call fastp as origin_fastp{
            input:
            read1 = read1,
            read2 = read2,
            sample_name = sample_name,
            do_trim = false,
            SRC_DIR = SRC_DIR
        }
        }

        if (do_cut_for_5) {
            call cut_5 {
                input:
                origin_read1 = read1,
                origin_read2 = read2,
                adapter_5_len = adapter_5_len,
                sample_name = sample_name,
                cut_5_mode = cut_5_mode,
                SRC_DIR = SRC_DIR

    }
    }
    File cutted_read1 = select_first([cut_5.cut_read1, read1])
    File cutted_read2 = select_first([cut_5.cut_read2, read2])
    
    call fastp as trim_fastp {
        input:
        read1 = cutted_read1,
        read2 = cutted_read2,
        sample_name = sample_name,
        do_trim = true,
        SRC_DIR = SRC_DIR
    }

    File read1_for_align = select_first([trim_fastp.trimmed_read1, read1])
    File read2_for_align = select_first([trim_fastp.trimmed_read2, read2])

    call do_align {
        input:
        read1 = read1_for_align,
        read2 = read2_for_align,
        sample_name = sample_name,
        bowtie2_index = bowtie2_index,
        picard_path = picard_path,
        preserve_sam = preserve_sam,
        SRC_DIR = SRC_DIR
    }
    call do_align_filter {
        input:
        raw_bam = do_align.align_rawbam_dupmarked,
        sample_name = sample_name,
        remove_chrM = select_first([remove_chrM, true]),
        remove_dup = select_first([remove_dup, true]),
        remove_unmapped = select_first([remove_unmapped, true]),
        remove_unproper_pair = select_first([remove_unproper_pair, true]),
        mpQ_threhold = select_first([mpQ_threhold, "30"]),
        SRC_DIR = SRC_DIR

    }

    call gen_bed_and_bigwig {
        input:
        bam_file = do_align_filter.filtered_bam,
        do_ATAC_shift = do_ATAC_shift,
        black_list_bed = black_list,
        chrom_sizes = chrom_sizes,
        genome_bed = chrom_sizes,
        SRC_DIR = SRC_DIR
    }
    call macs2_callpeaks{
        input:
        bed_file = gen_bed_and_bigwig.bed_file,
        genome_size = genome_size,
        ATAC_DEFAULT = ATAC_DEFAULT,
        peak_control = peak_control,
        SRC_DIR = SRC_DIR
    }
    if (do_QCcheck) {
    call seperate_NFR {
        input:
        bed_file = gen_bed_and_bigwig.bed_file,
        SRC_DIR = SRC_DIR
    }

    call do_final_QC {
        input:
        no_filter_bam_file = do_align.align_rawbam_dupmarked,
        accept_bam_file = do_align_filter.filtered_bam,
        align_summary = do_align.align_summary,
        picard_metric = do_align.picard_report,
        experiment_name = sample_name,
        full_bed_file = gen_bed_and_bigwig.bed_file,
        part_bed_file = [seperate_NFR.bed_NFR, seperate_NFR.bed_Nucleosome],
        SRC_DIR = SRC_DIR
    }
    }
}
    call merge_bam{
        input:
        bam_files = do_align_filter.filtered_bam,
        experiment_name = experiment_name,
        SRC_DIR = SRC_DIR

}
    call do_align_filter  as do_align_filter_bulk{
        input:
        raw_bam = merge_bam.aggregate_bam,
        sample_name = experiment_name,
        remove_chrM = select_first([remove_chrM, true]),
        remove_dup = select_first([remove_dup, true]),
        remove_unmapped = select_first([remove_unmapped, true]),
        remove_unproper_pair = select_first([remove_unproper_pair, true]),
        mpQ_threhold = select_first([mpQ_threhold, "30"]),
        SRC_DIR = SRC_DIR

    }
    call gen_bed_and_bigwig as gen_bed_and_bigwig_bulk {
        input:
        bam_file = do_align_filter_bulk.filtered_bam,
        do_ATAC_shift = do_ATAC_shift,
        black_list_bed = black_list,
        chrom_sizes = chrom_sizes,
        genome_bed = chrom_sizes,
        SRC_DIR = SRC_DIR

    }
    call macs2_callpeaks as macs2_callpeaks_bulk{
        input:
        bed_file = gen_bed_and_bigwig_bulk.bed_file,
        genome_size = genome_size,
        ATAC_DEFAULT = ATAC_DEFAULT,
        peak_control = peak_control,
        SRC_DIR = SRC_DIR

    }
    call gen_count_matrix {
        input:
        experiment_name = experiment_name,
        peak_ref_bed = macs2_callpeaks_bulk.macs2_peak_bed,
        queue_bams = do_align_filter.filtered_bam,
        SRC_DIR = SRC_DIR
    }

}

task get_sample_json{
    input {
        String dir
        String SRC_DIR
    }
    command {
        bash ${SRC_DIR}/snATAC_organize_file.sh ${dir} > sample_data.json
    }
    output {
        Array[String] sample = read_json("sample_data.json").sample
        Array[File] read1 = read_json("sample_data.json").Read1
        Array[File] read2 = read_json("sample_data.json").Read2
    }
}





task fastp {
    input{
        File read1
        File read2
        String sample_name
        String? prefix
        Boolean do_trim
        String SRC_DIR        
    }

    command{
        python3 ${SRC_DIR}/snATAC_do_fastp.py \
        --r1 ${read1} \
        ${'--r2 ' + read2} \
        --sample_name ${sample_name} \
        ${if do_trim then '--do_trim ' else ''} \
        ${if defined(prefix) then prefix else ''} \
        --sample_name ${sample_name}
        }



    output{
        File trimmed_read1 = if do_trim then
         glob("*_r1_trimed.fq*")[0] else read1
        File trimmed_read2 = if do_trim then
         glob("*_r2_trimed.fq*")[0] else read2
        File qc_html_report = glob("*fastp_report.html")[0]
    }
}


task cut_5 {
    input {
        File origin_read1
        File? origin_read2
        Array[Int]? adapter_5_len
        String cut_5_mode
        String? prefix
        String sample_name
        String SRC_DIR
    }

    command {
        python3 ${SRC_DIR}/snATAC_cut_strand5.py \
        --r1 ${origin_read1} \
        --r2 ${origin_read2} \
        --mode ${cut_5_mode} \
        ${if defined(prefix) then prefix else ''} \
        --sample_name ${sample_name} \
        --len ${sep=' ' adapter_5_len}
    }
    output {
        File? cut_read1 = glob("*_r1_*_cutted.fq*")[0]
        File? cut_read2 = glob("*_r2_*_cutted.fq*")[0]
    }
    }

task do_align {
    input{
    File read1
    File read2
    String bowtie2_index
    File? picard_path
    Boolean preserve_sam
    String sample_name
    String SRC_DIR
    }

    command {
        touch ./null
        python3 ${SRC_DIR}/snATAC_align.py \
        --r1 ${read1} \
        --r2 ${read2} \
        --sample_name ${sample_name} \
        --reference_path ${bowtie2_index} \
        ${'--picard_path ' + picard_path} \
        ${if preserve_sam then '--preserve_sam ' else ''}

    }

    output {
        
        File? align_sam = if preserve_sam then glob("*.sam")[0] else 
        "./null"
        # File
        File align_rawbam_dupmarked = glob("*dup.bam")[0]
        
        # Report
        File align_summary = glob("*.alignsummary")[0]
        File picard_report = glob("*.picardDupsummary")[0]
    }
}

task do_align_filter {
    input {
        # Files
        File raw_bam

        # Infos
        String sample_name

        # options
        String mpQ_threhold
        Boolean remove_chrM
        Boolean remove_unmapped
        Boolean remove_unproper_pair
        Boolean remove_dup
        String SRC_DIR
    
    }
    
    
    command {
    python3 ${SRC_DIR}/snATAC_align_filter.py \
    --sample_name ${sample_name} \
    --bam_file ${raw_bam} \
    --mpQ ${mpQ_threhold} \
    ${if remove_dup then "--rmdup " else ''} \
    ${if remove_unproper_pair then "--rmunpro_pair " else ''} \
    ${if remove_chrM then "--rmchrM " else ''} \
    ${if remove_unmapped then "--rmunmap " else ''}
    }
    
    output {
        # File filter_summary = glob("*.filtersummary")[0]

        File filtered_bam = glob("*.bam")[0]
        File filtered_bam_bai = glob("*.bai")[0]
    }

}

task gen_bed_and_bigwig {
    input {
        File bam_file
        # options
        Boolean do_ATAC_shift
        # refs
        File? black_list_bed
        File? chrom_sizes
        File? genome_bed

        String SRC_DIR

    }
    

    command {
        touch ./null
        python3 ${SRC_DIR}/snATAC_bam2bed.py \
        --bam_file ${bam_file} \
        ${'--black_list ' + black_list_bed} \
        ${if do_ATAC_shift then "--do_atac_shift " else ''} \
        ${'--chrom_sizes ' + chrom_sizes} \
        ${'--genome_bed ' + genome_bed}
    }

    output {
        File bed_file = glob("*.bed")[0]
        File bedgraph_file = if defined(chrom_sizes) then glob("*.bedgraph")[0] else "./null"
        File bw_file = if defined(genome_bed) then glob("*.bw")[0] else "./null"
    }
}

task macs2_callpeaks {
    input {
        File bed_file
        String genome_size
        Boolean? ATAC_DEFAULT
        File? peak_control
        String SRC_DIR

    }

    command {
        python3 ${SRC_DIR}/snATAC_callpeaks.py \
        --bed_file ${bed_file} \
        --genome_size ${genome_size} \
        ${if defined(ATAC_DEFAULT) then '--ATAC_default_mode' else ''} \
        ${'--control' + peak_control}
    }

    output {
        File macs2_peak_narrowPeak = glob("*.narrowPeak")[0]
        File macs2_peak_xls = glob("*.xls")[0]
        File macs2_peak_bed = glob("*.bed")[0]

    }
}


task merge_bam{
    input {
        Array[File] bam_files
        String experiment_name
        String SRC_DIR
    }
    command {
        python3 ${SRC_DIR}/snATAC_merge_cells.py \
        --bam_files ${sep=' ' bam_files} \
        --experiment_name ${experiment_name}
    }
    output {
        File aggregate_bam = glob('*.bam')[0]
    }
}


task gen_count_matrix{
    input{
        String experiment_name
        File peak_ref_bed
        Array[File] queue_bams
        String SRC_DIR
    }

    command{
        python3 ${SRC_DIR}/snATAC_calculate_counts.py \
        --experiment_name ${experiment_name} \
        --peak_reference_bed ${peak_ref_bed} \
        --queue_bam ${sep=' ' queue_bams}
    }

    output{
        File count_matrix = glob("*.mtx")[0]
        File count_row = glob("*.rows")[0]
        File count_col = glob("*.cols")[0]

    }
}


task seperate_NFR {
    input {
        File bed_file
        String SRC_DIR
    }
    command {
        bash ${SRC_DIR}/snATAC_extract_NFR_and_Nucleosome.sh ${bed_file}

    }
    output {
        File bed_NFR = glob("*_NFR.bed")[0]
        File bed_Nucleosome = glob("*_Nucleosome.bed")[0]
    }
}

task do_final_QC {
    input {
        File no_filter_bam_file
        File accept_bam_file
        File align_summary
        File picard_metric
        String experiment_name
        File full_bed_file
        Array[File] part_bed_file
        String SRC_DIR
    }

    command {
        python3 ${SRC_DIR}/snATAC_QC_function.py \
        --no_filter_bam_file ${no_filter_bam_file} \
        --accept_bam_file ${accept_bam_file} \
        --align_summary ${align_summary} \
        --picard_metric ${picard_metric} \
        --sample_name ${experiment_name} \
        --full_bed_file ${full_bed_file} \
        --part_bed_file ${sep=' ' part_bed_file}

    }

    output {
        File final_QC_report = glob("QCs/*.html")[0]
    }

}


