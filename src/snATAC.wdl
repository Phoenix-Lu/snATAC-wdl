version 1.0

workflow snATAC {
    # meta {
    #     author: 'Phoenix_Lu'
    #     parameter_group: {
    #         source_data: {
    #         }
    #         source_reference: {
    #         }
    #         source_parameters: {
                
    #         }

    #     }
    # }

    input {
        File origin_read1
        File origin_read2
        Boolean is_PE
        Boolean do_trim_for_5
        Array[Int]? adapter_5_len
        Boolean do_QCcheck
        String cut_5_mode
        String sample_name
    }
    

    
    call Cut_5 {
        input:
        origin_read1 = origin_read1,
        origin_read2 = origin_read2,
        adapter_5_len = adapter_5_len,
        sample_name = sample_name,
        cut_5_mode = cut_5_mode
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
    }

    command {
        python3 $(which snATAC_cut_adapter.py) \
        --r1 ${origin_read1} \
        --r2 ${origin_read2} \
        --mode ${cut_5_mode} \
        ${--prefix' ${prefix}} \ 
        --sample_name ${sample_name} \
        --len (sep(' ', adapter_5_len))
    }

    output {
        File read1_cut5 = 
        File read2_cut5 = 
    }
    runtime {
            conda: "work"
            cpu:'5'
            memory: '32' + "GB"

    }
}



task do_fastp{
    input {
        File fastq_1
        File? fastq_2
        Boolean watch_qc_only
        Boolean is_pair_end
    }
    command {
        
    }

    output {

    }
}
# task Check_sequence_QC {

# }


# task Align {

# }

# task Rm_duplicate_single {

# }

# task Check_insert_QC {

# }

# task Merge {

# }
