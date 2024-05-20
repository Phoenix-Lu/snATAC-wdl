# snATAC-wdl



# Step1: get dir

    cp /work/xulab/xulab-seq/lsr/PhoenixPipeline/snATAC-release YOURPATH
    or 
    clone this repository it from github

# Step2: set input_json

This is an example for wdl json profile
Since no comments are allowed in json, we explain the agrs here:

    # MUST set!
    ** in snATAC_v12.json:

    {
    # 必须设置的参数
    "snATAC.SRC_DIR": "/work/xulab/lushaorong_work/snATAC-release/src",
    "snATAC.picard_path": null,
    "snATAC.bowtie2_index": "/work/xulab/xulab-seq/reference/Bowtie2/Human/hg38/index/GRCh38_noalt_as/GRCh38_noalt_as",
    "snATAC.black_list": "/work/xulab/xulab-seq/reference/annotation/encode/blacklist/hg38-blacklist.v2.bed",

    # 基本数据相关参数
    "snATAC.data_dir": "/work/xulab/lushaorong_work/snATAC-release/snATAC/raw",
    "snATAC.experiment_name": "THP1_test",
    "snATAC.do_QCcheck": true,

    # 是否要进行5’处理
    "snATAC.do_cut_for_5": true,
    "snATAC.adapter_5_len": [34], # 这是一个Array，最多接受两个数值，分别对应R1和R2
    "snATAC.cut_5_mode": "R1",  # 可选:"R1","R2","R1R2"

    # align和filter参数，
    "snATAC.preserve_sam": true
    "snATAC.remove_dup": true,
    "snATAC.remove_unproper_pair": true,
    "snATAC.remove_chrM": true,
    "snATAC.remove_unmapped": true,
    "snATAC.mpQ_threhold": "30",

    # bed文件和macs2参数
    "snATAC.do_ATAC_shift": true,
    "snATAC.ATAC_DEFAULT": true,
    "snATAC.peak_control": null,
    "snATAC.genome_size": "hs",

    # 不推荐改动
    "snATAC.threads": null,
    "snATAC.ZEN_MODE": false,
    "snATAC.chrom_sizes": null,
    "snATAC.trim_fastp.prefix": null,
    "snATAC.origin_fastp.prefix": null,
    "snATAC.cut_5.prefix": null
    }





# Step3:配置环境。

    # croo 环境建议单独设置!!!
    conda create -n croo;
    conda activate croo;
    conda install python=3.11;
    pip install croo;
    conda install graphviz;

通常来说，应该和大家的环境差一个caper 和 fastp，大家单独在自己环境中装一个conda install caper, fastp 应该就可以运行。我的环境如下

    # 运行环境配置文件:* snATAC-release/ATAC.yaml

**IMPORTANT**: Read Caper's [README](https://github.com/ENCODE-DCC/caper/blob/master/README.md) carefully to choose a backend for your system. Follow the instruction in the configuration file.
	```bash
	# backend: local or your HPC type (e.g. slurm, sge, pbs, lsf). read Caper's README carefully.
	$ caper init [YOUR_BACKEND]

	# IMPORTANT: edit the conf file and follow commented instructions in there
	$ vi ~/.caper/default.conf
	```

# Step4:运行
    
    # 运行wdl, caper运行路径就是产生结果的dir
    $ conda activate ENV
	$ caper run ${WDL_FILE} -i "${INPUT_JSON}" 
    
    # 运行croo 组织输出
    $ conda activate croo
    $ croo ${metadata.json} --out-def-json snATAC.json --out-dir ./croo_out

