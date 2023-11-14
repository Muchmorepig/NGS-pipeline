#!/usr/bin/bash
get_info() {
  case "$1" in
    tair10)
      idx="/data5/wmc_data/index/bowtie2/Tair10/tair10"
      s=1.1e8
      ;;
    irgsp)
      idx="/data5/wmc_data/index/bowtie2/Oryza_sativa_RAP/Oryza_sativa"
      s=3.7e8
      ;;
    mpTak1)
    idx="/data5/wmc_data/index/bowtie2/Mpolymorpha/V5/Mpolymorpha_V5"
    s=2.2e8
        ;;
    mpTak2)
    idx="/data5/wmc_data/index/bowtie2/Mpolymorpha/V6/Mpolymorpha_V6"
    s=2.4e8
        ;;
    *)
      echo -e "\n\n\033[0;31mUnknown genome. Please Check the Usage\033[0m"
      exit
      ;;
  esac

  export bt2idx=$idx
  export gs=$s
}