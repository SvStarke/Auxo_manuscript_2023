cd /mnt/nuuk/2023/Validation_isolates/text_file

cat Isolates_Validation.txt | while read -r acc ; do
    esearch -db assembly -query $acc </dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_RefSeq \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget "$url/$fname" ;
        done ;
    done

mmv -- '*_*_*_genomic.fna.gz' '#1_#2.genomic.fna.gz'


