cd /home/svenja/Price_genomes_GCF

cat Price_127.txt | while read -r acc ; do
    esearch -db assembly -query $acc </dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_RefSeq \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget "$url/$fname" ;
        done ;
    done

mmv -- '*_*_*_*.genomic.fna.gz' '#1_#2.genomic.fna.gz'

