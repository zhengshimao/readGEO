[[ -f ./GSE106172/matrix/GSE106172_series_matrix.txt.gz ]] && echo "OK: ./GSE106172/matrix/GSE106172_series_matrix.txt.gz" || { mkdir -p ./GSE106172/matrix/ && ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE106nnn/GSE106172//matrix/GSE106172_series_matrix.txt.gz// ./GSE106172/matrix/ ; }
[[ -f ./GSE106172/miniml/GSE106172_family.xml.tgz ]] && echo "OK: ./GSE106172/miniml/GSE106172_family.xml.tgz" || { mkdir -p ./GSE106172/miniml/ && ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE106nnn/GSE106172//miniml/GSE106172_family.xml.tgz// ./GSE106172/miniml/ ; }
[[ -f ./GSE106172/soft/GSE106172_family.soft.gz ]] && echo "OK: ./GSE106172/soft/GSE106172_family.soft.gz" || { mkdir -p ./GSE106172/soft/ && ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE106nnn/GSE106172//soft/GSE106172_family.soft.gz// ./GSE106172/soft/ ; }
[[ -f ./GSE106172/suppl/GSE106172_RAW.tar ]] && echo "OK: ./GSE106172/suppl/GSE106172_RAW.tar" || { mkdir -p ./GSE106172/suppl/ && ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE106nnn/GSE106172//suppl/GSE106172_RAW.tar// ./GSE106172/suppl/ ; }
[[ -f ./GSE106172/suppl/filelist.txt ]] && echo "OK: ./GSE106172/suppl/filelist.txt" || { mkdir -p ./GSE106172/suppl/ && ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/geo/series/GSE106nnn/GSE106172//suppl/filelist.txt// ./GSE106172/suppl/ ; }