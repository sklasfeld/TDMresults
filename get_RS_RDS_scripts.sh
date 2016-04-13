ls output/test_RS_p??_brca_log_cm.RDS > RDS_RS_content_log.sh
ls output/test_RS_p???_brca_log_cm.RDS >> RDS_RS_content_log.sh
ls output/test_RS_p??_brca_npn_cm.RDS > RDS_RS_content_npn.sh
ls output/test_RS_p???_brca_npn_cm.RDS >> RDS_RS_content_npn.sh
ls output/test_RS_p??_brca_qn_cm.RDS > RDS_RS_content_qn.sh
ls output/test_RS_p???_brca_qn_cm.RDS >> RDS_RS_content_qn.sh
ls output/test_RS_p??_brca_tdm_cm.RDS > RDS_RS_content_tdm.sh
ls output/test_RS_p???_brca_tdm_cm.RDS >> RDS_RS_content_tdm.sh
ls output/test_RS_p??_brca_un_cm.RDS > RDS_RS_content_un.sh
ls output/test_RS_p???_brca_un_cm.RDS >> RDS_RS_content_un.sh
sed -e "s/out/Rscript\ read_RDS.R\ out/g" RDS_RS_content_log.sh > 1
mv 1 RDS_RS_content_log.sh
sed -e "s/out/Rscript\ read_RDS.R\ out/g" RDS_RS_content_npn.sh > 1
mv 1 RDS_RS_content_npn.sh
sed -e "s/out/Rscript\ read_RDS.R\ out/g" RDS_RS_content_qn.sh > 1
mv 1 RDS_RS_content_qn.sh
sed -e "s/out/Rscript\ read_RDS.R\ out/g" RDS_RS_content_tdm.sh > 1
mv 1 RDS_RS_content_tdm.sh
sed -e "s/out/Rscript\ read_RDS.R\ out/g" RDS_RS_content_un.sh > 1
mv 1 RDS_RS_content_un.sh
