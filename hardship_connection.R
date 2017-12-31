library(RODBC)

dbcname <- "aipprod.unix.anz"
database <- "PRD_CAA_CRE_DDWSP_PI7_Kung"
uid <- "shil5"
tdwallet <- "$tdwallet(shil5_pwd)"

in_data_conn_str <- paste0("DRIVER=Teradata;DBCNAME=",database,";UID=",uid,";PWD=",tdwallet,";AUTHENTICATION=LDAP")

channel <- odbcDriverConnect(in_data_conn_str) 

data_in_table <- "PRD_CAA_CRE_DDWSP_PI7_KUNG.POPU_CUST_INFO_V2"
sql <- paste0("select top 10* from", data_in_table)
scoringTable <- sqlQuerry(channel, sql)
head(scoringTable)
dim(scoringTable)

tbQuery <- "
select y.* from
(
select x.*, row_number() over (partition by cust_id order by bad desc) as row_no from
(
select T3.CUST_ID as cust_id,
T4.ACCT_ID,
T1.GNDR_CD_BMAP,
T1.MRTL_STAT_CD_BMAP,
T1.MRTL_STAT_CD_SRC,
T2.OCPTN_CD_BMAP,
extract(Year from CURRENT_DATE) - extract(Year from T1.BRTH_DT) as AGE,
T5.BLCK_CD_1_CD_SRC,
T5.BLCK_CD_2_CD_SRC,
(
case when T5.BLCK_CD_1_CD_SRC in ('Y','Z','A','B','I','Q','P','W','R')
or T5.BLCK_CD_2_CD_SRC in ('Y','Z','A','B','I','Q','P','W','R') then 1 else 0 end) as BAD
from PRD_ADS_IL_VR.VR_S_INDVDL_CUST T1
join PRD_ADS_IL_VR.VR_S_CUST_MAIN T2
on T1.CUST_ID = T2.CUST_ID
join PRD_ADS_IL_VR.VR_H_CUST T3
on T1.CUST_ID = T3.CUST_ID
join PRD_ADS_IL_VR.VR_L_CUST_ACCT_RLSHP T4
on T1.CUST_ID = T4.CUST_ID
join PRD_ADS_IL_VR.VR_S_ACCT_CRD T5
on T1.CUST_ID = T5.CUST_ID
where T1.END_DATE = '9999-12-31' (DATE)
and T2.END_DATE = '9999-12-31' (DATE)
and T4.END_DATE = '9999-12-31' (DATE)
and T5.END_DATE = '9999-12-31' (DATE)
and T1.RECORD_DELETED_FLAG = 0
and T2.RECORD_DELETED_FLAG = 0
and T3.RECORD_DELETED_FLAG = 0
and T4.RECORD_DELETED_FLAG = 0
) X
) Y
where Y.ROW_NO = 1
"
