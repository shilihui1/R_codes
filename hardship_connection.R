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
