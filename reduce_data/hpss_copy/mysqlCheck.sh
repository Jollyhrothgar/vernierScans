#!/bin/sh
#mysql -hphnxdata1.rhic.bnl.gov -uphenix_c_user -pBrass_Ring phenix_carousel --execute="select User, Status from Entries where  User='abhisek';"
mysql -hphnxdata1.rhic.bnl.gov -uphenix_c_user -pBrass_Ring phenix_carousel --execute="select * from Entries where  User='beaumim';"
