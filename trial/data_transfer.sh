#!/usr/bin/bash


aws s3 sync genesis_data_share s3://genesis-share-actrec-tata/
while [ $? -ne 0 ]; do
	#aws s3 sync genesis_data_share s3://genesis-share-actrec-tata/
	#echo "transferring data"
	aws s3 sync tmp s3://genesis-share-actrec-tata/
done

