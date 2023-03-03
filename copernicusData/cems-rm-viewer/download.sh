#!/bin/bash
for i in $(seq 573 1 619)
do
   aws s3 sync s3://cems-rm-viewer/EMSR$i ./EMSR$i
done
