#!/bin/bash

# script to compile and open the .dot GraphViz Dot script output for the flowcharts

# local dot: dot - graphviz version 2.40.1 (20161225.0304)
# server dot: dot - graphviz version 2.26.0 (20091210.2329)

# download:
# http://www.graphviz.org/Download_source.php
# https://github.com/ellson/MOTHBALLED-graphviz
# https://gitlab.com/graphviz/graphviz
# wget http://www.graphviz.org/pub/graphviz/stable/SOURCES/graphviz-2.40.1.tar.gz
# tar -zxvf graphviz-2.40.1.tar.gz
# cd graphviz-2.40.1
# ./configure
# make
# couldnt get it to work oh well 


# dot pipline_workflow.dot  -T"$filetype" -o pipline_workflow."$filetype" && open pipline_workflow."$filetype"
# dot data_workflow.dot -Tpdf -o data_workflow.pdf && open data_workflow.pdf
# dot test.dot  -Tpng -o test.png
# dot placement.dot -Tpng -o placement3.png && open placement3.png

dot pipline_workflow.dot -Tpdf -o pipline_workflow.pdf && firefox pipline_workflow.pdf & 