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

# use this command on the server
# dot pipline_workflow.dot -Tpdf -o pipline_workflow.pdf && firefox pipline_workflow.pdf &

# use this command on local desktop
dot pipline_workflow.dot -Tpdf -o pipline_workflow.pdf # && open pipline_workflow.pdf &
dot pipline_workflow.dot -Tpng -o pipline_workflow.png && open pipline_workflow.png


dot data_workflow.dot -Tpdf -o data_workflow.pdf
dot data_workflow.dot -Tpng -o data_workflow.png # && open data_workflow.png &
