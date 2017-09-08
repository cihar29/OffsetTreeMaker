#!/bin/bash

# execute as ./run_stack.sh

  mc="MC_R4.root"
  data="Data_R4.root"

  var="nPU"
  ratio="true"
  label="Run 2016 - 35.9 fb^{-1} (13 TeV)"

  cmds=( "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "all", $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "ne",  $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "hfe", $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "nh",  $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "hfh", $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "chu", $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "chm", $ratio, \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"$mc\", \"$data\", \"$var\", "untrk", $ratio, \"$label\")'"
       )

  for cmd in "${cmds[@]}"
  do
    eval $cmd
  done
