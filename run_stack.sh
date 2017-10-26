#!/bin/bash

# execute as ./run_stack.sh

  args=("$@")
  if [ $# -lt 2 ] ; then
    echo "Please provide run era (BCD,EF,FG,H) and Recorded Luminosity"
  fi

  if [ $# -eq 2 ] ; then

    run_era=${args[0]}
    luminosity=${args[1]}

  fi
  mc="root_files_R48/SingleNeutrino_MC_R4.root"
  data="root_files_R48/Legacy_${run_era}_R4.root"

  var="nPU"
  ratio="true"
  label="Run 2016${run_era} - ${luminosity} fb^{-1} (13 TeV)"

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
