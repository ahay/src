#!/bin/bash

    # scontrol is a slurm command. in the next list it produces a list
    # of the individual nodes eg c557-[204,301] changes to:
    # c557-204
    # c557-301
    declare -a hostlist=(`scontrol show hostname $SLURM_NODELIST `)
    declare -a node_clusters=(`echo $SLURM_TASKS_PER_NODE | sed -e's/,/ /g'`)

    if [ $? -ne 0  ];then
        echo "TACC: Error -> slurm host list unavailable"
        exit 1
    fi
    echo "DEBUG: ${hostlist[@]} "

    #Initialize the hostlist index
    host_id=0
    #Set up the PE task index to store wayness for each task
    # this is for tacc_affinity
    task_id=0

    export output=""
    #Build the hostlist using the SLURM_TASKS_PER_NODE syntax
    for nodes in ${node_clusters[@]}; do
      echo nodes=$nodes
      #Get the task count and node count for each node cluster
      task_count=`echo $nodes | awk -F '(' '{print $1}'`
      if [[ `echo $nodes | grep x` ]]; then
        node_count=`echo $nodes | sed -e's/.*x\([0-9]\+\).*/\1/'`
      else
        node_count=1
      fi
     echo ${hostlist[${host_id}]} $task_count
     echo $host_id
      echo "DEBUG: nodes=$nodes task_count=$task_count  node_count=$node_count"

      #Build the host list to match tasks per node
      # and set up PE env variable for each task
      echo node_count=$((node_count-1))
      for i in `seq 0 $((node_count-1))`; do
         output="${output} ${hostlist[${host_id}]} $task_count" 
      ((host_id++))
      done
    done

echo $output
