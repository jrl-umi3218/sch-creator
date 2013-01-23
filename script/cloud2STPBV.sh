#! /bin/sh

local_dir=`dirname $0`

if [ $# -lt 2 ]
then
  echo "Usage: `basename $0` input_path output_path [R=300, r=0.2]"
  exit 1
fi

in_path=$1
out_path=$2
R=300
r=0.2

[ $3 ] && R=$3
[ $4 ] && r=$4

if [ ! -d "$in_path" ]
then
  echo "'$in_path' is not a valid directory."
  exit 2
fi

echo "$out_path"
if [ ! -d "$out_path" ]
then
  echo "'$out_path' is not a valid directory."
  exit 2
fi

for file in $in_path/*.qc
do
  echo $file
  file_name=`basename $file .qc`

  STPBV_generator -r $r -R $R $file "$out_path/$file_name.txt"
done

